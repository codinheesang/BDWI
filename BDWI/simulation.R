# ============================================================
# Simulation study for ROI-specific d-NN model — FIXED worker version
#   - X and ROI are fixed at the observed data
#   - y is generated from chosen true beta and sigma (= rho in code)
#   - model is refitted for each simulated data set
#   - bias, RMSE, 95% HPD interval width, and coverage are summarized
# ============================================================
library(Rcpp)
library(coda)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(coda)

Rcpp::sourceCpp("fnirs_bdwi.cpp")

input_file <- "fnirs_roi8_old2_igt_resting_basis19.RData"
load(input_file)

set.seed(2)
ind<-sample(c(1:nrow(datt)), nrow(datt)*0.2)

K <- ncol(datt) - 2   

train<-datt[ind,]
test<-datt[-ind,]

pp<-dim(train)[2]
kk<-length(unique(train[,(K+1)]))
basiss <- ifelse(grepl("basis19", input_file), 21, 101)
J <- length(unique(train[,basiss]))

dnn_mcmc_roi_theta <- function(
    training_data,
    rho,
    J,
    loop,
    burn,
    loop_aux,
    type = c("norm", "exp"),
    mu_alpha = -1,
    s2_alpha = 0.2,
    mu_rho   = 0,
    s2_rho   = 0.25,
    a_tau = 8,
    b_tau = 2,
    theta_init = rep(mu_alpha, J),
    tau2_init = 1,
    sd_theta   = 0.05,
    sd_log_rho = 0.05,
    target_acc = 0.234,
    print_every = 200
){
  type <- match.arg(type)
  
  n  <- nrow(training_data)
  p2 <- ncol(training_data)
  
  x     <- training_data[, 1:(p2 - 2), drop = FALSE]
  group <- as.integer(training_data[, p2 - 1])
  roi   <- as.integer(training_data[, p2])
  
  k <- max(group)
  
  stopifnot(length(theta_init) == J)
  stopifnot(all(roi >= 1), all(roi <= J))
  stopifnot(rho > 0, tau2_init > 0)
  
  neigh_dist <- nearest_neighs_k_cpp(x, group)
  wt <- dist2wt_cpp(neigh_dist, rho, n, type)
  
  theta <- as.numeric(theta_init)
  tau2  <- tau2_init
  
  rho_chain   <- numeric(loop)
  tau2_chain  <- numeric(loop)
  theta_chain <- matrix(NA_real_, nrow = loop, ncol = J)
  
  # derived quantities for convenience
  alpha_chain <- numeric(loop)
  beta_chain  <- numeric(loop)
  u_chain     <- matrix(NA_real_, nrow = loop, ncol = J)
  beta_roi_chain <- matrix(NA_real_, nrow = loop, ncol = J)
  
  acc_rho   <- 0
  acc_theta <- integer(J)
  
  for(l in 1:loop){
    upd_rho <- exchange_move_rho_theta_cpp_adapt(
      k = k,
      rho = rho,
      theta = theta,
      neigh_dist = neigh_dist,
      group = group,
      roi = roi,
      n = n,
      loop_aux = loop_aux,
      accept = acc_rho,
      wt = wt,
      type = type,
      sd_log_rho = sd_log_rho,
      mu_rho = mu_rho,
      s2_rho = s2_rho,
      iter = l,
      burnin = burn,
      target_acc = target_acc
    )
    
    rho <- upd_rho$rho
    wt  <- upd_rho$wt
    acc_rho <- upd_rho$accept
    sd_log_rho <- upd_rho$sd_log_rho
    for(j in 1:J){
      upd_theta <- exchange_move_theta_single_cpp_adapt(
        k = k,
        theta = theta,
        j_update_1based = j,
        tau2 = tau2,
        mu_alpha = mu_alpha,
        s2_alpha = s2_alpha,
        group = group,
        roi = roi,
        n = n,
        loop_aux = loop_aux,
        accept = acc_theta[j],
        wt = wt,
        sd_theta = sd_theta,
        iter = l,
        burnin = burn,
        target_acc = target_acc
      )
      theta <- upd_theta$theta
      acc_theta[j] <- upd_theta$accept
      sd_theta <- upd_theta$sd_theta
    }
    alpha_now <- mean(theta)
    u_now <- theta - alpha_now
    
    shape <- a_tau + (J - 1) / 2
    rate  <- b_tau + 0.5 * sum(u_now^2)
    tau2  <- 1 / rgamma(1, shape = shape, rate = rate)
    
    # store
    theta_chain[l, ] <- theta
    rho_chain[l]     <- rho
    tau2_chain[l]    <- tau2
    
    alpha_chain[l] <- alpha_now
    beta_chain[l]  <- (exp(alpha_now))
    u_chain[l, ]   <- u_now
    beta_roi_chain[l, ] <- (exp(theta))
    
    if(l %% print_every == 0){
      cat(sprintf(
        "loop=%d beta=%.4f rho=%.4f tau2=%.4f beta_roi_mean=%.4f acc_rho=%.3f acc_theta_mean=%.3f\n",
        l,
        (exp(alpha_now)),
        rho,
        tau2,
        mean((exp(theta))),
        acc_rho / l,
        mean(acc_theta / l)
      ))
    }
  }
  
  list(
    theta = theta_chain,
    rho   = rho_chain,
    tau2  = tau2_chain,
    alpha = alpha_chain,
    beta  = beta_chain,
    u     = u_chain,
    beta_roi = beta_roi_chain,
    accept_rho = acc_rho / loop,
    accept_theta = acc_theta / loop
  )
}

summarize_posterior_roi_theta <- function(fit, burn = 1000, prob = 0.95){
  
  idx <- seq(burn + 1, length(fit$rho))
  beta_post     <- fit$beta[idx]
  rho_post      <- fit$rho[idx]
  tau2_post     <- fit$tau2[idx]
  alpha_post    <- fit$alpha[idx]
  theta_post    <- fit$theta[idx, , drop = FALSE]
  u_post        <- fit$u[idx, , drop = FALSE]
  beta_roi_post <- fit$beta_roi[idx, , drop = FALSE]
  
  beta_mean <- mean(beta_post)
  beta_hpd  <- coda::HPDinterval(coda::as.mcmc(beta_post), prob = prob)
  
  rho_mean <- mean(rho_post)
  rho_hpd  <- coda::HPDinterval(coda::as.mcmc(rho_post), prob = prob)
  
  tau2_mean <- mean(tau2_post)
  tau2_hpd  <- coda::HPDinterval(coda::as.mcmc(tau2_post), prob = prob)
  
  alpha_mean <- mean(alpha_post)
  alpha_hpd  <- coda::HPDinterval(coda::as.mcmc(alpha_post), prob = prob)
  
  theta_mean <- colMeans(theta_post)
  theta_hpd <- t(apply(theta_post, 2, function(x)
    coda::HPDinterval(coda::as.mcmc(x), prob = prob)
  ))
  colnames(theta_hpd) <- c("lower", "upper")
  
  u_mean <- colMeans(u_post)
  u_hpd <- t(apply(u_post, 2, function(x)
    coda::HPDinterval(coda::as.mcmc(x), prob = prob)
  ))
  colnames(u_hpd) <- c("lower", "upper")
  
  beta_roi_mean <- colMeans(beta_roi_post)
  beta_roi_hpd <- t(apply(beta_roi_post, 2, function(x)
    coda::HPDinterval(coda::as.mcmc(x), prob = prob)
  ))
  colnames(beta_roi_hpd) <- c("lower", "upper")
  
  list(
    beta = list(mean = beta_mean, hpd = beta_hpd),
    rho  = list(mean = rho_mean, hpd = rho_hpd),
    tau2 = list(mean = tau2_mean, hpd = tau2_hpd),
    alpha = list(mean = alpha_mean, hpd = alpha_hpd),
    theta = list(mean = theta_mean, hpd = theta_hpd),
    u = list(mean = u_mean, hpd = u_hpd),
    beta_roi = list(mean = beta_roi_mean, hpd = beta_roi_hpd)
  )
}

simulate_y_fixed_x <- function(
    x,
    roi,
    beta_true,
    rho_true,            # sigma in the paper, rho in the current code
    k = 2L,
    type = c("norm", "exp"),
    sweeps_gen = 1000,
    max_attempts = 50,
    seed = NULL,
    D_fixed = NULL,
    wt_true_fixed = NULL
) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)
  roi <- as.integer(roi)
  beta_true <- as.numeric(beta_true)

  n <- nrow(x)
  J <- length(beta_true)

  stopifnot(length(roi) == n)
  stopifnot(all(roi >= 1L & roi <= J))
  stopifnot(all(beta_true > 0))
  stopifnot(rho_true > 0)
  stopifnot(k >= 2L)

  # nearest_neighs_k_cpp() currently does not use group in the distance
  # calculation, but a dummy vector is supplied to match its interface.
  if (is.null(D_fixed)) {
    D <- nearest_neighs_k_cpp(
      x = x,
      group = rep(1L, n)
    )
  } else {
    D <- D_fixed
  }

  if (is.null(wt_true_fixed)) {
    wt_true <- dist2wt_cpp(
      D = D,
      rho = rho_true,
      n = n,
      type = type
    )
  } else {
    wt_true <- wt_true_fixed
  }

  for (attempt in seq_len(max_attempts)) {
    # Random starting configuration. After enough sweeps this is only an
    # initialization for the Markov chain used to draw y.
    z0 <- sample.int(k, n, replace = TRUE)

    y_sim <- gibbs_draw_k1_roi_shared_cpp(
      k = k,
      n = n,
      beta_roi = beta_true,
      sweeps = sweeps_gen,
      wt = wt_true,
      roi = roi,
      z0 = z0
    )

    # Current fitting code assumes a genuine k-class problem.
    if (length(unique(y_sim)) == k) {
      return(list(
        y = as.integer(y_sim),
        D = D,
        wt = wt_true,
        attempts = attempt
      ))
    }
  }

  stop(
    "Could not generate all k classes after ", max_attempts,
    " attempts. The chosen beta/rho may be causing near-degeneracy."
  )
}

run_one_simulation <- function(
    sim_id,
    x,
    roi,
    beta_true,
    rho_true,
    # generation settings
    k = 2L,
    type = "norm",
    sweeps_gen = 1000,
    # MCMC settings: matched to the current application code by default
    loop = 20000,
    burn = 10000,
    loop_aux = 10,
    mu_alpha = -1,
    s2_alpha = 0.2,
    mu_rho = 0.5,
    s2_rho = 0.2,
    a_tau = 2,
    b_tau = 1,
    theta_init = NULL,
    tau2_init = 0.15,
    sd_theta = 0.1,
    sd_log_rho = 0.1,
    target_acc = 0.20,
    prob = 0.95,
    seed_base = 20260912,
    D_fixed = NULL,
    wt_true_fixed = NULL
) {
  J <- length(beta_true)
  theta_true <- log(beta_true)
  total_start <- proc.time()[["elapsed"]]

  if (is.null(theta_init)) theta_init <- rep(0, J)

  # Separate reproducible seed for each replication
  set.seed(seed_base + sim_id)

  # ---- Generate y ----
  gen_start <- proc.time()[["elapsed"]]
  gen <- simulate_y_fixed_x(
    x = x,
    roi = roi,
    beta_true = beta_true,
    rho_true = rho_true,
    k = k,
    type = type,
    sweeps_gen = sweeps_gen,
    D_fixed = D_fixed,
    wt_true_fixed = wt_true_fixed
  )
  generation_time_sec <- proc.time()[["elapsed"]] - gen_start

  sim_data <- cbind(
    x,
    label = gen$y,
    roi = roi
  )

  # ---- Refit the same model ----
  fit_start <- proc.time()[["elapsed"]]
  fit <- dnn_mcmc_roi_theta(
    training_data = sim_data,
    rho = rho_true,      # starting value only; rho is still estimated
    J = J,
    loop = loop,
    burn = burn,
    loop_aux = loop_aux,
    type = type,
    mu_alpha = mu_alpha,
    s2_alpha = s2_alpha,
    mu_rho = mu_rho,
    s2_rho = s2_rho,
    a_tau = a_tau,
    b_tau = b_tau,
    theta_init = theta_init,
    tau2_init = tau2_init,
    sd_theta = sd_theta,
    sd_log_rho = sd_log_rho,
    target_acc = target_acc,
    print_every = loop + 1L   # suppress usual iteration printing
  )

  post <- summarize_posterior_roi_theta(
    fit = fit,
    burn = burn,
    prob = prob
  )
  fitting_time_sec <- proc.time()[["elapsed"]] - fit_start
  total_time_sec <- proc.time()[["elapsed"]] - total_start

  beta_hat <- as.numeric(post$beta_roi$mean)
  beta_ci <- as.matrix(post$beta_roi$hpd)
  rho_hat <- as.numeric(post$rho$mean)
  rho_ci <- as.numeric(post$rho$hpd)

  # theta summaries are also useful because beta = exp(theta)
  theta_hat <- as.numeric(post$theta$mean)
  theta_ci <- as.matrix(post$theta$hpd)

  beta_cover <- beta_ci[, 1] <= beta_true & beta_true <= beta_ci[, 2]
  rho_cover <- rho_ci[1] <= rho_true & rho_true <= rho_ci[2]
  theta_cover <- theta_ci[, 1] <= theta_true & theta_true <= theta_ci[, 2]

  list(
    sim_id = sim_id,
    result = data.frame(
      sim = sim_id,
      parameter = c(paste0("beta", seq_len(J)), "rho"),
      true = c(beta_true, rho_true),
      estimate = c(beta_hat, rho_hat),
      lower = c(beta_ci[, 1], rho_ci[1]),
      upper = c(beta_ci[, 2], rho_ci[2]),
      cover = c(beta_cover, rho_cover),
      stringsAsFactors = FALSE
    ),
    theta_result = data.frame(
      sim = sim_id,
      parameter = paste0("theta", seq_len(J)),
      true = theta_true,
      estimate = theta_hat,
      lower = theta_ci[, 1],
      upper = theta_ci[, 2],
      cover = theta_cover,
      stringsAsFactors = FALSE
    ),
    diagnostics = data.frame(
      sim = sim_id,
      class1_prop = mean(gen$y == 1L),
      class2_prop = mean(gen$y == 2L),
      generation_attempts = gen$attempts,
      accept_rho = fit$accept_rho,
      accept_theta_mean = mean(fit$accept_theta),
      generation_time_sec = generation_time_sec,
      fitting_time_sec = fitting_time_sec,
      total_time_sec = total_time_sec,
      stringsAsFactors = FALSE
    ),
    timing = data.frame(
      sim = sim_id,
      generation_time_sec = generation_time_sec,
      fitting_time_sec = fitting_time_sec,
      total_time_sec = total_time_sec,
      generation_time_min = generation_time_sec / 60,
      fitting_time_min = fitting_time_sec / 60,
      total_time_min = total_time_sec / 60,
      stringsAsFactors = FALSE
    )
  )
}

run_simulation_study <- function(
    datt,
    beta_true,
    rho_true,
    nsim = 100,
    type = "norm",
    sweeps_gen = 1000,
    loop = 20000,
    burn = 10000,
    loop_aux = 10,
    mu_alpha = -1,
    s2_alpha = 0.2,
    mu_rho = 0.5,
    s2_rho = 0.2,
    a_tau = 2,
    b_tau = 1,
    tau2_init = 0.15,
    sd_theta = 0.1,
    sd_log_rho = 0.1,
    target_acc = 0.20,
    prob = 0.95,
    seed_base = 20260912,
    checkpoint_file = "simulation_checkpoint.rds"
) {
  datt <- as.matrix(datt)

  p2 <- ncol(datt)
  x <- as.matrix(datt[, 1:(p2 - 2), drop = FALSE])
  roi <- as.integer(datt[, p2])

  J <- length(beta_true)
  if (length(unique(roi)) != J) {
    stop("length(beta_true) must match the number of ROI levels in datt.")
  }

  # X and true rho are fixed, so generation D/W only need to be computed once.
  D_fixed_worker <- nearest_neighs_k_cpp(x, rep(1L, nrow(x)))
  wt_true_fixed_worker <- dist2wt_cpp(D_fixed_worker, rho_true, nrow(x), type)

  all_res <- vector("list", nsim)

  for (s in seq_len(nsim)) {
    cat("\n==============================\n")
    cat("Simulation", s, "of", nsim, "\n")
    cat("==============================\n")

    one <- tryCatch(
      run_one_simulation(
        sim_id = s,
        x = x,
        roi = roi,
        beta_true = beta_true,
        rho_true = rho_true,
        k = 2L,
        type = type,
        sweeps_gen = sweeps_gen,
        loop = loop,
        burn = burn,
        loop_aux = loop_aux,
        mu_alpha = mu_alpha,
        s2_alpha = s2_alpha,
        mu_rho = mu_rho,
        s2_rho = s2_rho,
        a_tau = a_tau,
        b_tau = b_tau,
        tau2_init = tau2_init,
        sd_theta = sd_theta,
        sd_log_rho = sd_log_rho,
        target_acc = target_acc,
        prob = prob,
        seed_base = seed_base,
        D_fixed = D_fixed_worker,
        wt_true_fixed = wt_true_fixed_worker
      ),
      error = function(e) {
        message("Simulation ", s, " failed: ", conditionMessage(e))
        list(sim_id = s, error = conditionMessage(e))
      }
    )

    all_res[[s]] <- one

    # checkpoint: useful because 100 x 20,000-iteration fits can be expensive
    saveRDS(all_res, checkpoint_file)
  }

  ok <- vapply(all_res, function(z) !is.null(z$result), logical(1))
  if (!any(ok)) stop("All simulation replications failed.")

  result_long <- do.call(rbind, lapply(all_res[ok], `[[`, "result"))
  theta_long <- do.call(rbind, lapply(all_res[ok], `[[`, "theta_result"))
  diagnostics <- do.call(rbind, lapply(all_res[ok], `[[`, "diagnostics"))
  timing <- do.call(rbind, lapply(all_res[ok], `[[`, "timing"))

  timing_summary <- data.frame(
    metric = c("generation", "fitting", "total"),
    mean_minutes = c(mean(timing$generation_time_min), mean(timing$fitting_time_min), mean(timing$total_time_min)),
    median_minutes = c(median(timing$generation_time_min), median(timing$fitting_time_min), median(timing$total_time_min)),
    min_minutes = c(min(timing$generation_time_min), min(timing$fitting_time_min), min(timing$total_time_min)),
    max_minutes = c(max(timing$generation_time_min), max(timing$fitting_time_min), max(timing$total_time_min)),
    stringsAsFactors = FALSE
  )

  # Frequentist operating characteristics across simulated data sets
  parameter_summary <- do.call(
    rbind,
    lapply(split(result_long, result_long$parameter), function(d) {
      data.frame(
        parameter = d$parameter[1],
        true = d$true[1],
        mean_estimate = mean(d$estimate),
        bias = mean(d$estimate - d$true),
        rmse = sqrt(mean((d$estimate - d$true)^2)),
        coverage_95 = mean(d$cover),
        mean_CI_width = mean(d$upper - d$lower),
        n_success = nrow(d),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(parameter_summary) <- NULL

  theta_summary <- do.call(
    rbind,
    lapply(split(theta_long, theta_long$parameter), function(d) {
      data.frame(
        parameter = d$parameter[1],
        true = d$true[1],
        mean_estimate = mean(d$estimate),
        bias = mean(d$estimate - d$true),
        rmse = sqrt(mean((d$estimate - d$true)^2)),
        coverage_95 = mean(d$cover),
        mean_CI_width = mean(d$upper - d$lower),
        n_success = nrow(d),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(theta_summary) <- NULL

  out <- list(
    settings = list(
      nsim = nsim,
      beta_true = beta_true,
      theta_true = log(beta_true),
      rho_true = rho_true,
      type = type,
      sweeps_gen = sweeps_gen,
      loop = loop,
      burn = burn,
      loop_aux = loop_aux,
      prob = prob,
      seed_base = seed_base
    ),
    parameter_summary = parameter_summary,
    theta_summary = theta_summary,
    results = result_long,
    theta_results = theta_long,
    diagnostics = diagnostics,
    timing = timing,
    timing_summary = timing_summary,
    raw = all_res
  )

  out
}



# Helper: safely read one replication file.
# Returns NULL if the file is missing/corrupted/incomplete.
read_completed_replication <- function(file) {
  if (!file.exists(file)) return(NULL)
  z <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(z)) return(NULL)
  if (is.null(z$result) || is.null(z$theta_result) || is.null(z$diagnostics)) {
    return(NULL)
  }
  z
}

# Helper: combine individually saved replication files.
collect_simulation_results <- function(
    result_dir,
    nsim,
    beta_true,
    rho_true,
    type = "norm",
    sweeps_gen = NA_integer_,
    loop = NA_integer_,
    burn = NA_integer_,
    loop_aux = NA_integer_,
    prob = 0.95,
    seed_base = NA_integer_
) {
  result_files <- file.path(
    result_dir,
    sprintf("sim_%03d.rds", seq_len(nsim))
  )

  all_res <- lapply(result_files, read_completed_replication)
  ok <- !vapply(all_res, is.null, logical(1))

  if (!any(ok)) {
    stop("No completed simulation replications were found in: ", result_dir)
  }

  result_long <- do.call(rbind, lapply(all_res[ok], `[[`, "result"))
  theta_long <- do.call(rbind, lapply(all_res[ok], `[[`, "theta_result"))
  diagnostics <- do.call(rbind, lapply(all_res[ok], `[[`, "diagnostics"))

  timing_list <- lapply(all_res[ok], function(z) {
    if (!is.null(z$timing)) return(z$timing)
    if (all(c("generation_time_sec", "fitting_time_sec", "total_time_sec") %in% names(z$diagnostics))) {
      data.frame(
        sim = z$diagnostics$sim,
        generation_time_sec = z$diagnostics$generation_time_sec,
        fitting_time_sec = z$diagnostics$fitting_time_sec,
        total_time_sec = z$diagnostics$total_time_sec,
        generation_time_min = z$diagnostics$generation_time_sec / 60,
        fitting_time_min = z$diagnostics$fitting_time_sec / 60,
        total_time_min = z$diagnostics$total_time_sec / 60
      )
    } else NULL
  })
  timing_list <- Filter(Negate(is.null), timing_list)
  timing <- if (length(timing_list)) do.call(rbind, timing_list) else data.frame()

  timing_summary <- if (nrow(timing)) {
    data.frame(
      metric = c("generation", "fitting", "total"),
      mean_minutes = c(mean(timing$generation_time_min), mean(timing$fitting_time_min), mean(timing$total_time_min)),
      median_minutes = c(median(timing$generation_time_min), median(timing$fitting_time_min), median(timing$total_time_min)),
      min_minutes = c(min(timing$generation_time_min), min(timing$fitting_time_min), min(timing$total_time_min)),
      max_minutes = c(max(timing$generation_time_min), max(timing$fitting_time_min), max(timing$total_time_min)),
      stringsAsFactors = FALSE
    )
  } else data.frame()

  parameter_summary <- do.call(
    rbind,
    lapply(split(result_long, result_long$parameter), function(d) {
      data.frame(
        parameter = d$parameter[1],
        true = d$true[1],
        mean_estimate = mean(d$estimate),
        bias = mean(d$estimate - d$true),
        rmse = sqrt(mean((d$estimate - d$true)^2)),
        coverage_95 = mean(d$cover),
        mean_CI_width = mean(d$upper - d$lower),
        n_success = nrow(d),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(parameter_summary) <- NULL

  theta_summary <- do.call(
    rbind,
    lapply(split(theta_long, theta_long$parameter), function(d) {
      data.frame(
        parameter = d$parameter[1],
        true = d$true[1],
        mean_estimate = mean(d$estimate),
        bias = mean(d$estimate - d$true),
        rmse = sqrt(mean((d$estimate - d$true)^2)),
        coverage_95 = mean(d$cover),
        mean_CI_width = mean(d$upper - d$lower),
        n_success = nrow(d),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(theta_summary) <- NULL

  list(
    settings = list(
      nsim = nsim,
      beta_true = beta_true,
      theta_true = log(beta_true),
      rho_true = rho_true,
      type = type,
      sweeps_gen = sweeps_gen,
      loop = loop,
      burn = burn,
      loop_aux = loop_aux,
      prob = prob,
      seed_base = seed_base,
      result_dir = result_dir
    ),
    completion = data.frame(
      sim = seq_len(nsim),
      completed = ok,
      file = result_files,
      stringsAsFactors = FALSE
    ),
    parameter_summary = parameter_summary,
    theta_summary = theta_summary,
    results = result_long,
    theta_results = theta_long,
    diagnostics = diagnostics,
    timing = timing,
    timing_summary = timing_summary,
    raw = all_res[ok]
  )
}


run_simulation_study_parallel <- function(
    datt,
    beta_true,
    rho_true,
    nsim = 100,
    n_workers = 4,
    result_dir = "sim_results",
    # C++ source needed by each fresh worker process on Windows.
    # model_r_file is retained only for backward compatibility and is NOT sourced.
    model_r_file = NULL,
    cpp_file = NULL,
    # generation settings
    type = "norm",
    sweeps_gen = 1000,
    # MCMC settings
    loop = 20000,
    burn = 10000,
    loop_aux = 10,
    mu_alpha = -1,
    s2_alpha = 0.2,
    mu_rho = 0.5,
    s2_rho = 0.2,
    a_tau = 2,
    b_tau = 1,
    tau2_init = 0.15,
    sd_theta = 0.1,
    sd_log_rho = 0.1,
    target_acc = 0.20,
    prob = 0.95,
    seed_base = 20260912,
    retry_failed = TRUE,
    save_error_files = TRUE
) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The base R 'parallel' package is required.")
  }

  datt <- as.matrix(datt)
  p2 <- ncol(datt)
  x <- as.matrix(datt[, 1:(p2 - 2), drop = FALSE])
  roi <- as.integer(datt[, p2])

  J <- length(beta_true)
  if (length(unique(roi)) != J) {
    stop("length(beta_true) must match the number of ROI levels in datt.")
  }

  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }

  # Which replications are already safely completed?
  result_files <- file.path(
    result_dir,
    sprintf("sim_%03d.rds", seq_len(nsim))
  )
  existing <- lapply(result_files, read_completed_replication)
  completed <- !vapply(existing, is.null, logical(1))
  pending_ids <- which(!completed)

  cat("Completed:", sum(completed), "of", nsim, "\n")
  cat("Pending  :", length(pending_ids), "\n")

  if (length(pending_ids) == 0L) {
    cat("All replications are already complete. Reconstructing summary...\n")
    return(collect_simulation_results(
      result_dir = result_dir,
      nsim = nsim,
      beta_true = beta_true,
      rho_true = rho_true,
      type = type,
      sweeps_gen = sweeps_gen,
      loop = loop,
      burn = burn,
      loop_aux = loop_aux,
      prob = prob,
      seed_base = seed_base
    ))
  }

  n_workers <- max(1L, min(as.integer(n_workers), length(pending_ids)))
  cat("Using", n_workers, "workers.\n")

  # PSOCK works on Windows, macOS, and Linux.
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Load packages / source code / compile C++ separately inside each worker.
  parallel::clusterEvalQ(cl, {
    library(Rcpp)
    library(coda)
    NULL
  })

  if (!is.null(model_r_file)) {
    message("NOTE: model_r_file is intentionally not sourced. The full analysis script contains top-level fitting code. Clean model functions are embedded in this simulation script.")
  }

  if (!is.null(cpp_file)) {
    cpp_file <- normalizePath(cpp_file, winslash = "/", mustWork = TRUE)
    parallel::clusterExport(cl, "cpp_file", envir = environment())
    parallel::clusterEvalQ(cl, Rcpp::sourceCpp(cpp_file))
  }

  # Export R functions defined in this script and common fixed objects once.
  parallel::clusterExport(
    cl,
    varlist = c(
      "dnn_mcmc_roi_theta",
      "summarize_posterior_roi_theta",
      "simulate_y_fixed_x",
      "run_one_simulation",
      "x", "roi", "beta_true", "rho_true",
      "type", "sweeps_gen",
      "loop", "burn", "loop_aux",
      "mu_alpha", "s2_alpha", "mu_rho", "s2_rho",
      "a_tau", "b_tau", "tau2_init",
      "sd_theta", "sd_log_rho", "target_acc",
      "prob", "seed_base", "result_dir",
      "save_error_files"
    ),
    envir = environment()
  )

  # Precompute the fixed distance matrix and generation weight matrix ONCE per worker.
  # They are reused for every simulation assigned to that worker.
  parallel::clusterEvalQ(cl, {
    D_fixed_worker <- nearest_neighs_k_cpp(x, rep(1L, nrow(x)))
    wt_true_fixed_worker <- dist2wt_cpp(
      D_fixed_worker, rho_true, nrow(x), type
    )
    NULL
  })

  worker_fun <- function(s) {
    outfile <- file.path(result_dir, sprintf("sim_%03d.rds", s))
    errfile <- file.path(result_dir, sprintf("sim_%03d_ERROR.rds", s))

    # A second check makes resume safe even if this function is called again.
    old <- tryCatch(readRDS(outfile), error = function(e) NULL)
    if (!is.null(old) &&
        !is.null(old$result) &&
        !is.null(old$theta_result) &&
        !is.null(old$diagnostics)) {
      return(list(sim = s, status = "skipped_existing"))
    }

    ans <- tryCatch({
      one <- run_one_simulation(
        sim_id = s,
        x = x,
        roi = roi,
        beta_true = beta_true,
        rho_true = rho_true,
        k = 2L,
        type = type,
        sweeps_gen = sweeps_gen,
        loop = loop,
        burn = burn,
        loop_aux = loop_aux,
        mu_alpha = mu_alpha,
        s2_alpha = s2_alpha,
        mu_rho = mu_rho,
        s2_rho = s2_rho,
        a_tau = a_tau,
        b_tau = b_tau,
        tau2_init = tau2_init,
        sd_theta = sd_theta,
        sd_log_rho = sd_log_rho,
        target_acc = target_acc,
        prob = prob,
        seed_base = seed_base,
        D_fixed = D_fixed_worker,
        wt_true_fixed = wt_true_fixed_worker
      )

      # Write to a temporary file first, then rename.
      # This reduces the chance that a crash leaves a half-written RDS.
      tmpfile <- paste0(outfile, ".tmp_", Sys.getpid())
      saveRDS(one, tmpfile)

      ok <- file.rename(tmpfile, outfile)
      if (!ok) {
        # file.rename may fail across unusual filesystems; fallback to copy.
        ok2 <- file.copy(tmpfile, outfile, overwrite = TRUE)
        unlink(tmpfile)
        if (!ok2) stop("Could not finalize result file: ", outfile)
      }

      if (file.exists(errfile)) unlink(errfile)
      list(sim = s, status = "completed")
    }, error = function(e) {
      if (isTRUE(save_error_files)) {
        saveRDS(
          list(
            sim = s,
            error = conditionMessage(e),
            time = Sys.time()
          ),
          errfile
        )
      }
      list(sim = s, status = "failed", error = conditionMessage(e))
    })

    ans
  }

  # Dynamic scheduling is preferable because some MCMC fits may take longer.
  status <- parallel::parLapplyLB(cl, pending_ids, worker_fun)

  status_df <- do.call(
    rbind,
    lapply(status, function(z) {
      data.frame(
        sim = z$sim,
        status = z$status,
        error = if (!is.null(z$error)) z$error else NA_character_,
        stringsAsFactors = FALSE
      )
    })
  )

  print(status_df)

  # Optional one-time retry of replications that raised ordinary R errors.
  # Catastrophic worker/R-session termination is naturally handled by rerunning
  # this whole function; completed RDS files will be skipped.
  if (isTRUE(retry_failed) && any(status_df$status == "failed")) {
    failed_ids <- status_df$sim[status_df$status == "failed"]
    cat("Retrying", length(failed_ids), "failed replication(s)...\n")
    retry_status <- parallel::parLapplyLB(cl, failed_ids, worker_fun)
    retry_df <- do.call(
      rbind,
      lapply(retry_status, function(z) {
        data.frame(
          sim = z$sim,
          status = paste0("retry_", z$status),
          error = if (!is.null(z$error)) z$error else NA_character_,
          stringsAsFactors = FALSE
        )
      })
    )
    status_df <- rbind(status_df, retry_df)
    print(retry_df)
  }

  # Rebuild the complete output only from valid per-replication files.
  out <- collect_simulation_results(
    result_dir = result_dir,
    nsim = nsim,
    beta_true = beta_true,
    rho_true = rho_true,
    type = type,
    sweeps_gen = sweeps_gen,
    loop = loop,
    burn = burn,
    loop_aux = loop_aux,
    prob = prob,
    seed_base = seed_base
  )

  out$run_status <- status_df

  # Also keep a lightweight summary/checkpoint in the result directory.
  saveRDS(out, file.path(result_dir, "simulation_summary_latest.rds"))
  write.csv(
    out$parameter_summary,
    file.path(result_dir, "parameter_summary_latest.csv"),
    row.names = FALSE
  )
  write.csv(
    out$completion,
    file.path(result_dir, "completion_status.csv"),
    row.names = FALSE
  )
  if (nrow(out$timing)) {
    write.csv(
      out$timing,
      file.path(result_dir, "timing_by_replication.csv"),
      row.names = FALSE
    )
    write.csv(
      out$timing_summary,
      file.path(result_dir, "timing_summary.csv"),
      row.names = FALSE
    )
  }

  out
}


# IMPORTANT FOR WINDOWS:
# Only cpp_file is needed by fresh PSOCK workers.
# top-level data loading and a 20,000-iteration model fit.
#
beta_true <- c(4.3, 3.0, 3.6, 3.4, 4.1, 3.1, 3.3, 3.6)
rho_true  <- 0.2

sim_res <- run_simulation_study_parallel(
   datt = train,
   beta_true = beta_true,
   rho_true = rho_true,
   nsim = 100,
   n_workers = 8,            # start conservatively for RAM safety
   result_dir = "sim_results2",
   cpp_file = "fnirs_bdwi.cpp",
   type = "norm",
   sweeps_gen = 100,
   loop = 20000,
   burn = 10000,
   loop_aux = 20,
   mu_alpha = -1,
   s2_alpha = 0.2,
   mu_rho = 0.5,
   s2_rho = 0.2,
   a_tau = 2,
   b_tau = 1,
   tau2_init = 0.15,
   sd_theta = 0.1,
   sd_log_rho = 0.1,
   target_acc = 0.20,
   prob = 0.95,
   seed_base = 20260912
)
#
# # If R/PC crashes, simply run EXACTLY the same call again.
# # Existing sim_001.rds, sim_002.rds, ... are detected and skipped.
save(sim_res, file="old_igt_final_simu_loop20.RData")


