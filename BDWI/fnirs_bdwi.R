rm(list = ls())

library(Rcpp)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(coda)
Rcpp::sourceCpp("fnirs_bdwi.cpp")

input_file <- "fnirs_roi8_young_igt_resting_basis19.RData"
input_file <- "fnirs_roi8_young_stroop_resting_basis19.RData"
input_file <- "fnirs_roi8_old2_igt_resting_basis19.RData"
input_file <- "fnirs_roi8_old2_stroop_resting_basis19.RData"

##choose one file
load(input_file)


#############################################
cum_p2map_p <- function(cum_p, k) {
  size_p <- nrow(cum_p)
  c1 <- matrix(0, size_p, k)
  for (i in 1:size_p) {
    row <- cum_p[i, ]
    c1[i, ] <- as.integer(row == max(row))
    if (sum(c1[i, ]) == k) {
      c1[i, ] <- 0
      cen <- ceiling(k * runif(1))
      c1[i, cen] <- 1
    }
  }
  c_tmp <- matrix(rep(1:k, each = size_p), ncol = k)
  c_tmp <- c_tmp * c1
  c <- rowSums(c_tmp)
  c
}


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
    # ---------------------------------------------------------
    # (1) exchange update for rho | theta
    # ---------------------------------------------------------
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
    # ---------------------------------------------------------
    # (2) one-at-a-time update for theta_j
    # ---------------------------------------------------------
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
    # ---------------------------------------------------------
    # (3) Gibbs update for tau2 | theta
    # u = theta - mean(theta), dimension J-1 effectively
    # ---------------------------------------------------------
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


dnn_predict_mode_roi_theta <- function(
    fit,
    training_data,
    test_data,
    burn,
    type = c("norm", "exp"),
    thin = 1,
    sweep = 100,
    return_prob = TRUE
){
  type <- match.arg(type)
  n_train <- nrow(training_data)
  p2_tr   <- ncol(training_data)
  x     <- training_data[, 1:(p2_tr - 2), drop = FALSE]
  group <- as.integer(training_data[, p2_tr - 1])
  p2_te  <- ncol(test_data)
  y      <- test_data[, 1:(p2_te - 2), drop = FALSE]
  test_y <- as.integer(test_data[, p2_te - 1])
  roi_y  <- as.integer(test_data[, p2_te])
  k <- max(group)
  neigh_disty <- nearest_neighsy_k_cpp(y, x, group)
  idx <- seq(burn + 1, length(fit$rho), by = thin)
  T_eff <- length(idx)
  m <- nrow(y)
  label_chain <- matrix(NA_integer_, nrow = m, ncol = T_eff)
  prob_acc <- if(return_prob) matrix(0, nrow = m, ncol = k) else NULL
  iter_id <- 1
  for(t in idx){
    wty <- dist2wt_y_cpp(
      neigh_disty,
      rho = fit$rho[t],
      n   = n_train,
      type = type
    )
    beta_roi_y <- (exp(fit$theta[t, roi_y]))
    z <- sample.int(k, m, replace = TRUE)
    for(sw in 1:sweep){
      for(i in 1:m){
        zi <- z[i]
        prop <- sample.int(k, 1)
        while(prop == zi) prop <- sample.int(k, 1)
        cur <- sum(wty[i, group == zi])
        pro <- sum(wty[i, group == prop])
        log_alpha <- beta_roi_y[i] * (pro - cur)
        if(log(runif(1)) <= log_alpha) z[i] <- prop
      }
    }
    label_chain[, iter_id] <- z
    if(return_prob){
      for(c in 1:k) prob_acc[, c] <- prob_acc[, c] + as.integer(z == c)
    }
    iter_id <- iter_id + 1
  }
  c_hat_mode <- apply(label_chain, 1, function(v){
    tab <- tabulate(v, nbins = k)
    which.max(tab)
  })
  ok <- !is.na(test_y)
  misclass <- if(any(ok)) 1 - mean(test_y[ok] == c_hat_mode[ok]) else NA_real_
  out <- list(
    label_chain = label_chain,
    class = c_hat_mode,
    misclass = misclass
  )
  if(return_prob){
    prob <- prob_acc / T_eff
    colnames(prob) <- paste0("class", 1:k)
    out$prob <- prob
    out$uncertainty <- 1 - apply(prob, 1, max)
  }
  out
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


diagnose_mcmc_roi_theta <- function(
    fit,
    burn = 1000,
    roi_show = 1:min(4, ncol(fit$theta))
){
  idx <- seq(burn + 1, length(fit$rho))
  beta_post  <- fit$beta[idx]
  rho_post   <- fit$rho[idx]
  tau2_post  <- fit$tau2[idx]
  theta_post <- fit$theta[idx, , drop = FALSE]
  u_post     <- fit$u[idx, , drop = FALSE]
  ess_beta  <- coda::effectiveSize(beta_post)
  ess_rho   <- coda::effectiveSize(rho_post)
  ess_tau2  <- coda::effectiveSize(tau2_post)
  ess_theta <- apply(theta_post, 2, coda::effectiveSize)
  ess_u     <- apply(u_post, 2, coda::effectiveSize)
  cat("ESS beta :", ess_beta, "\n")
  cat("ESS rho  :", ess_rho, "\n")
  cat("ESS tau2 :", ess_tau2, "\n")
  cat("ESS theta min/median/max:",
      min(ess_theta), median(ess_theta), max(ess_theta), "\n")
  cat("ESS u min/median/max:",
      min(ess_u), median(ess_u), max(ess_u), "\n")
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(2, 3))
  plot(beta_post, type = "l", main = "beta trace", xlab = "iter", ylab = "beta")
  plot(rho_post,  type = "l", main = "rho trace",  xlab = "iter", ylab = "rho")
  plot(tau2_post, type = "l", main = "tau2 trace", xlab = "iter", ylab = "tau2")
  for(j in roi_show[1:min(3, length(roi_show))]){
    plot(theta_post[, j], type = "l",
         main = paste0("theta_", j, " trace"),
         xlab = "iter", ylab = "theta")
  }
  invisible(list(
    ess_beta = ess_beta,
    ess_rho  = ess_rho,
    ess_tau2 = ess_tau2,
    ess_theta = ess_theta,
    ess_u = ess_u
  ))
}


set.seed(2)
ind<-sample(c(1:nrow(datt)), nrow(datt)*0.7)

K <- ncol(datt) - 2   # number of feature

train<-datt[ind,]
test<-datt[-ind,]

pp<-dim(train)[2]
kk<-length(unique(train[,(K+1)]))
basiss <- ifelse(grepl("basis19", input_file), 21, 101)
J <- length(unique(train[,basiss]))

start <- proc.time()
fit <- dnn_mcmc_roi_theta(
  training_data = train,  # [features..., label, roi]
  rho  = 1,
  J = J,
  loop = 20000,
  burn = 10000,
  loop_aux = 10,
  type = "norm",
  mu_alpha = -1,
  s2_alpha = 0.2,
  mu_rho = 0.5,
  s2_rho = 0.2,
  a_tau = 2,
  b_tau = 1,
  theta_init = rep(0, J),
  tau2_init = 0.15,
  sd_theta = 0.1,
  sd_log_rho  = 0.1,
  target_acc  = 0.20,
  print_every = 1000
)
mcmc_time<-proc.time() - start


start <- proc.time()
pred <- dnn_predict_mode_roi_theta(
  fit = fit,
  training_data = train,
  test_data = test,
  burn = 10000,
  type = "norm",
  thin = 10,
  sweep = 100,
  return_prob = TRUE
)
pred_time<-proc.time() - start

pred$misclass
mcmc_time
pred_time

posterior_summary <- summarize_posterior_roi_theta(fit, burn = 10000)
posterior_summary$rho
posterior_summary$beta_roi

diag_res <- diagnose_mcmc_roi_theta(fit, burn = 10000, roi_show = 1:8)

output_file <- sub("\\.RData$", "_result.RData", input_file)
save(fit, mcmc_time, pred, pred_time, posterior_summary, diag_res, file=output_file)


#####posterior predictive checks
# ------------------------------------------------------------
# 0. label re-index helpers
# ------------------------------------------------------------

make_label_index <- function(group_raw) {
  label_set <- sort(unique(as.integer(group_raw)))
  group_idx <- match(as.integer(group_raw), label_set)
  list(
    group_idx = as.integer(group_idx),   # 1,2,...,K
    label_set = as.integer(label_set)    # original labels
  )
}
recover_original_labels <- function(group_idx, label_set) {
  as.integer(label_set[as.integer(group_idx)])
}
# ------------------------------------------------------------
# 1. statistics
# ------------------------------------------------------------
weighted_agreement_stat <- function(group, wt) {
  group <- as.integer(group)
  same_mat <- outer(group, group, "==") * 1
  sum(wt * same_mat) / sum(wt)
}
class_prop_stat <- function(group, label_set = sort(unique(group))) {
  group <- as.integer(group)
  label_set <- as.integer(label_set)
  out <- numeric(length(label_set))
  names(out) <- as.character(label_set)
  for (i in seq_along(label_set)) {
    out[i] <- mean(group == label_set[i])
  }
  out
}

roi_agreement_stat <- function(group, wt, roi, J = max(roi)) {
  group <- as.integer(group)
  roi   <- as.integer(roi)
  same_mat <- outer(group, group, "==") * 1
  out <- rep(NA_real_, J)
  for (j in 1:J) {
    idx <- which(roi == j)
    if (length(idx) == 0) {
      out[j] <- NA_real_
    } else {
      num <- sum(wt[idx, , drop = FALSE] * same_mat[idx, , drop = FALSE])
      den <- sum(wt[idx, , drop = FALSE])
      out[j] <- ifelse(den > 0, num / den, NA_real_)
    }
  }
  out
}

misclassification_rate_stat <- function(group_obs, group_rep) {
  mean(as.integer(group_obs) != as.integer(group_rep))
}

simulate_rep_group_dnn_roi_theta <- function(theta, rho, x, roi, group_raw,
                                             sweeps = 300,
                                             type = c("norm", "exp")) {
  type <- match.arg(type)
  x <- as.matrix(x)
  roi <- as.integer(roi)
  group_raw <- as.integer(group_raw)
  lab <- make_label_index(group_raw)
  group_idx <- lab$group_idx
  label_set <- lab$label_set
  n <- nrow(x)
  k <- length(label_set)
  # distance / weight
  D <- nearest_neighs_k_cpp(x, group_idx)
  wt <- dist2wt_cpp(D, rho, n, type)
  beta_roi <- make_beta_roi_from_theta(theta)
  # replicated labels in indexed space {1,...,k}
  z_rep_idx <- gibbs_draw_k1_roi_shared_cpp(
    k = k,
    n = n,
    beta_roi = beta_roi,
    sweeps = sweeps,
    wt = wt,
    roi = roi,
    z0 = group_idx
  )
  # recover original labels, e.g. {1,3}
  z_rep_raw <- recover_original_labels(z_rep_idx, label_set)
  list(
    z_rep = as.integer(z_rep_raw),     # original label space
    z_rep_idx = as.integer(z_rep_idx), # contiguous label space
    group_idx = group_idx,
    label_set = label_set,
    wt = wt,
    D = D,
    beta_roi = beta_roi
  )
}


roi_energy_stat <- function(
    group,
    wt,
    beta_roi,
    roi,
    J = max(roi)
) {
  group <- as.integer(group)
  roi <- as.integer(roi)
  beta_roi <- as.numeric(beta_roi)
  same_mat <- outer(group, group, "==") * 1
  out <- rep(NA_real_, J)
  names(out) <- paste0("ROI", 1:J)
  for (j in 1:J) {
    idx <- which(roi == j)
    if (length(idx) > 0) {
      out[j] <- beta_roi[j] *
        sum(
          wt[idx, , drop = FALSE] *
            same_mat[idx, , drop = FALSE]
        )
    }
  }
  out
}


ppc_dnn_mcmc_roi_theta <- function(fit,
                                   training_data,
                                   burn = 0,
                                   thin = 1,
                                   n_draws = 200,
                                   sweeps = 300,
                                   type = c("norm", "exp"),
                                   seed = NULL,
                                   verbose = TRUE) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  n_total <- nrow(fit$theta)
  if (burn >= n_total) stop("burn must be smaller than nrow(fit$theta)")
  p2 <- ncol(training_data)
  x <- as.matrix(training_data[, 1:(p2 - 2), drop = FALSE])
  group_obs <- as.integer(training_data[, p2 - 1])
  roi <- as.integer(training_data[, p2])
  J <- max(roi)
  label_info <- make_label_index(group_obs)
  label_set <- label_info$label_set
  k <- length(label_set)
  keep_all <- seq(from = burn + 1, to = n_total, by = thin)
  if (length(keep_all) == 0) stop("No posterior draws left after burn/thin.")
  if (length(keep_all) > n_draws) {
    keep <- sort(sample(keep_all, n_draws))
  } else {
    keep <- keep_all
  }
  M <- length(keep)
  rep_groups <- vector("list", M)
  rep_groups_idx <- vector("list", M)
  rep_wt     <- vector("list", M)
  rep_beta   <- vector("list", M)
  rep_agree  <- numeric(M)
  obs_agree  <- numeric(M)
  rep_energy <- numeric(M)
  obs_energy <- numeric(M)
  rep_miscls <- numeric(M)
  rep_class_prop <- matrix(NA_real_, nrow = M, ncol = k)
  obs_class_prop <- matrix(NA_real_, nrow = M, ncol = k)
  colnames(rep_class_prop) <- label_set
  colnames(obs_class_prop) <- label_set
  rep_roi_agree <- matrix(NA_real_, nrow = M, ncol = J)
  obs_roi_agree <- matrix(NA_real_, nrow = M, ncol = J)
  if (verbose) {
    cat("Running PPC with", M, "posterior draws...\n")
    cat("Observed labels:", paste(label_set, collapse = ", "), "\n")
  }
  for (m in seq_along(keep)) {
    ii <- keep[m]
    theta_m <- fit$theta[ii, ]
    rho_m   <- fit$rho[ii]
    sim_m <- simulate_rep_group_dnn_roi_theta(
      theta = theta_m,
      rho = rho_m,
      x = x,
      roi = roi,
      group_raw = group_obs,
      sweeps = sweeps,
      type = type
    )
    z_rep   <- sim_m$z_rep
    z_rep_idx <- sim_m$z_rep_idx
    wt_m    <- sim_m$wt
    beta_m  <- sim_m$beta_roi
    group_idx <- sim_m$group_idx
    rep_groups[[m]] <- z_rep
    rep_groups_idx[[m]] <- z_rep_idx
    rep_wt[[m]]     <- wt_m
    rep_beta[[m]]   <- beta_m
    obs_agree[m] <- weighted_agreement_stat(group_obs, wt_m)
    rep_agree[m] <- weighted_agreement_stat(z_rep, wt_m)
    obs_class_prop[m, ] <- class_prop_stat(group_obs, label_set)
    rep_class_prop[m, ] <- class_prop_stat(z_rep, label_set)
    obs_roi_agree[m, ] <- roi_agreement_stat(group_obs, wt_m, roi, J)
    rep_roi_agree[m, ] <- roi_agreement_stat(z_rep, wt_m, roi, J)
    obs_energy[m] <- piofx_roi_shared_cpp(
      z = as.integer(group_idx),
      beta_roi = beta_m,
      roi = as.integer(roi),
      wt = wt_m
    )
    rep_energy[m] <- piofx_roi_shared_cpp(
      z = as.integer(z_rep_idx),
      beta_roi = beta_m,
      roi = as.integer(roi),
      wt = wt_m
    )
    rep_miscls[m] <- misclassification_rate_stat(group_obs, z_rep)
    if (verbose && (m %% 20 == 0 || m == M)) {
      cat("  finished", m, "/", M, "\n")
    }
  }
  pval_agree  <- mean(rep_agree >= obs_agree)
  pval_energy <- mean(rep_energy >= obs_energy)
  pval_class_prop <- numeric(k)
  names(pval_class_prop) <- label_set
  for (c in 1:k) {
    pval_class_prop[c] <- mean(rep_class_prop[, c] >= obs_class_prop[, c])
  }
  pval_roi_agree <- numeric(J)
  for (j in 1:J) {
    idx <- !is.na(rep_roi_agree[, j]) & !is.na(obs_roi_agree[, j])
    pval_roi_agree[j] <- mean(rep_roi_agree[idx, j] >= obs_roi_agree[idx, j])
  }
  out <- list(
    call = match.call(),
    keep = keep,
    label_set = label_set,
    observed = list(
      group = group_obs,
      weighted_agreement = obs_agree,
      class_prop = obs_class_prop,
      roi_agreement = obs_roi_agree,
      energy = obs_energy
    ),
    replicated = list(
      groups = rep_groups,              # original labels
      groups_idx = rep_groups_idx,      # internal contiguous labels
      wt = rep_wt,
      beta_roi = rep_beta,
      weighted_agreement = rep_agree,
      class_prop = rep_class_prop,
      roi_agreement = rep_roi_agree,
      energy = rep_energy,
      misclassification_rate_vs_obs = rep_miscls
    ),
    ppp = list(
      weighted_agreement = pval_agree,
      energy = pval_energy,
      class_prop = pval_class_prop,
      roi_agreement = pval_roi_agree
    ),
    summary = list(
      obs_agreement_mean = mean(obs_agree),
      rep_agreement_mean = mean(rep_agree),
      obs_energy_mean = mean(obs_energy),
      rep_energy_mean = mean(rep_energy),
      rep_miscls_mean = mean(rep_miscls)
    )
  )
  class(out) <- "ppc_dnn_mcmc_roi_theta"
  out
}

print.ppc_dnn_mcmc_roi_theta <- function(x, ...) {
  cat("Posterior Predictive Check for dnn_mcmc_roi_theta\n")
  cat("Number of posterior draws used:", length(x$keep), "\n")
  cat("Observed label set:", paste(x$label_set, collapse = ", "), "\n\n")
  cat("Posterior predictive p-values:\n")
  cat("  weighted_agreement :", x$ppp$weighted_agreement, "\n")
  cat("  energy             :", x$ppp$energy, "\n")
  cat("  class_prop         :", 
      paste(paste(names(x$ppp$class_prop), round(x$ppp$class_prop, 4), sep = ":"), collapse = ", "),
      "\n")
  cat("  roi_agreement      :", paste(round(x$ppp$roi_agreement, 4), collapse = ", "), "\n\n")
  cat("Summary:\n")
  cat("  mean observed agreement   :", x$summary$obs_agreement_mean, "\n")
  cat("  mean replicated agreement :", x$summary$rep_agreement_mean, "\n")
  cat("  mean observed energy      :", x$summary$obs_energy_mean, "\n")
  cat("  mean replicated energy    :", x$summary$rep_energy_mean, "\n")
  cat("  mean misclassification    :", x$summary$rep_miscls_mean, "\n")
  invisible(x)
}

ppc_dnn_mcmc_roi_theta2 <- function(
    fit,
    training_data,
    burn = 0,
    thin = 1,
    n_draws = 200,
    sweeps = 300,
    type = c("norm", "exp"),
    seed = NULL,
    verbose = TRUE
) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  n_total <- nrow(fit$theta)
  if (burn >= n_total)
    stop("burn must be smaller than nrow(fit$theta)")
  p2 <- ncol(training_data)
  x <- as.matrix(training_data[, 1:(p2 - 2), drop = FALSE])
  group_obs <- as.integer(training_data[, p2 - 1])
  roi <- as.integer(training_data[, p2])
  J <- max(roi)
  label_info <- make_label_index(group_obs)
  label_set <- label_info$label_set
  k <- length(label_set)
  
  keep_all <- seq(from = burn + 1, to = n_total, by = thin)
  if (length(keep_all) == 0)
    stop("No posterior draws left after burn/thin.")
  if (length(keep_all) > n_draws) {
    keep <- sort(sample(keep_all, n_draws))
  } else {
    keep <- keep_all
  }
  M <- length(keep)
  
  rep_groups <- vector("list", M)
  rep_groups_idx <- vector("list", M)
  rep_wt <- vector("list", M)
  rep_beta <- vector("list", M)
  rep_agree <- numeric(M)
  obs_agree <- numeric(M)
  rep_energy <- numeric(M)
  obs_energy <- numeric(M)
  rep_miscls <- numeric(M)
  rep_class_prop <- matrix(NA_real_, nrow = M, ncol = k)
  obs_class_prop <- matrix(NA_real_, nrow = M, ncol = k)
  colnames(rep_class_prop) <- label_set
  colnames(obs_class_prop) <- label_set
  rep_roi_agree <- matrix(NA_real_, nrow = M, ncol = J)
  obs_roi_agree <- matrix(NA_real_, nrow = M, ncol = J)
  colnames(rep_roi_agree) <-
    colnames(obs_roi_agree) <-
    paste0("ROI", 1:J)
  
  rep_roi_energy <- matrix(NA_real_, nrow = M, ncol = J)
  obs_roi_energy <- matrix(NA_real_, nrow = M, ncol = J)
  colnames(rep_roi_energy) <-
    colnames(obs_roi_energy) <-
    paste0("ROI", 1:J)
  obs_energy_decomp_diff <- numeric(M)
  rep_energy_decomp_diff <- numeric(M)
  if (verbose) {
    cat("Running PPC with", M, "posterior draws...\n")
    cat(
      "Observed labels:",
      paste(label_set, collapse = ", "),
      "\n"
    )
    cat("Number of ROIs:", J, "\n")
  }
  
  for (m in seq_along(keep)) {
    ii <- keep[m]
    theta_m <- fit$theta[ii, ]
    rho_m <- fit$rho[ii]
    sim_m <- simulate_rep_group_dnn_roi_theta(
      theta = theta_m,
      rho = rho_m,
      x = x,
      roi = roi,
      group_raw = group_obs,
      sweeps = sweeps,
      type = type
    )
    z_rep <- sim_m$z_rep
    z_rep_idx <- sim_m$z_rep_idx
    wt_m <- sim_m$wt
    beta_m <- sim_m$beta_roi
    group_idx <- sim_m$group_idx
    rep_groups[[m]] <- z_rep
    rep_groups_idx[[m]] <- z_rep_idx
    rep_wt[[m]] <- wt_m
    rep_beta[[m]] <- beta_m
    obs_agree[m] <- weighted_agreement_stat(group_obs, wt_m)
    rep_agree[m] <- weighted_agreement_stat(z_rep, wt_m)
    obs_class_prop[m, ] <- class_prop_stat(group_obs, label_set)
    rep_class_prop[m, ] <- class_prop_stat(z_rep, label_set)
    obs_roi_agree[m, ] <- roi_agreement_stat(group_obs, wt_m, roi, J)
    rep_roi_agree[m, ] <- roi_agreement_stat(z_rep, wt_m, roi,  J)
    obs_energy[m] <- piofx_roi_shared_cpp(
      z = as.integer(group_idx),
      beta_roi = beta_m,
      roi = as.integer(roi),
      wt = wt_m
    )
    rep_energy[m] <- piofx_roi_shared_cpp(
      z = as.integer(z_rep_idx),
      beta_roi = beta_m,
      roi = as.integer(roi),
      wt = wt_m
    )
    obs_roi_energy[m, ] <- roi_energy_stat(
      group = group_obs,
      wt = wt_m,
      beta_roi = beta_m,
      roi = roi,
      J = J
    )
    rep_roi_energy[m, ] <- roi_energy_stat(
      group = z_rep,
      wt = wt_m,
      beta_roi = beta_m,
      roi = roi,
      J = J
    )
    obs_energy_decomp_diff[m] <-
      obs_energy[m] -
      sum(obs_roi_energy[m, ], na.rm = TRUE)
    rep_energy_decomp_diff[m] <-
      rep_energy[m] -
      sum(rep_roi_energy[m, ], na.rm = TRUE)
    rep_miscls[m] <- misclassification_rate_stat(
      group_obs,
      z_rep
    )
    if (verbose &&(m %% 20 == 0 || m == M)) {
      cat("  finished", m, "/", M, "\n")
    }
  }
  pval_agree <- mean(
    rep_agree >= obs_agree,
    na.rm = TRUE
  )
  pval_energy <- mean(
    rep_energy >= obs_energy,
    na.rm = TRUE
  )
  pval_class_prop <- numeric(k)
  names(pval_class_prop) <- label_set
  for (c in 1:k) {
    pval_class_prop[c] <- mean(
      rep_class_prop[, c] >=
        obs_class_prop[, c],
      na.rm = TRUE
    )
  }
  pval_roi_agree <- numeric(J)
  names(pval_roi_agree) <- paste0("ROI", 1:J)
  for (j in 1:J) {
    idx <- !is.na(rep_roi_agree[, j]) &
      !is.na(obs_roi_agree[, j])
    pval_roi_agree[j] <- if (any(idx)) {
      mean(rep_roi_agree[idx, j] >= obs_roi_agree[idx, j])
    } else {
      NA_real_
    }
  }
  
  pval_roi_energy <- numeric(J)
  names(pval_roi_energy) <- paste0("ROI", 1:J)
  for (j in 1:J) {
    idx <- !is.na(rep_roi_energy[, j]) &
      !is.na(obs_roi_energy[, j])
    pval_roi_energy[j] <- if (any(idx)) {
      mean(rep_roi_energy[idx, j] >= obs_roi_energy[idx, j])
    } else {
      NA_real_
    }
  }
  out <- list(
    call = match.call(),
    keep = keep,
    label_set = label_set,
    observed = list(
      group = group_obs,
      weighted_agreement = obs_agree,
      class_prop = obs_class_prop,
      roi_agreement = obs_roi_agree,
      energy = obs_energy,
      roi_energy = obs_roi_energy),
    replicated = list(
      groups = rep_groups,
      groups_idx = rep_groups_idx,
      wt = rep_wt,
      beta_roi = rep_beta,
      weighted_agreement = rep_agree,
      class_prop = rep_class_prop,
      roi_agreement = rep_roi_agree,
      energy = rep_energy,
      roi_energy = rep_roi_energy,
      misclassification_rate_vs_obs = rep_miscls),
    ppp = list(
      weighted_agreement = pval_agree,
      energy = pval_energy,
      class_prop = pval_class_prop,
      roi_agreement = pval_roi_agree,
      roi_energy = pval_roi_energy),
    energy_decomposition_check = list(
      observed_difference = obs_energy_decomp_diff,
      replicated_difference = rep_energy_decomp_diff,
      max_abs_observed_difference =
        max(abs(obs_energy_decomp_diff), na.rm = TRUE),
      max_abs_replicated_difference =
        max(abs(rep_energy_decomp_diff), na.rm = TRUE)
    ),
    summary = list(
      obs_agreement_mean = mean(obs_agree, na.rm = TRUE),
      rep_agreement_mean = mean(rep_agree, na.rm = TRUE),
      obs_energy_mean = mean(obs_energy, na.rm = TRUE),
      rep_energy_mean = mean(rep_energy, na.rm = TRUE),
      obs_roi_energy_mean = colMeans(obs_roi_energy, na.rm = TRUE),
      rep_roi_energy_mean = colMeans(rep_roi_energy, na.rm = TRUE),
      rep_miscls_mean = mean(rep_miscls, na.rm = TRUE))
  )
  class(out) <- "ppc_dnn_mcmc_roi_theta"
  out
}

plot_ppc_dnn_mcmc_roi_theta <- function(ppc_res) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(2, 2))
  hist(ppc_res$replicated$weighted_agreement,
       breaks = 30,
       main = "PPC: weighted agreement",
       xlab = "replicated")
  abline(v = mean(ppc_res$observed$weighted_agreement), col = 2, lwd = 2)
  hist(ppc_res$replicated$energy,
       breaks = 30,
       main = "PPC: energy",
       xlab = "replicated")
  abline(v = mean(ppc_res$observed$energy), col = 2, lwd = 2)
  boxplot(as.data.frame(ppc_res$replicated$roi_agreement),
          main = "PPC: ROI agreement",
          xlab = "ROI")
  points(1:ncol(ppc_res$replicated$roi_agreement),
         colMeans(ppc_res$observed$roi_agreement, na.rm = TRUE),
         col = 2, pch = 19)
  matplot(t(ppc_res$replicated$class_prop),
          type = "l", lty = 1,
          main = "PPC: class proportions",
          xlab = "class index", ylab = "proportion")
  points(1:ncol(ppc_res$replicated$class_prop),
         colMeans(ppc_res$observed$class_prop),
         col = 2, pch = 19)
  axis(1, at = 1:ncol(ppc_res$replicated$class_prop), labels = colnames(ppc_res$replicated$class_prop))
}

ppc_res <- ppc_dnn_mcmc_roi_theta2(
  fit = fit,
  training_data = train,
  burn = 10000,
  thin = 10,
  n_draws = 1000,
  sweeps = 100,
  type = "norm",
  seed = 1
)

ppc_res$replicated$wt <- NULL
ppc_res$replicated$groups <- NULL
ppc_res$replicated$groups_idx <- NULL

ppc_vals <- sapply(1:8, function(i) {
  mean(ppc_res$replicated$roi_energy[, i] >
         mean(ppc_res$observed$roi_energy[, i]))
})

sqrt(sum((ppc_vals - 0.5)^2))

save(fit, mcmc_time, pred, pred_time, posterior_summary, diag_res, ppc_res, file=output_file)


####### connectivity analysis
compute_posterior_connectivity <- function(
    fit,
    training_data,
    burn = 10000,
    thin = 10,
    type = "norm",
    symmetrize = TRUE
){
  p2 <- ncol(training_data)
  x   <- as.matrix(training_data[, 1:(p2 - 2), drop = FALSE])
  y   <- as.integer(training_data[, p2 - 1])
  roi <- as.integer(training_data[, p2])
  J <- max(roi)
  n <- nrow(training_data)
  idx <- seq(burn + 1, length(fit$rho), by = thin)
  conn_arr <- array(NA_real_, dim = c(length(idx), J, J))
  D <- nearest_neighs_k_cpp(x, y)
  for(s in seq_along(idx)){
    t <- idx[s]
    wt <- dist2wt_cpp(
      D = D,
      rho = fit$rho[t],
      n = n,
      type = type
    )
    C <- matrix(0, J, J)
    for(j in 1:J){
      idx_j <- which(roi == j)
      for(k in 1:J){
        idx_k <- which(roi == k)
        C[j, k] <- mean(wt[idx_j, idx_k, drop = FALSE])
      }
    }
    if(symmetrize){
      C <- (C + t(C)) / 2
    }
    conn_arr[s, , ] <- C
  }
  conn_mean <- apply(conn_arr, c(2, 3), mean, na.rm = TRUE)
  conn_low  <- apply(conn_arr, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
  conn_high <- apply(conn_arr, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)
  rownames(conn_mean) <- colnames(conn_mean) <- paste0("ROI", 1:J)
  rownames(conn_low)  <- colnames(conn_low)  <- paste0("ROI", 1:J)
  rownames(conn_high) <- colnames(conn_high) <- paste0("ROI", 1:J)
  list(
    mean = conn_mean,
    lower = conn_low,
    upper = conn_high,
    draws = conn_arr
  )
}


conn_res <- compute_posterior_connectivity(
  fit = fit,
  training_data = train,
  burn = 10000,
  thin = 10,
  type = "norm",
  symmetrize = TRUE
)


save(fit, mcmc_time, pred, pred_time, posterior_summary, diag_res, ppc_res, conn_res, file=output_file)
 

plot_connectivity_upper <- function(
    conn_mat,
    title = "Standardized ROI similarity",
    standardize = TRUE,
    digits = 2,
    text_size = 4.8,
    zlim = NULL
){
  J <- nrow(conn_mat)
  roi_labels <- c(
    "1L", "2L", "3L", "4L",
    "1R", "2R", "3R", "4R"
  )
  if(J != length(roi_labels)){
    stop("The number of ROIs must match the length of roi_labels.")
  }
  ij <- which(
    upper.tri(conn_mat),
    arr.ind = TRUE
  )
  df <- data.frame(
    ROI_from = ij[, 1],
    ROI_to   = ij[, 2],
    similarity = conn_mat[ij]
  )
  if(standardize){
    mu <- mean(df$similarity, na.rm = TRUE)
    s  <- sd(df$similarity, na.rm = TRUE)
    df$similarity <- (df$similarity - mu) / s
  }
  df$ROI_from <- factor(
    roi_labels[df$ROI_from],
    levels = roi_labels
  )
  df$ROI_to <- factor(
    roi_labels[df$ROI_to],
    levels = roi_labels
  )
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = ROI_to,
      y = ROI_from,
      fill = similarity
    )
  ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf(
          paste0("%.", digits, "f"),
          similarity
        )
      ),
      size = text_size
    )
  if(standardize){
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        limits = zlim,
        name = "Standardized\nsimilarity"
      )
  } else {
    p <- p +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "red",
        name = "Similarity"
      )
  }
  p +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        size = 14
      ),
      axis.text.y = ggplot2::element_text(
        size = 14
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 14
      )
    )
}

plot_connectivity_upper(
  conn_res$mean
)


plot_connectivity_networks <- function(
    conn_res,
    roi_names = c(
      "1L", "2L", "3L", "4L",
      "1R", "2R", "3R", "4R"
    ),
    top_n = 5,
    digits = 2
){
  
  J <- nrow(conn_res$mean)
  
  if(length(roi_names) != J){
    stop("length(roi_names) must equal the number of ROIs.")
  }
  pair_idx <- which(
    upper.tri(conn_res$mean),
    arr.ind = TRUE
  )
  pair_df <- data.frame(
    j = pair_idx[, 1],
    k = pair_idx[, 2]
  )
  pair_df$ROI1 <- roi_names[pair_df$j]
  pair_df$ROI2 <- roi_names[pair_df$k]
  pair_df$mean <- mapply(
    function(j, k) conn_res$mean[j, k],
    pair_df$j,
    pair_df$k
  )
  pair_df$lower <- mapply(
    function(j, k) conn_res$lower[j, k],
    pair_df$j,
    pair_df$k
  )
  pair_df$upper <- mapply(
    function(j, k) conn_res$upper[j, k],
    pair_df$j,
    pair_df$k
  )
  scale_mean <- mean(pair_df$mean, na.rm = TRUE)
  scale_sd <- sd(pair_df$mean, na.rm = TRUE)
  if(!is.finite(scale_sd) || scale_sd == 0){
    stop("SD of ROI-pair similarities is zero or non-finite.")
  }
  pair_df$mean_raw  <- pair_df$mean
  pair_df$lower_raw <- pair_df$lower
  pair_df$upper_raw <- pair_df$upper
  pair_df$mean <- (pair_df$mean_raw - scale_mean) / scale_sd
  pair_df$lower <- (pair_df$lower_raw - scale_mean) / scale_sd
  pair_df$upper <- (pair_df$upper_raw - scale_mean) / scale_sd
  top_df <- pair_df[
    order(pair_df$mean, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  top_df <- head(
    top_df,
    top_n
  )
  lr_pairs <- data.frame(
    j = 1:4,
    k = 5:8
  )
  lr_df <- do.call(
    rbind,
    lapply(seq_len(nrow(lr_pairs)), function(i){
      pair_df[
        pair_df$j == lr_pairs$j[i] &
          pair_df$k == lr_pairs$k[i],
        ,
        drop = FALSE
      ]
      
    })
  )
  nodes <- data.frame(
    ROI = roi_names,
    x = c(
      # 1L, 2L, 3L, 4L
      0.9,  1.7,  1.7,  0.9,
      # 1R, 2R, 3R, 4R
      -0.9, -1.7, -1.7, -0.9
    ),
    y = c(
      # L
      1.5,  0.6, -0.6, -1.5,
      # R
      1.5,  0.6, -0.6, -1.5
    )
  )
  make_network <- function(edge_df, title){
    edge_df$x1 <- nodes$x[
      match(edge_df$ROI1, nodes$ROI)
    ]
    edge_df$y1 <- nodes$y[
      match(edge_df$ROI1, nodes$ROI)
    ]
    edge_df$x2 <- nodes$x[
      match(edge_df$ROI2, nodes$ROI)
    ]
    edge_df$y2 <- nodes$y[
      match(edge_df$ROI2, nodes$ROI)
    ]
    edge_df$xmid <- (
      edge_df$x1 + edge_df$x2
    ) / 2
    edge_df$ymid <- (
      edge_df$y1 + edge_df$y2
    ) / 2
    edge_df$label <- sprintf(
      paste0(
        "%.", digits, "f\n",
        "(%.", digits, "f, %.", digits, "f)"
      ),
      edge_df$mean,
      edge_df$lower,
      edge_df$upper
    )
    p <- ggplot() +
      geom_segment(
        data = edge_df,
        aes(
          x = x1,
          y = y1,
          xend = x2,
          yend = y2
        ),
        linewidth = 0.9,
        color = "gray35"
      ) +
      geom_label(
        data = edge_df,
        aes(
          x = xmid,
          y = ymid,
          label = label
        ),
        size = 3.4,
        lineheight = 0.9,
        label.size = 0.15,
        fill = "white"
      ) +
      geom_point(
        data = nodes,
        aes(
          x = x,
          y = y
        ),
        size = 11,
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 1
      ) +
      geom_text(
        data = nodes,
        aes(
          x = x,
          y = y,
          label = ROI
        ),
        size = 4,
        fontface = "bold"
      ) +
      coord_equal(
        xlim = c(-2.3, 2.3),
        ylim = c(-2.0, 2.0),
        clip = "off"
      ) +
      labs(
        title = title
      ) +
      theme_void(
        base_size = 13
      ) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 14
        ),
        plot.margin = margin(
          10, 20, 10, 20
        )
      )
    list(
      plot = p,
      data = edge_df
    )
  }
  top_result <- make_network(
    edge_df = top_df,
    title = paste0(
      "Top ", top_n,
      " standardized ROI-pair similarities"
    )
  )
  lr_result <- make_network(
    edge_df = lr_df,
    title = "Standardized left-right homologous ROI similarities"
  )
  return(
    list(
      top5_plot = top_result$plot,
      leftright_plot = lr_result$plot,
      top5_data = top_result$data,
      leftright_data = lr_result$data,
      all_pairs = pair_df,
      scale_mean = scale_mean,
      scale_sd = scale_sd,
      nodes = nodes
    )
  )
}

network_res <- plot_connectivity_networks(
  conn_res = conn_res,
  top_n = 5,
  digits = 2
)

network_res$top5_plot

