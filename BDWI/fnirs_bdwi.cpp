// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ============================================================
// utilities
// ============================================================

static inline double log_norm_density(double x, double mu, double s2){
  return -0.5 * std::log(2.0 * M_PI * s2) - 0.5 * (x - mu) * (x - mu) / s2;
}

static inline double mean_cpp(const NumericVector& x){
  int n = x.size();
  double s = 0.0;
  for(int i=0;i<n;++i) s += x[i];
  return s / n;
}

static inline NumericVector center_cpp(const NumericVector& x){
  int n = x.size();
  double m = mean_cpp(x);
  NumericVector out(n);
  for(int i=0;i<n;++i) out[i] = x[i] - m;
  return out;
}

static inline double ssq_cpp(const NumericVector& x){
  int n = x.size();
  double s = 0.0;
  for(int i=0;i<n;++i) s += x[i] * x[i];
  return s;
}

// prior induced by
// alpha = mean(theta) ~ N(mu_alpha, s2_alpha)
// u = theta - mean(theta), centered Gaussian with tau2
static inline double log_prior_theta_cpp(
    const NumericVector& theta,
    double mu_alpha,
    double s2_alpha,
    double tau2
){
  double alpha = mean_cpp(theta);
  NumericVector u = center_cpp(theta);
  
  double lp = log_norm_density(alpha, mu_alpha, s2_alpha)
    - 0.5 * ssq_cpp(u) / tau2;
  
  return lp;
}

// ============================================================
// distance functions
// ============================================================

// [[Rcpp::export]]
NumericMatrix nearest_neighs_k_cpp(const NumericMatrix& x,
                                   const IntegerVector& group){
  int n = x.nrow(), p = x.ncol();
  NumericMatrix D(n,n);
  const double BIG = 1e5;
  
  for(int i=0;i<n;++i){
    for(int j=i+1;j<n;++j){
      double s=0.0;
      for(int t=0;t<p;++t){
        double d = x(i,t)-x(j,t);
        s += d*d;
      }
      double dist = std::sqrt(s);
      D(i,j)=D(j,i)=dist;
    }
    D(i,i)=BIG;
  }
  return D;
}

// [[Rcpp::export]]
NumericMatrix nearest_neighsy_k_cpp(const NumericMatrix& y,
                                    const NumericMatrix& x,
                                    const IntegerVector& group){
  int m=y.nrow(), n=x.nrow(), p=x.ncol();
  NumericMatrix D(m,n);
  
  for(int i=0;i<m;++i){
    for(int j=0;j<n;++j){
      double s=0.0;
      for(int t=0;t<p;++t){
        double d=y(i,t)-x(j,t);
        s+=d*d;
      }
      D(i,j)=std::sqrt(s);
    }
  }
  return D;
}

// ============================================================
// weights
// ============================================================

// [[Rcpp::export]]
NumericMatrix dist2wt_cpp(const NumericMatrix& D,
                          double rho, int n,
                          std::string type){
  if (D.nrow() != n || D.ncol() != n)
    stop("dist2wt_cpp: D dim mismatch");
  if (rho <= 0.0) stop("dist2wt_cpp: rho must be positive");
  
  NumericMatrix W(n, n);
  double inv_rho = 1.0 / rho;
  
  for(int i=0;i<n;++i){
    std::vector<double> scr(n);
    double mx = -INFINITY;
    
    for(int j=0;j<n;++j){
      double d = D(i,j);
      double v;
      if (type == "exp") {
        v = -d * inv_rho;
      } else {
        double z = d * inv_rho;
        v = -0.5 * z * z;
      }
      scr[j] = v;
      if (v > mx) mx = v;
    }
    
    double Z = 0.0;
    for(int j=0;j<n;++j){
      double val = std::exp(scr[j] - mx);
      W(i,j) = val;
      Z += val;
    }
    
    double invZ = (Z > 0.0 && R_finite(Z)) ? 1.0 / Z : 1.0 / n;
    for(int j=0;j<n;++j) W(i,j) *= invZ;
  }
  
  return W;
}

// [[Rcpp::export]]
NumericMatrix dist2wt_y_cpp(const NumericMatrix& D,
                            double rho, int n, std::string type){
  int m=D.nrow();
  if (rho <= 0.0) stop("dist2wt_y_cpp: rho must be positive");
  
  NumericMatrix W(m,n);
  
  for(int i=0;i<m;++i){
    std::vector<double> scr(n);
    double mx=-INFINITY;
    
    for(int j=0;j<n;++j){
      double v=(type=="exp") ? -D(i,j)/rho
      : -0.5*(D(i,j)/rho)*(D(i,j)/rho);
      scr[j]=v;
      if(v>mx) mx=v;
    }
    
    double Z=0.0;
    for(int j=0;j<n;++j){
      W(i,j)=std::exp(scr[j]-mx);
      Z+=W(i,j);
    }
    
    double inv=(Z>0)?1.0/Z:1.0/n;
    for(int j=0;j<n;++j) W(i,j)*=inv;
  }
  
  return W;
}

// ============================================================
// ROI-shared helpers under theta parameterization
// beta_j = exp(theta_j)
// ============================================================
// [[Rcpp::export]]
NumericVector make_beta_roi_from_theta(const NumericVector& theta){
  int J = theta.size();
  NumericVector beta(J);
  
  for(int j = 0; j < J; ++j){
    beta[j] = std::exp(theta[j]);
  }
  
  return beta;
}



// [[Rcpp::export]]
double piofx_roi_shared_cpp(
    const IntegerVector& z,
    const NumericVector& beta_roi,
    const IntegerVector& roi,
    const NumericMatrix& wt
){
  int n = z.size();
  if (roi.size() != n) stop("roi length mismatch");
  if (wt.nrow()!=n || wt.ncol()!=n) stop("wt dim mismatch");
  
  int J = beta_roi.size();
  double s = 0.0;
  
  for(int i=0;i<n;++i){
    int zi = z[i];
    int r  = roi[i];
    if(r < 1 || r > J) stop("roi index out of range");
    
    double bi = beta_roi[r - 1];
    
    for(int j=0;j<n;++j){
      if(i==j) continue;
      if(z[j] == zi) s += bi * wt(i,j);
    }
  }
  
  return s;
}

// [[Rcpp::export]]
IntegerVector gibbs_draw_k1_roi_shared_cpp(
    int k,
    int n,
    const NumericVector& beta_roi,
    int sweeps,
    const NumericMatrix& wt,
    const IntegerVector& roi,
    const IntegerVector& z0
){
  RNGScope rng;
  
  if(z0.size() != n)  stop("z0 length mismatch");

  if(roi.size() != n)  stop("roi length mismatch");
  
  if(wt.nrow() != n || wt.ncol() != n)  stop("wt dim mismatch");
  
  int J = beta_roi.size();
  
  IntegerVector z = clone(z0);
  
  
  for(int it = 0; it < sweeps; ++it){
    
    for(int i = 0; i < n; ++i){
      
      // --------------------------------------------------
      // propose a different class
      // --------------------------------------------------
      
      int prop;
      
      do{
        prop =  1 +  (int)std::floor(k * ::unif_rand());
        
      } while(prop == z[i]);
      
      int zi = z[i];
      
      int ri = roi[i];
      
      if(ri < 1 || ri > J)
        stop("roi index out of range");
      
      double beta_i =
        beta_roi[ri - 1];
      // --------------------------------------------------
      // Exact change in the joint energy
      //
      // Delta H =
      //
      // beta_{r_i} sum_j w_ij Delta I_j
      //
      // +
      //
      // sum_j beta_{r_j} w_ji Delta I_j
      //
      // --------------------------------------------------
      
      double delta_out = 0.0;
      double delta_in  = 0.0;
      
      
      for(int j = 0; j < n; ++j){
        
        if(i == j)
          continue;
        
        
        // Delta indicator:
        //
        // I(z_j = proposed)
        // -
        // I(z_j = current)
        
        double diff = 0.0;
        
        if(z[j] == prop)
          diff += 1.0;
        
        if(z[j] == zi)
          diff -= 1.0;
        if(diff == 0.0)
          continue;
        // ----------------------------------------------
        // outgoing:
        //
        // beta_{r_i} * w_ij
        // ----------------------------------------------
        
        delta_out +=  beta_i *  wt(i, j) *  diff;
        // ----------------------------------------------
        // incoming:
        //
        // beta_{r_j} * w_ji
        // ----------------------------------------------
        
        int rj = roi[j];
        
        if(rj < 1 || rj > J)
          stop("roi index out of range");
        
        double beta_j =  beta_roi[rj - 1];
        delta_in +=  beta_j *  wt(j, i) *  diff;
      }
      // --------------------------------------------------
      // Exact MH log acceptance ratio
      // --------------------------------------------------
      double log_alpha = delta_out +  delta_in;
      
      if(
        std::log(::unif_rand())<= log_alpha
      ){
        z[i] = prop;
      }
    }
  }
  return z;
}
// ============================================================
// exchange update for rho | theta
// prior:
//   log rho ~ N(mu_rho, s2_rho)
// ============================================================

// [[Rcpp::export]]
List exchange_move_rho_theta_cpp_adapt(
    int k,
    double rho,
    const NumericVector& theta,
    const NumericMatrix& neigh_dist,
    const IntegerVector& group,
    const IntegerVector& roi,
    int n,
    int loop_aux,
    int accept,
    NumericMatrix wt,
    std::string type,
    double sd_log_rho,
    double mu_rho,
    double s2_rho,
    int iter,
    int burnin,
    double target_acc = 0.25
){
  RNGScope rng;
  
  if(rho <= 0.0) stop("rho must be positive");
  
  double log_rho  = std::log(rho);
  double log_rho2 = log_rho + R::rnorm(0.0, sd_log_rho);
  double rho2     = std::exp(log_rho2);
  
  if(!R_finite(rho2) || rho2 <= 0.0) rho2 = 1e-12;
  
  NumericMatrix wt2 = dist2wt_cpp(neigh_dist, rho2, n, type);
  
  NumericVector beta_roi = make_beta_roi_from_theta(theta);
  
  IntegerVector g2 =
    gibbs_draw_k1_roi_shared_cpp(
      k, n, beta_roi, loop_aux, wt2, roi, group
    );
  
  double log_like_ratio =
    piofx_roi_shared_cpp(group, beta_roi, roi, wt2)
    + piofx_roi_shared_cpp(g2,    beta_roi, roi, wt )
    - piofx_roi_shared_cpp(group, beta_roi, roi, wt )
    - piofx_roi_shared_cpp(g2,    beta_roi, roi, wt2);
    
    double log_prior_ratio =
    log_norm_density(log_rho2, mu_rho, s2_rho)
      - log_norm_density(log_rho,  mu_rho, s2_rho);
    
    double log_alpha = log_like_ratio + log_prior_ratio;
    
    int accepted = 0;
    
    if(log_alpha >= std::log(::unif_rand())){
      rho = rho2;
      wt  = wt2;
      accept++;
      accepted = 1;
    }
    
    if(iter < burnin){
      double eta = 1.0 / std::sqrt(iter + 10.0);
      sd_log_rho =
        std::exp(std::log(sd_log_rho) + eta * (accepted - target_acc));
    }
    
    return List::create(
      _["rho"]         = rho,
      _["wt"]          = wt,
      _["accept"]      = accept,
      _["sd_log_rho"]  = sd_log_rho,
      _["accepted"]    = accepted
    );
}

// ============================================================
// one-at-a-time update for theta_j | rho, tau2
// induced prior:
//   alpha = mean(theta) ~ N(mu_alpha, s2_alpha)
//   u = theta - mean(theta), centered Gaussian with tau2
// ============================================================

// [[Rcpp::export]]
List exchange_move_theta_single_cpp_adapt(
    int k,
    NumericVector theta,
    int j_update_1based,
    double tau2,
    double mu_alpha,
    double s2_alpha,
    const IntegerVector& group,
    const IntegerVector& roi,
    int n,
    int loop_aux,
    int accept,
    NumericMatrix wt,
    double sd_theta,
    int iter,
    int burnin,
    double target_acc = 0.25
){
  RNGScope rng;
  
  int J = theta.size();
  int j = j_update_1based - 1;
  
  if(j < 0 || j >= J) stop("j_update out of range");
  if(tau2 <= 0.0) stop("tau2 must be positive");
  if(s2_alpha <= 0.0) stop("s2_alpha must be positive");
  
  NumericVector theta2 = clone(theta);
  theta2[j] = theta[j] + R::rnorm(0.0, sd_theta);
  
  NumericVector beta_roi  = make_beta_roi_from_theta(theta);
  NumericVector beta2_roi = make_beta_roi_from_theta(theta2);
  
  IntegerVector g2 =
    gibbs_draw_k1_roi_shared_cpp(
      k, n, beta2_roi, loop_aux, wt, roi, group
    );
  
  double log_like_ratio =
    piofx_roi_shared_cpp(group, beta2_roi, roi, wt)
    + piofx_roi_shared_cpp(g2,    beta_roi,  roi, wt)
    - piofx_roi_shared_cpp(group, beta_roi,  roi, wt)
    - piofx_roi_shared_cpp(g2,    beta2_roi, roi, wt);
    
    double log_prior_ratio =
    log_prior_theta_cpp(theta2, mu_alpha, s2_alpha, tau2)
      - log_prior_theta_cpp(theta,  mu_alpha, s2_alpha, tau2);
    
    double log_alpha = log_like_ratio + log_prior_ratio;
    
    int accepted = 0;
    
    if(log_alpha >= std::log(::unif_rand())){
      theta = theta2;
      accept++;
      accepted = 1;
    }
    
    if(iter < burnin){
      double eta = 1.0 / std::sqrt(iter + 10.0);
      sd_theta =
        std::exp(std::log(sd_theta) + eta * (accepted - target_acc));
    }
    
    return List::create(
      _["theta"]      = theta,
      _["accept"]     = accept,
      _["sd_theta"]   = sd_theta,
      _["accepted"]   = accepted
    );
}

// ============================================================
// Regularized Horseshoe helpers
// theta_j = alpha + u_j, alpha = mean(theta), u = centered(theta)
// u_j | tau, lambda_j, c ~ N(0, tau^2 * \tilde lambda_j^2)
// \tilde lambda_j^2 = c^2 lambda_j^2 / (c^2 + tau^2 lambda_j^2)
// ------------------------------------------------------------
// lambda_j ~ Half-Cauchy(0, 1)
// tau      ~ Half-Cauchy(0, scale_tau)
// c^2 is treated as fixed here
// ============================================================

static inline double log_half_cauchy_density(double x, double scale){
  if(x <= 0.0 || scale <= 0.0) return R_NegInf;
  double z = x / scale;
  return std::log(2.0) - std::log(M_PI) - std::log(scale) - std::log1p(z * z);
}

static inline double rhs_lambda_tilde2(double lambda2, double tau2, double c2){
  if(lambda2 <= 0.0) stop("rhs_lambda_tilde2: lambda2 must be positive");
  if(tau2    <= 0.0) stop("rhs_lambda_tilde2: tau2 must be positive");
  if(c2      <= 0.0) stop("rhs_lambda_tilde2: c2 must be positive");
  
  return (c2 * lambda2) / (c2 + tau2 * lambda2);
}

// prior kernel under regularized horseshoe
// alpha = mean(theta) ~ N(mu_alpha, s2_alpha)
// u_j = theta_j - mean(theta)
// u_j | tau, lambda_j, c ~ N(0, tau2 * lambda_tilde_j^2)
//
// NOTE:
// As in the original centered-theta code, this is used as a practical
// prior kernel on centered u. For MH ratios this is usually what matters.
static inline double log_prior_theta_rhs_cpp(
    const NumericVector& theta,
    double mu_alpha,
    double s2_alpha,
    double tau2,
    const NumericVector& lambda2,
    double c2
){
  int J = theta.size();
  if(lambda2.size() != J) stop("log_prior_theta_rhs_cpp: lambda2 length mismatch");
  if(s2_alpha <= 0.0) stop("log_prior_theta_rhs_cpp: s2_alpha must be positive");
  if(tau2 <= 0.0)     stop("log_prior_theta_rhs_cpp: tau2 must be positive");
  if(c2 <= 0.0)       stop("log_prior_theta_rhs_cpp: c2 must be positive");
  
  double alpha = mean_cpp(theta);
  NumericVector u = center_cpp(theta);
  
  double lp = log_norm_density(alpha, mu_alpha, s2_alpha);
  
  for(int j=0; j<J; ++j){
    if(lambda2[j] <= 0.0) stop("log_prior_theta_rhs_cpp: lambda2 must be positive");
    
    double ltil2 = rhs_lambda_tilde2(lambda2[j], tau2, c2);
    double s2j   = tau2 * ltil2;
    
    lp += -0.5 * std::log(s2j) - 0.5 * u[j] * u[j] / s2j;
  }
  
  return lp;
}

// [[Rcpp::export]]
NumericVector rhs_effective_var_cpp(
    double tau2,
    const NumericVector& lambda2,
    double c2
){
  int J = lambda2.size();
  NumericVector out(J);
  
  for(int j=0; j<J; ++j){
    double ltil2 = rhs_lambda_tilde2(lambda2[j], tau2, c2);
    out[j] = tau2 * ltil2;
  }
  
  return out;
}

// ============================================================
// one-at-a-time theta_j update under regularized horseshoe
// ============================================================

// [[Rcpp::export]]
List exchange_move_theta_single_rhs_cpp_adapt(
    int k,
    NumericVector theta,
    int j_update_1based,
    double tau2,
    const NumericVector& lambda2,
    double c2,
    double mu_alpha,
    double s2_alpha,
    const IntegerVector& group,
    const IntegerVector& roi,
    int n,
    int loop_aux,
    int accept,
    NumericMatrix wt,
    double sd_theta,
    int iter,
    int burnin,
    double target_acc = 0.25
){
  RNGScope rng;
  
  int J = theta.size();
  int j = j_update_1based - 1;
  
  if(j < 0 || j >= J) stop("j_update out of range");
  if(lambda2.size() != J) stop("lambda2 length mismatch");
  if(tau2 <= 0.0) stop("tau2 must be positive");
  if(c2 <= 0.0)   stop("c2 must be positive");
  if(s2_alpha <= 0.0) stop("s2_alpha must be positive");
  
  NumericVector theta2 = clone(theta);
  theta2[j] = theta[j] + R::rnorm(0.0, sd_theta);
  
  NumericVector beta_roi  = make_beta_roi_from_theta(theta);
  NumericVector beta2_roi = make_beta_roi_from_theta(theta2);
  
  IntegerVector g2 =
    gibbs_draw_k1_roi_shared_cpp(
      k, n, beta2_roi, loop_aux, wt, roi, group
    );
  
  double log_like_ratio =
    piofx_roi_shared_cpp(group, beta2_roi, roi, wt)
    + piofx_roi_shared_cpp(g2,    beta_roi,  roi, wt)
    - piofx_roi_shared_cpp(group, beta_roi,  roi, wt)
    - piofx_roi_shared_cpp(g2,    beta2_roi, roi, wt);
    
    double log_prior_ratio =
    log_prior_theta_rhs_cpp(theta2, mu_alpha, s2_alpha, tau2, lambda2, c2)
      - log_prior_theta_rhs_cpp(theta,  mu_alpha, s2_alpha, tau2, lambda2, c2);
    
    double log_alpha = log_like_ratio + log_prior_ratio;
    
    int accepted = 0;
    
    if(log_alpha >= std::log(::unif_rand())){
      theta = theta2;
      accept++;
      accepted = 1;
    }
    
    if(iter < burnin){
      double eta = 1.0 / std::sqrt(iter + 10.0);
      sd_theta =
        std::exp(std::log(sd_theta) + eta * (accepted - target_acc));
    }
    
    return List::create(
      _["theta"]      = theta,
      _["accept"]     = accept,
      _["sd_theta"]   = sd_theta,
      _["accepted"]   = accepted
    );
}

// ============================================================
// one-at-a-time lambda_j update under regularized horseshoe
// proposal is on log(lambda_j), target includes Jacobian +log(lambda_j)
// ============================================================

// [[Rcpp::export]]
List exchange_move_lambda_single_rhs_cpp_adapt(
    NumericVector lambda2,
    int j_update_1based,
    const NumericVector& theta,
    double tau2,
    double c2,
    double mu_alpha,
    double s2_alpha,
    int accept,
    double sd_log_lambda,
    int iter,
    int burnin,
    double target_acc = 0.25
){
  RNGScope rng;
  
  int J = lambda2.size();
  int j = j_update_1based - 1;
  
  if(j < 0 || j >= J) stop("j_update out of range");
  if(theta.size() != J) stop("theta and lambda2 length mismatch");
  if(tau2 <= 0.0) stop("tau2 must be positive");
  if(c2 <= 0.0)   stop("c2 must be positive");
  
  double lambda_cur = std::sqrt(lambda2[j]);
  if(lambda_cur <= 0.0 || !R_finite(lambda_cur))
    stop("current lambda must be positive and finite");
  
  double log_lambda_cur = std::log(lambda_cur);
  double log_lambda_new = log_lambda_cur + R::rnorm(0.0, sd_log_lambda);
  double lambda_new     = std::exp(log_lambda_new);
  double lambda2_new_j  = lambda_new * lambda_new;
  
  if(!R_finite(lambda_new) || lambda_new <= 0.0) lambda_new = 1e-12;
  
  NumericVector lambda2_prop = clone(lambda2);
  lambda2_prop[j] = lambda2_new_j;
  
  double log_prior_theta_ratio =
    log_prior_theta_rhs_cpp(theta, mu_alpha, s2_alpha, tau2, lambda2_prop, c2)
    - log_prior_theta_rhs_cpp(theta, mu_alpha, s2_alpha, tau2, lambda2,     c2);
  
  // Half-Cauchy(0,1) prior on lambda_j + Jacobian for log-transform
  double log_prior_lambda_ratio =
  (log_half_cauchy_density(lambda_new, 1.0) + std::log(lambda_new))
    - (log_half_cauchy_density(lambda_cur, 1.0) + std::log(lambda_cur));
  
  double log_alpha = log_prior_theta_ratio + log_prior_lambda_ratio;
  
  int accepted = 0;
  
  if(log_alpha >= std::log(::unif_rand())){
    lambda2 = lambda2_prop;
    accept++;
    accepted = 1;
  }
  
  if(iter < burnin){
    double eta = 1.0 / std::sqrt(iter + 10.0);
    sd_log_lambda =
      std::exp(std::log(sd_log_lambda) + eta * (accepted - target_acc));
  }
  
  return List::create(
    _["lambda2"]        = lambda2,
    _["accept"]         = accept,
    _["sd_log_lambda"]  = sd_log_lambda,
    _["accepted"]       = accepted
  );
}

// ============================================================
// global tau update under regularized horseshoe
// proposal is on log(tau), target includes Jacobian +log(tau)
// tau ~ Half-Cauchy(0, scale_tau)
// ============================================================

// [[Rcpp::export]]
List exchange_move_tau_rhs_cpp_adapt(
    double tau2,
    const NumericVector& theta,
    const NumericVector& lambda2,
    double c2,
    double mu_alpha,
    double s2_alpha,
    double scale_tau,
    int accept,
    double sd_log_tau,
    int iter,
    int burnin,
    double target_acc = 0.25
){
  RNGScope rng;
  
  int J = theta.size();
  if(lambda2.size() != J) stop("theta and lambda2 length mismatch");
  if(tau2 <= 0.0) stop("tau2 must be positive");
  if(c2 <= 0.0)   stop("c2 must be positive");
  if(scale_tau <= 0.0) stop("scale_tau must be positive");
  
  double tau_cur = std::sqrt(tau2);
  if(tau_cur <= 0.0 || !R_finite(tau_cur))
    stop("current tau must be positive and finite");
  
  double log_tau_cur = std::log(tau_cur);
  double log_tau_new = log_tau_cur + R::rnorm(0.0, sd_log_tau);
  double tau_new     = std::exp(log_tau_new);
  double tau2_new    = tau_new * tau_new;
  
  if(!R_finite(tau_new) || tau_new <= 0.0) tau_new = 1e-12;
  
  double log_prior_theta_ratio =
    log_prior_theta_rhs_cpp(theta, mu_alpha, s2_alpha, tau2_new, lambda2, c2)
    - log_prior_theta_rhs_cpp(theta, mu_alpha, s2_alpha, tau2,     lambda2, c2);
  
  // Half-Cauchy(0, scale_tau) prior on tau + Jacobian for log-transform
  double log_prior_tau_ratio =
  (log_half_cauchy_density(tau_new, scale_tau) + std::log(tau_new))
    - (log_half_cauchy_density(tau_cur, scale_tau) + std::log(tau_cur));
  
  double log_alpha = log_prior_theta_ratio + log_prior_tau_ratio;
  
  int accepted = 0;
  
  if(log_alpha >= std::log(::unif_rand())){
    tau2 = tau2_new;
    accept++;
    accepted = 1;
  }
  
  if(iter < burnin){
    double eta = 1.0 / std::sqrt(iter + 10.0);
    sd_log_tau =
      std::exp(std::log(sd_log_tau) + eta * (accepted - target_acc));
  }
  
  return List::create(
    _["tau2"]        = tau2,
    _["accept"]      = accept,
    _["sd_log_tau"]  = sd_log_tau,
    _["accepted"]    = accepted
  );
}

#include <Rcpp.h>
using namespace Rcpp;

// ------------------------------------------------------------
// small helpers
// ------------------------------------------------------------
static inline double mean_cpp2(const NumericVector& x){
  int n = x.size();
  double s = 0.0;
  for(int i=0; i<n; ++i) s += x[i];
  return s / n;
}

static inline NumericVector center_cpp2(const NumericVector& x){
  int n = x.size();
  double m = mean_cpp2(x);
  NumericVector out(n);
  for(int i=0; i<n; ++i) out[i] = x[i] - m;
  return out;
}

static inline double log_norm_density2(double x, double mu, double s2){
  return -0.5 * std::log(2.0 * M_PI * s2) - 0.5 * (x - mu) * (x - mu) / s2;
}

static inline double log_half_cauchy_density2(double x, double scale){
  if(x <= 0.0 || scale <= 0.0) return R_NegInf;
  double z = x / scale;
  return std::log(2.0) - std::log(M_PI) - std::log(scale) - std::log1p(z * z);
}

static inline double rhs_lambda_tilde2_cpp(double lambda2, double tau2, double c2){
  return (c2 * lambda2) / (c2 + tau2 * lambda2);
}

// ------------------------------------------------------------
// These two MUST already exist in your current cpp codebase:
//   make_beta_roi_from_theta(theta)
//   gibbs_draw_k1_roi_shared_cpp(...)
//   piofx_roi_shared_cpp(...)
// ------------------------------------------------------------

// ============================================================
// VERSION 1: SIMPLE SPIKE-AND-SLAB
// u_j = theta_j - mean(theta)
// u_j | gamma_j ~ N(0, tau2_small) if gamma_j=0
//                 N(0, tau2_large) if gamma_j=1
// ============================================================

static inline double log_prior_theta_spike_slab_cpp(
    const NumericVector& theta,
    const IntegerVector& gamma,
    double tau2_small,
    double tau2_large,
    double mu_alpha,
    double s2_alpha
){
  int J = theta.size();
  if(gamma.size() != J) stop("gamma length mismatch");
  if(tau2_small <= 0.0 || tau2_large <= 0.0) stop("tau2_small/tau2_large must be positive");
  if(s2_alpha <= 0.0) stop("s2_alpha must be positive");
  
  double alpha = mean_cpp2(theta);
  NumericVector u = center_cpp2(theta);
  
  double lp = log_norm_density2(alpha, mu_alpha, s2_alpha);
  
  for(int j=0; j<J; ++j){
    double s2 = (gamma[j] == 1 ? tau2_large : tau2_small);
    lp += -0.5 * std::log(s2) - 0.5 * u[j] * u[j] / s2;
  }
  return lp;
}

// [[Rcpp::export]]
List exchange_move_theta_single_spike_slab_cpp_adapt(
    int k,
    NumericVector theta,
    int j_update_1based,
    const IntegerVector& gamma,
    double tau2_small,
    double tau2_large,
    double mu_alpha,
    double s2_alpha,
    const IntegerVector& group,
    const IntegerVector& roi,
    int n,
    int loop_aux,
    int accept,
    NumericMatrix wt,
    double sd_theta,
    int iter,
    int burnin,
    double target_acc = 0.25
){
  RNGScope rng;
  
  int J = theta.size();
  int j = j_update_1based - 1;
  if(j < 0 || j >= J) stop("j_update out of range");
  if(gamma.size() != J) stop("gamma length mismatch");
  
  NumericVector theta2 = clone(theta);
  theta2[j] = theta[j] + R::rnorm(0.0, sd_theta);
  
  NumericVector beta_roi  = make_beta_roi_from_theta(theta);
  NumericVector beta2_roi = make_beta_roi_from_theta(theta2);
  
  IntegerVector g2 =
    gibbs_draw_k1_roi_shared_cpp(
      k, n, beta2_roi, loop_aux, wt, roi, group
    );
  
  double log_like_ratio =
    piofx_roi_shared_cpp(group, beta2_roi, roi, wt)
    + piofx_roi_shared_cpp(g2,    beta_roi,  roi, wt)
    - piofx_roi_shared_cpp(group, beta_roi,  roi, wt)
    - piofx_roi_shared_cpp(g2,    beta2_roi, roi, wt);
    
    double log_prior_ratio =
    log_prior_theta_spike_slab_cpp(theta2, gamma, tau2_small, tau2_large, mu_alpha, s2_alpha)
      - log_prior_theta_spike_slab_cpp(theta,  gamma, tau2_small, tau2_large, mu_alpha, s2_alpha);
    
    double log_alpha = log_like_ratio + log_prior_ratio;
    
    int accepted = 0;
    if(log_alpha >= std::log(::unif_rand())){
      theta = theta2;
      accept++;
      accepted = 1;
    }
    
    if(iter < burnin){
      double eta = 1.0 / std::sqrt(iter + 10.0);
      sd_theta = std::exp(std::log(sd_theta) + eta * (accepted - target_acc));
    }
    
    return List::create(
      _["theta"]    = theta,
      _["accept"]   = accept,
      _["sd_theta"] = sd_theta,
      _["accepted"] = accepted
    );
}

// [[Rcpp::export]]
IntegerVector update_gamma_spike_slab_cpp(
    const NumericVector& theta,
    double tau2_small,
    double tau2_large,
    double pi
){
  RNGScope rng;
  
  int J = theta.size();
  if(tau2_small <= 0.0 || tau2_large <= 0.0) stop("tau2_small/tau2_large must be positive");
  if(pi <= 0.0 || pi >= 1.0) stop("pi must be in (0,1)");
  
  double alpha = mean_cpp2(theta);
  NumericVector u = center_cpp2(theta);
  
  IntegerVector gamma(J);
  
  for(int j=0; j<J; ++j){
    double log_p1 =
      -0.5 * std::log(tau2_large) - 0.5 * u[j] * u[j] / tau2_large + std::log(pi);
      double log_p0 =
      -0.5 * std::log(tau2_small) - 0.5 * u[j] * u[j] / tau2_small + std::log(1.0 - pi);
      
      double m = std::max(log_p1, log_p0);
      double p1 = std::exp(log_p1 - m);
      double p0 = std::exp(log_p0 - m);
      double prob1 = p1 / (p1 + p0);
      
      gamma[j] = (R::runif(0, 1) < prob1) ? 1 : 0;
  }
  
  return gamma;
}

