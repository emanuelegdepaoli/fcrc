#' Fit a sparse concurrent functional regression model with compositional covariates
#' 
#' Compute the estimates of the model 
#' \deqn{y(t) = Z_c(t)*beta_c(t) + Z(t)*beta(t) + e(t)}
#' subject to 
#' \deqn{L*beta(t) = 0}
#' with a Group Lasso penalty, as described in the paper, for a path of \eqn{lambdas}.
#' 
#' The model is estimated using the augmented Lagrangian method (ALM), with associated 
#' Lagrangian
#' \deqn{L_ρ(beta, u) = -t(beta)*J + 0.5*t(beta)*K*beta + lambda*sum(||theta_j||_2) 
#'  + t(u)*L*beta + rho/2*(||L*beta||_2)^2},
#' as described in the paper.
#'  
#' The Alternating Direction Method of Multipliers (ADMM) is used to solve the first 
#' step of ALM, which is a standard Group Lasso problem. The implementation follows Boyd et al. (2010) 
#' \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}. The penalty 
#' parameter of ADMM can vary for a maximum of 50 iterations maximum and then at 
#' least 50 more iterations are done.
#' @param matrices A list from \link{mat_comp}.
#' @param L A matrix from \link{comp_L}.
#' @param lambda_path The path of \eqn{lambda}s for which compute the solution. The default
#'                    value is the output of \link{path_comp_cgl}.
#' @param beta_start optional. A starting point for \eqn{beta}. Default value is zero. 
#' @param u_start optional. A starting point for \eqn{u}. Default value is zero. 
#' @param rho optional. The penalty parameter for ALM. 
#' @param max_rho optional. The upper bound for \eqn{rho}. 
#' @param varying_gamma logical, optional. Should the penalty parameter for internal ADMM vary at each iteration?. 
#'                      If FALSE, it varies only at the first internal ADMM iteration. 
#' @param gamma optional. The penalty parameter for ADMM.
#' @param eps optional. Tolerance for ALM.
#' @param abs_tol_int optional. Tolerance for primal feasibility condition of ADMM. 
#' @param rel_tol_int optional. Tolerance for dual feasibility condition of ADMM.
#' @param max_iter optional. Maximum number of iterations for ALM. 
#' @param max_iter_int optional. Maximum number of iterations for ADMM. 
#' @param zero_tol optional. Tolerance for setting a coefficient to zero. 
#' @return A list containing
#'         \enumerate{
#'            \item \eqn{beta} A matrix whose rows are \eqn{beta} for each \eqn{lambda}.
#'            \item \eqn{betac} A matrix whose rows are \eqn{beta_c} for each \eqn{lambda}.
#'            \item \emph{index_null} A matrix whose rows are a vector of (0,1) for each predictor,
#'                  one row for each \eqn{lambda}. If 1, the predictor is selected. 
#'            \item \emph{npred} Total number of predictors selected for each \eqn{lambda}.
#'            \item \emph{log_lambda_path} The path of \eqn{lambda} in log scale. 
#'            \item \emph{niter} Number of iterations performed by ALM for each \eqn{lambda}.
#'         }
#' @export
alm_cgl_path = function(matrices, L, lambda_path = NULL, beta_start = NULL, u_start = NULL, 
                        rho = 1, max_rho = 1e6, gamma = 1, varying_gamma = T, eps = 1e-5, 
                        max_iter_int = 1e3, max_iter = 100, abs_tol_int = 1e-6, rel_tol_int = 1e-4, 
                        zero_tol = 1e-4){
  
  # initialization
  k = matrices$k
  n = matrices$n
  p = matrices$p
  pc = matrices$pc
  ncoef = p*k
  idx = rep(1:p, each = k)
  Ltilda = kronecker(L, diag(1, k))
  K = matrices$K
  J = matrices$J
  grps = unique(idx)
  if (is.null(beta_start)) beta = matrix(rep(0, ncoef), ncol = 1) else beta = beta_start
  if (is.null(u_start)) u = matrix(rep(0, k*dim(L)[1]), ncol = 1) else u = u_start
  if (is.null(lambda_path)) lambda_path = path_comp_cgl(matrices)
  nlambdas = length(lambda_path)
  u_int = matrix(rep(0, ncoef), ncol = 1)
  M = matrices$M 
  P = matrices$P 
  Q = matrices$Q 
  Minv = solve(M)
  qtminvq = t(Q)%*%Minv%*%Q
  qtminvp = t(Q)%*%Minv%*%P
  # subtract control variables 
  K = K - qtminvq
  J = J - qtminvp
  xtxr = K + rho*t(Ltilda)%*%Ltilda
  
  res_betac = array(0, dim = c(nlambdas, pc*k)) 
  res_idxnull = array(0, dim = c(nlambdas, p))
  res_npred = rep(0, nlambdas)
  res = .Call(`_fcrc_alm_cgl_path_int`, J, K, Ltilda, idx, u, u_int, beta, lambda_path, 
              rho, max_rho, gamma, varying_gamma, eps, max_iter, max_iter_int, 
              abs_tol_int, rel_tol_int)
  res_beta = t(res$res_beta)
  
  for(i in 1:nlambdas){
    # index null coef
    idx_null = which(sapply(1:p, function(j) max(abs(res_beta[i,idx == j])) < zero_tol))
    res_idxnull[i,idx_null] = 1
    # intercept and control variables 
    res_betac[i,] = matrix(Minv%*%(P-Q%*%res_beta[i,]), nrow = 1)
    res_npred[i] = p-length(idx_null)
  }
  
  return(list(beta = res_beta, betac = res_betac, index_null = res_idxnull, 
              npred = res_npred, log_lambda_path = log(lambda_path), 
              niter = res$res_niter))
}