#' Compute constrained group Lasso 
#' 
#' tba
#' 
#' tba 
#'
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