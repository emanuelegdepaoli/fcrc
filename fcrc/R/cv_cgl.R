#' n-fold cross validation for a sparse concurrent functional regression model with compositional covariates
#' 
#' Perform n-fold cross validation for the model estimated by \link{alm_cgl_path} to
#' select the number of basis functions and \eqn{lambda}. 
#' 
#' @param Y A \eqn{n_t*p*n} array, containing \eqn{n} samples of \eqn{p} functional
#'          predictors observed at \eqn{n_t} points. 
#' @param Z A \eqn{n_t*n} array, containing \eqn{n} samples of the functional 
#'          response at \eqn{n_t} points. 
#' @param Zc optional. An array of \eqn{p_c} control variables, with same structure of Z. The
#'           first predictor must be a vector of ones. The default value is an array of ones,
#'           to estimate the intercept. 
#' @param L A matrix from \link{comp_L}.
#' @param t A vector containing the values at which functions are observed.
#' @param nfold optional. Number of folds. 
#' @param k_range optional. Number of basis functions.
#' @param ncores optional. Number of cores to be used. 
#' @param eps optional. Tolerance for ALM.
#' @param max_iter optional. Maximum number of iterations for ALM. 
#' @param max_iter_int optional. Maximum number of iterations for ADMM. 
#' @param abs_tol_int optional. Tolerance for primal feasibility condition of ADMM. 
#' @param rel_tol_int optional. Tolerance for dual feasibility condition of ADMM.
#' @param zero_tol optional. Tolerance for setting a coefficient to zero. 
#' @import parallel
#' @import foreach 
#' @importFrom doParallel registerDoParallel
#' @export
# n-fold cross validation 
# for parallelizing cpp code see: 
# https://stackoverflow.com/questions/6074310/using-rcpp-within-parallel-code-via-snow-to-make-a-cluster
cv_cgl = function(Y, Z, L, t, nfold = 10, k_range = 4:10, Zc = NULL, ncores = NULL,  
                  eps = 1e-5, max_iter_int = 1e3, max_iter = 100, 
                  abs_tol_int = 1e-6, rel_tol_int = 1e-4, zero_tol = 1e-4){
  
  n = dim(Z)[3]
  p = dim(Z)[2]
  n_k = length(k_range)
  n_t = length(t)
  if(is.null(ncores)) ncores = detectCores()-1
  nyears = dim(Z)[1]
  # if Zc is missing, use just a vector of one for the intercept
  if (is.null(Zc)) Zc = array(1, dim = c(n_t, 1, n))
  # randomly shuffle the data
  reord = sample(n)
  Z = Z[,, reord, drop = F]
  Zc = Zc[,, reord, drop = F]
  Y = Y[, reord, drop = F]
  # create 10 equally size folds
  folds = cut(seq(1, n), breaks = nfold, labels = FALSE)
  
  # fixed objects for each k
  # parallelization
  cl =  makeCluster(min(n_k, ncores))
  registerDoParallel(cl)
  res_data = foreach(j = 1:n_k, .packages = c("fda", "fcrc")) %dopar% {
    mat = mat_comp(Z, Y, t, k_range[j], Zc)
    # res raw data
    res = alm_cgl_path(matrices = mat, L = L, max_iter_int = max_iter_int, 
                       max_iter = max_iter, abs_tol_int = abs_tol_int, 
                       rel_tol_int = rel_tol_int, zero_tol = zero_tol)
    res
  }
  stopCluster(cl)
  lambda_path = lapply(res_data, function(x) exp(x$log_lambda_path))
  
  # perform n-fold cross validation 
  # parallelization
  outfile = "cv_outfile.txt"
  if (file.exists(outfile)) {
    # delete file if it exists
    file.remove(outfile)
  }
  cl =  makeCluster(min(n_k*nfold, ncores), outfile = outfile)
  registerDoParallel(cl)
  res_cv = foreach(i = 1:nfold, .packages = c("fda", "fcrc")) %:% 
    foreach(j = 1:n_k, .combine = "cbind") %dopar% {
      idx_test = which(folds == i, arr.ind = TRUE)
      Z_test = Z[,, idx_test, drop = F]
      Zc_test = Zc[,, idx_test, drop = F]
      Y_test = Y[, idx_test, drop = F]
      Z_train = Z[,, -idx_test, drop = F]
      Zc_train = Zc[,, -idx_test, drop = F]
      Y_train = Y[,-idx_test, drop = F]
      lambda_path_j = lambda_path[[j]]
      mat = mat_comp(Z_train, Y_train, t, k_range[j], Zc_train)
      Phi = mat$Phi
      res = alm_cgl_path(matrices = mat, L = L, lambda_path = lambda_path_j,  
                         max_iter_int = max_iter_int, max_iter = max_iter, 
                         abs_tol_int = abs_tol_int, rel_tol_int = rel_tol_int, 
                         zero_tol = zero_tol)
      apsse = sapply(1:100, function(l) error_comp(res$beta[l,], res$betac[l,], 
                                                   Y_test, Z_test, Phi, t, Zc_test)$APMSE*length(idx_test))
      print(paste("Fold", i, "k", k_range[j], "completed"))
      apsse
    }
  stopCluster(cl)
  
  # reorder results
  apsse = array(simplify2array(res_cv), dim = c(100, n_k, nfold))
  lambda_k = array(dim = c(100, n_k))
  for(j in 1:n_k) lambda_k[,j] = lambda_path[[j]]
  
  return(list(apsse = apsse, lambda_k = lambda_k, res_data = res_data,
              k_range = k_range, reord = reord, folds = folds))
}

