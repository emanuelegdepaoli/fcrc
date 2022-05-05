################################################################################
# Sparse functional concurrent log-constrast regression
# final script: no type (no Lasso, no penalty Mayer)
# added number of predictors alm_path 
################################################################################

# Least squares with constraint L*beta=c
# Lagrangian L(beta, u) = -t(y)*x*beta+t(beta)*t(x)*x*beta+t(u)*(Ltilda*beta-c)
# functional regression
# CONTROL VARIABLES TO BE ADDED 
fcls = function(matrices, Ltilda, c){
  k = matrices$k
  n = matrices$n
  p = matrices$p
  ncoef = p*k
  nconstr = nrow(Ltilda)
  J = matrices$J
  K = matrices$K
  M = matrices$M
  P = matrices$P
  Q = matrices$Q
  xtx = K - t(Q)%*%solve(M)%*%Q
  mat2inv = cbind(rbind(xtx, Ltilda), 
                  rbind(t(Ltilda), matrix(0, ncol = nconstr, nrow = nconstr))
  )
  xtys = rbind(J - t(Q)%*%solve(M)%*%P, 
               matrix(c, ncol = 1, nrow = nconstr))
  sol = solve(mat2inv)%*%xtys
  beta = sol[1:ncoef]
  # intercept and control variables 
  betac = (solve(M)%*%(P-Q%*%beta))
  u = sol[-(1:ncoef)]
  
  return(list(beta = beta, beta0 = betac[1:k], betac = betac, u = u))
} 

# ordinary regression
# CONTROL VARIABLES TO BE ADDED 
cls = function(y, Z, L, c, intercept = T){
  n = nrow(Z)
  p = ncol(Z)
  nconstr = nrow(L)
  ymean = mean(y)
  Zmeans = colMeans(Z)
  if (intercept) {
    yc = y-ymean
    Zc = t(apply(Z, 1, function(x) x-Zmeans))
  } else {
    yc = y
    Zc = Z
  }
  mat2inv = cbind(rbind(t(Zc)%*%Zc, L), 
                  rbind(t(L), matrix(0, ncol = nconstr, nrow = nconstr))
  )
  xtys = rbind(t(Zc)%*%yc, matrix(c, ncol = 1, nrow = nconstr))
  sol = solve(mat2inv)%*%xtys
  beta = sol[1:p]
  if (intercept) beta0 = ymean-sum(Zmeans*beta) else beta0 = NULL
  u = sol[-(1:p)]
  
  # SSE and F test H0:beta_j=0 forall j
  res = y - beta0 - Z%*%beta
  SSE = sum(res^2)
  SSE0 = sum((y-ymean)^2)
  Ftest_null = (SSE0-SSE)/p*(n-p-1)/SSE
  if (intercept) {
    pval = pf(Ftest_null, p, n-p-1, lower.tail = F)
  } else {
    pval = NULL
  }
  
  return(list(beta = beta, beta0 = beta0, u = u, residuals = res, pval_Ftest0 = pval))
} 

alm_cgl = function(matrices, L, lambda, beta_start = NULL, u_start = NULL, rho = 1, 
                   max_rho = 1e6, gamma = 1, varying_gamma = T, eps = 1e-5, max_iter_int = 1e3, 
                   max_iter = 100, abs_tol_int = 1e-6, rel_tol_int = 1e-4, 
                   zero_tol = 1e-4){
  
  # initialization
  k = matrices$k
  n = matrices$n
  p = matrices$p
  ncoef = p*k
  idx = rep(1:p, each = k)
  Ltilda = kronecker(L, diag(1, k))
  K = matrices$K
  J = matrices$J
  grps = unique(idx)
  if (is.null(beta_start)) beta = matrix(rep(0, ncoef), ncol = 1) else beta = beta_start
  if (is.null(u_start)) u = matrix(rep(0, k*dim(L)[1]), ncol = 1) else u = u_start
  u_int = matrix(rep(0, ncoef), ncol = 1)
  err = array(0, max_iter)
  rho_list = array(0, max_iter)
  rho_old = rho
  niter_int = array(0, max_iter)
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
  
  # iterations
  res = alm_cgl_int(J, K, Ltilda, idx, u, u_int, beta, lambda, rho, max_rho, 
                    gamma, varying_gamma, eps, max_iter, max_iter_int, abs_tol_int,
                    rel_tol_int)
  
  # index null coef
  index_null = which(sapply(1:p, function(j) max(abs(res$beta[idx == j])) < zero_tol))
  
  # intercept and control variables 
  betac = (Minv%*%(P-Q%*%res$beta))
  
  return(list(beta = res$beta, beta0 = betac[1:k], betac = betac, u = res$u, 
              index_null = index_null, rho = res$rho_list, 
              err = res$err, niter = res$iter, niter_int = res$niter_int, 
              gamma_list = res$gamma_hist))
}

# Regularization path
path_comp_cgl = function(matrices){
  p = matrices$p
  k = matrices$k
  ncoef = p*k
  idx = rep(1:p, each = k)
  grps = unique(idx)
  n = matrices$n
  # scale K and J
  qtminvp = t(matrices$Q)%*%solve(matrices$M)%*%matrices$P
  xtystar = matrices$J - qtminvp
  lambda_max = max(sapply(grps, function(j) norm(xtystar[idx == j], "2")))
  # grid as in Regularization Paths for Generalized Linear Models
  # Friedman et al, 2010
  lambda_path = exp(seq(log(lambda_max), log(lambda_max*0.001),
                        length.out = 100))
  
  return(lambda_path)
}

# Augmented Lagrangian Method for a regularization path 
alm_cgl_path = function(matrices, L, lambda_path = NULL, beta_start = NULL, u_start = NULL, 
                        rho = 1, max_rho = 1e6, gamma = 1, varying_gamma = T, eps = 1e-5, 
                        max_iter_int = 1e3, max_iter = 100, abs_tol_int = 1e-6, rel_tol_int = 1e-4, 
                        zero_tol = 1e-4){
  p = matrices$p
  pc = matrices$pc
  k = matrices$k
  ncoef = p*k
  idx = rep(1:p, each = k)
  if (is.null(lambda_path)) lambda_path = path_comp_cgl(matrices)
  if (is.null(beta_start)) beta_start = matrix(rep(0, ncoef), ncol = 1) 
  if (is.null(u_start)) u_start = matrix(rep(0, k*dim(L)[1]), ncol = 1) 
  res_beta = array(0, dim = c(100, ncoef))
  res_betac = array(0, dim = c(100, pc*k)) 
  res_idxnull = array(0, dim = c(100, p))
  res_npred = rep(0, 100)
  
  for (i in 1:100) {
    obj =  alm_cgl(matrices, L, lambda_path[i], beta_start, u_start, rho, max_rho,
                   gamma, varying_gamma, eps, max_iter_int, max_iter, abs_tol_int, 
                   rel_tol_int, zero_tol)
    res_betac[i,] = obj$betac
    res_beta[i,] = obj$beta
    beta_start = obj$beta # warm start
    u_start = obj$u # warm start
    res_idxnull[i, obj$index_null] = 1
    res_npred[i] = p-length(obj$index_null)
  }
  
  return(list(beta = res_beta, betac = res_betac, idx_null = res_idxnull, 
              npred = res_npred, log_lambda_path = log(lambda_path)))
}

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
  cl =  makeCluster(ncores, outfile = outfile)
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

