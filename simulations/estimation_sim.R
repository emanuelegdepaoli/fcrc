################################################################################
# Sparse functional concurrent log-constrast regression
# estimation file for simulations
# model 1: constrained group LASSO as in estimation
# model 2: group LASSO with baseline level 
# model 3: log-contrast regression at each time t with LASSO regularization
# all the functions are implemented without control variables 
################################################################################
library(fcrc)
library(parallel)
library(foreach)

################################################################################
# Baseline group Lasso
################################################################################

# Regularization path for group Lasso with baseline level
path_comp_bgl = function(matrices){
  p = matrices$p
  k = matrices$k
  ncoef = p*k
  idx = rep(1:p, each = k)
  grps = unique(idx)
  n = matrices$n
  qtminvp = t(matrices$Q)%*%solve(matrices$M)%*%matrices$P
  xtystar = matrices$J - qtminvp
  lambda_max = max(sapply(grps, function(j) norm(xtystar[idx == j], "2")))
  # grid as in Regularization Paths for Generalized Linear Models
  # Friedman et al, 2010
  lambda_path = exp(seq(log(lambda_max), log(lambda_max*0.001),
                        length.out = 100))
  
  return(lambda_path)
}

# baseline group Lasso for a regularization path 
bgl_path = function(matrices, lambda_path = NULL, beta_start = NULL, u_start = NULL, 
                    rho = 1, varying_rho = T, abs_tol = 1e-6, rel_tol = 1e-4, 
                    max_iter = 1e4, zero_tol = 1e-4){
  
  # initialization
  p = matrices$p
  pc = matrices$pc
  k = matrices$k
  ncoef = p*k
  idx = rep(1:p, each = k)
  K = matrices$K
  J = matrices$J
  grps = unique(idx)
  M = matrices$M 
  P = matrices$P 
  Q = matrices$Q 
  Minv = solve(M)
  qtminvq = t(Q)%*%Minv%*%Q
  qtminvp = t(Q)%*%Minv%*%P
  # subtract intercept
  K = K - qtminvq
  J = J - qtminvp
  if (is.null(lambda_path)) {
    lambda_max = max(sapply(grps, function(j) norm(J[idx == j], "2")))
    # grid as in Regularization Paths for Generalized Linear Models
    # Friedman et al, 2010
    lambda_path = exp(seq(log(lambda_max), log(lambda_max*0.001),
                          length.out = 100))
  }
  nlambdas = length(lambda_path)
  if(is.null(beta_start)) beta_start = matrix(rep(0, ncoef), ncol = 1) 
  if(is.null(u_start)) u_start = matrix(rep(0, ncoef), ncol = 1) 

  res_beta0 = array(0, dim = c(nlambdas, pc*k)) 
  res_idxnull = array(0, dim = c(nlambdas, p))
  res_npred = rep(0, nlambdas)
  
  res = admm_grplasso_path_int(xty = J, xtx = K, index = idx, lambda_path = lambda_path,
                               beta_start = beta_start, u_start = u_start, 
                               varying_rho = varying_rho, rho = rho, abs_tol = abs_tol,
                               rel_tol = rel_tol, max_iter = max_iter)
  
  res_beta = t(res$res_beta)
  
  for(i in 1:nlambdas){
    # index null coef
    idx_null = which(sapply(1:p, function(j) max(abs(res_beta[i,idx == j])) < zero_tol))
    res_idxnull[i,idx_null] = 1
    # intercept and control variables 
    res_beta0[i,] = matrix(Minv%*%(P-Q%*%res_beta[i,]), nrow = 1)
    res_npred[i] = p-length(idx_null)
  }
  
  return(list(beta = res_beta, beta0 = res_beta0, idx_null = res_idxnull, 
              npred = res_npred, log_lambda_path = log(lambda_path)))
}

# n-fold cross validation 
# for parallelizing cpp code see: 
# https://stackoverflow.com/questions/6074310/using-rcpp-within-parallel-code-via-snow-to-make-a-cluster
cv_bgl = function(Y, Z, t, nfold = 10, k_range = 4:10, ncores = NULL,
                  max_iter = 1e4, abs_tol = 1e-6, rel_tol = 1e-4, zero_tol = 1e-4){
  
  n = dim(Z)[3]
  p = dim(Z)[2]
  n_k = length(k_range)
  n_t = length(t)
  if(is.null(ncores)) ncores = detectCores()-1
  nyears = dim(Z)[1]
  Zc = array(1, dim = c(n_t, 1, n))
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
  # call functions in each cluster 
  clusterEvalQ(cl, {source("estimation_sim.R", local = T)})
  res_data = foreach(j = 1:n_k, .packages = c("fda", "fcrc")) %dopar% {
    # regularization path
    mat = mat_comp(Z, Y, t, k_range[j], Zc)
    res = bgl_path(matrices = mat, max_iter = max_iter, abs_tol = abs_tol, 
                   rel_tol = rel_tol, zero_tol = zero_tol)
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
  cl =  makeCluster(min(nfold*n_k, ncores), outfile = outfile)
  registerDoParallel(cl)
  # call functions in each cluster 
  clusterEvalQ(cl, {source("estimation_sim.R", local = T)})
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
      res = bgl_path(matrices = mat, lambda_path = lambda_path_j, 
                     max_iter = max_iter, abs_tol = abs_tol, rel_tol = rel_tol, 
                     zero_tol = zero_tol)
      apsse = sapply(1:100, function(l) error_comp(res$beta[l,], res$beta0[l,], 
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

################################################################################
# log-contrast model at each time t with LASSO penalty and arbitrary number of 
# compositions
################################################################################

# Augmented Lagrangian Method (ALM) for log-contrast model at each time t with 
# LASSO penalty and arbitrary number of compositions

# Regularization path 
path_comp_clt = function(Z, Y){
  p = dim(Z)[2]
  idx = grps = 1:p
  n = dim(Z)[1]
  
  # center predictors and response
  Yc = Y-mean(Y)
  Zmeans = colMeans(Z)
  Zc = t(apply(Z, 1, function(x) x-Zmeans))
  
  K = t(Zc)%*%Zc
  J = t(Zc)%*%Y
  lambda_max = max(sapply(grps, function(j) norm(J[idx == j], "2")))
  # grid as in Regularization Paths for Generalized Linear Models
  # Friedman et al, 2010
  lambda_path = exp(seq(log(lambda_max), log(lambda_max*0.001),
                        length.out = 100))
  
  return(lambda_path)
}

# ALM for a regulariation path
alm_clt_path = function(Z, Y, L, lambda_path = NULL, beta_start = NULL, u_start = NULL, 
                        rho = 1, max_rho = 1e6, gamma = 1, varying_gamma = T, eps = 1e-5, 
                        max_iter_int = 1e3, max_iter = 100, abs_tol_int = 1e-6, 
                        rel_tol_int = 1e-4, zero_tol = 1e-4){
  
  # initialization
  n = dim(Z)[1]
  p = dim(Z)[2]
  # center predictors and response
  Yc = Y-mean(Y)
  Zmeans = colMeans(Z)
  Zc = t(apply(Z, 1, function(x) x-Zmeans))
  
  K = t(Zc)%*%Zc
  J = t(Zc)%*%Yc
  grps = idx = 1:p
  if (is.null(beta_start)) beta = matrix(rep(0, p), ncol = 1) else beta = beta_start
  if (is.null(u_start)) u = matrix(rep(0, dim(L)[1]), ncol = 1) else u = u_start
  u_int = matrix(rep(0, p), ncol = 1)
  if (is.null(lambda_path)) lambda_path = path_comp_clt(Z, Y)
  nlambdas = length(lambda_path)
  res_beta = array(0, dim = c(nlambdas, p))
  res_beta0 = rep(0, nlambdas)
  res_idxnull = array(0, dim = c(nlambdas, p))
  res_npred = rep(0, nlambdas)
  
  res = alm_cgl_path_int(J, K, L, idx, u, u_int, beta, lambda_path, rho, max_rho, 
                         gamma, varying_gamma, eps, max_iter, max_iter_int, abs_tol_int,
                         rel_tol_int)
  
  res_beta = t(res$res_beta)
  
  for(i in 1:nlambdas){
    # index null coef
    idx_null = which(sapply(1:p, function(j) max(abs(res_beta[i,idx == j])) < zero_tol))
    res_idxnull[i,idx_null] = 1
    # intercept and control variables 
    res_beta0[i] = mean(Y)-sum(Zmeans*res_beta[i,])
    res_npred[i] = p-length(idx_null)
  }
  
  return(list(beta = res_beta, beta0 = res_beta0, idx_null = res_idxnull, 
              npred = res_npred, log_lambda_path = log(lambda_path), 
              niter = res$res_niter))
}


# n-fold cross validation 
cv_clt_full = function(Z, Y, L, t, nfold = 10, ncores = NULL, eps = 1e-5, 
                       max_iter_int = 1e3, max_iter = 100, abs_tol_int = 1e-6, 
                       rel_tol_int = 1e-4, zero_tol = 1e-4){
  
  n_t = length(t)
  n = dim(Z)[3]
  p = dim(Z)[2]
  if(is.null(ncores)) ncores = detectCores()-1
  # randomly shuffle the data for each t
  reord = lapply(1:n_t, function(t) sample(n))
  # create 10 equally size folds
  folds = cut(seq(1, n), breaks = nfold, labels = FALSE)

  lambda_path = lapply(1:n_t, function(t) path_comp_clt(t(Z[t,,]), Y[t,]))
  
  # perform n-fold cross validation for each t
  # parallelization
  outfile = "cv_outfile.txt"
  if (file.exists(outfile)) {
    # delete file if it exists
    file.remove(outfile)
  }
  cl =  makeCluster(min(ncores, nfold*n_t), outfile = outfile) 
  registerDoParallel(cl)
  # call functions in each cluster 
  clusterEvalQ(cl, {source("estimation_sim.R", local = T)})
  res_cv = foreach(i = 1:nfold) %:% 
    foreach(j = 1:n_t, .combine = "cbind") %dopar% {
      idx_test = reord[[j]][folds == i]
      idx_train = reord[[j]][folds != i]
      Z_test = t(Z[j,,idx_test])
      Y_test = matrix(Y[j,idx_test], ncol = 1)
      Z_train = t(Z[j,,idx_train])
      Y_train = matrix(Y[j,idx_train], ncol = 1)
      res = alm_clt_path(Z_train, Y_train, L, lambda_path[[j]], eps = eps, max_iter_int = max_iter_int, 
                     max_iter = max_iter, abs_tol_int = abs_tol_int, 
                     rel_tol_int = rel_tol_int, zero_tol = zero_tol)
      
      Y_est = sweep(Z_test%*%t(res$beta), 2, res$beta0, "+")
      apsse = colSums(sweep(Y_est, 1, Y_test, "-")^2)
      print(paste("Fold", i, "t", (1:n_t)[j], "completed"))
      apsse 
    }
  stopCluster(cl)
  
  apmse = array(0, dim = c(100, n_t))
  idx_min = rep(0, n_t)
  idx_min1se = rep(0, n_t)
  log_lambda_min1se =  rep(0, n_t)
  log_lambda_min =   rep(0, n_t)
  for(j in 1:n_t) {
    apsse_j = array(0, dim = c(100, nfold))
    for(i in 1:nfold) {
      apsse_j[,i] = res_cv[[i]][,j]
    }
    apmse[,j] = rowSums(apsse_j)/n
    idx_min[j] = which(apmse[,j] == min(apmse[,j]))
    min = apmse[idx_min[j],j]
    apmse_j2 = sapply(1:nfold, function(i) apsse_j[,i]/sum(folds == i)) # mean of apmse of each folder for se
    cv_se = apply(apmse_j2, 1, sd)/sqrt(nfold)
    idx_min1se[j] = which(apmse[,j] < min+cv_se[idx_min[j]])[1]
    log_lambda_min1se[j] = log(lambda_path[[j]][idx_min1se[j]])
    log_lambda_min[j] =  log(lambda_path[[j]][idx_min[j]])
    
  }
  
  return(list(apmse = apmse, log_lambda_path = lapply(lambda_path, log), 
              log_lambda_min1se = log_lambda_min1se, log_lambda_min = log_lambda_min,
              reord = reord, folds = folds))
}


