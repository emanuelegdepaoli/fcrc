################################################################################
# Simulations
################################################################################
rm(list = ls())
source('estimation_sim.R')
library(fda)
library(parallel)
library(foreach)
library(doParallel)

# simulation setup
sim_nrep = 1
# n, p, ncomp
temp1 = matrix(c(50, 40, 1, 
                 50, 40, 4, 
                 50, 100, 4), ncol = 3, byrow = T)
corr_t = c(0.2, 0.6)
corr_X = c(0.2, 0.6)
snr = c(2, 4)
temp2 = expand.grid(corr_t, corr_X, snr)
grid_sim = do.call(rbind, apply(temp1, 1, function(x) cbind(matrix(rep(x, nrow(temp2)), ncol = ncol(temp1), byrow = T),
                                                            temp2)))
colnames(grid_sim) = c("n", "p", "ncomp", "corr_t", "corr_x", "snr")

res = array(dim = c(sim_nrep, 9, 3, nrow(grid_sim)))
dimnames(res)[[2]] = c("pred_err", "est_error", "FPR", "FNR", "FPR_t", "FNR_t", 
                       "est_error_t", "pred_error_test", "k")
dimnames(res)[[3]] = c("CGL", "BGL", "naive")
attr(res, "info") = grid_sim
# estimation error
# mean(int_{0,1} abs(beta_est(t)-beta(t))^2 dt)^(0.5)

# idx setup to simulate
idx_setup_sim = 1:2

# fixed parameters
range = c(0, 1) # T range 
n_t = 20 # number of points in T range 
ncores = 7
# set.seed(111)
# t = sort(runif(n_t)) # random t
t = seq(0, 1, length.out = n_t) # fixed grid t
delta = diff(t)
k = 5 # number of basis for beta coefficients 
beta_basis = create.bspline.basis(rangeval = range, nbasis = k)
Phi = eval.basis(evalarg = t, basisobj = beta_basis) # nyears x k 
beta_list = list()
beta_list[[1]] = c(1, -1, 0, 0, 0, 
                   0, 0, -0.5, 1, 0, 
                   -1, 1, 0.5, -1, 0)
beta_list[[2]] = c(0.5, 0, 0, -0.5, 1,
                   0, 1, -1, 0, -1, 
                   -0.5, -1, 1, 0.5, 0)
beta_list[[3]] = c(0.5, -1, -1, 1, 0,
                   0, 1, 1, 0, 0,
                   -0.5, 0, 0, -1, 0)
beta_list[[4]] = c(1, 0, 0.5, 0, -1,
                   0, 0, -0.5, 0, 0,
                   -1, 0, 0, 0, 1)

# intercept
intercept = Phi%*%matrix(c(0, 0, 0, 0, 0), ncol = 1)

set.seed(222)
for (s in idx_setup_sim) {
  
  n = grid_sim[s, "n"]
  n_test = 1000
  p = grid_sim[s, "p"]
  n_comp = grid_sim[s, "ncomp"]
  p_vec = rep(p/n_comp, n_comp)
  n_null = p_vec-rep(3, n_comp)
  L = comp_L(p_vec)
  
  if(n_comp == 1) { # 1 composition with same sparsity as model with q compositions
    B = matrix(0, nrow = p, ncol = k)
    B[1:(3*4),] = matrix(unlist(beta_list), ncol = k, nrow = 3*4, byrow = T)
  } else {
    B_list = list()
    for(i in 1:n_comp) {
      idx_nonull = 1:(p_vec[i]-n_null[i])
      B_list[[i]] = array(0, dim = c(p_vec[i], k))
      B_list[[i]][idx_nonull,1:k] = matrix(beta_list[[i]], ncol = k, byrow = T) 
    }
    B = do.call(rbind, B_list)
  }
  idx_null = which(rowSums(abs(B)) == 0)
  Bphi = Phi%*%t(B)
  Beta_fdobj = fd(coef = t(B), basisobj = beta_basis)
  
  corr_t = grid_sim[s, "corr_t"]
  corr_x = grid_sim[s, "corr_x"]
  Sigma_t = corr_t^abs(outer(t,t,"-")) # autoregressive correlation structure
  
  # simulations for setup s start
  Y_test = array(dim = c(n_t, n_test))
  X_test = array(dim = c(n_t, p, n_test))
  
  # simulate covariates test set
  count = 1
  for(j in 1:n_comp){
    Sigma_x = diag(1, p_vec[j]) # compound symmetry correlation structure X 
    Sigma_x[Sigma_x != diag(Sigma_x)] = corr_x
    sigma2_x = 9
    Sigma_sim = sigma2_x*kronecker(Sigma_t, Sigma_x)
    L_Sigmasim = t(chol(Sigma_sim)) # Sigma_sim = L*t(L)
    for(i in 1:n_test){
      # w_i = (w_i(t_1), wi(t_2),..., wi(t_nt)) ~ N(0, Sigma_sim)
      # w_i(t_l) R^{p_vec[j]}
      wi = L_Sigmasim%*%matrix(rnorm(dim(Sigma_sim)[1]), dim(Sigma_sim)[1], 1)
      wij = matrix(wi, nrow = n_t, ncol = p_vec[j], byrow = T)
      Xij = apply(wij, 1, function(x) exp(x)/sum(exp(x)))
      X_test[,count:(count+p_vec[j]-1),i] = t(Xij)
    }
    count = count + p_vec[j]
  }
  Z_test = log(X_test)
  
  signal_test = array(dim = c(n_t, n_test))
  for(i in 1:n_test){
    signal_test[,i] = rowSums(Z_test[,1:p,i]*Bphi) + intercept
  }
  
  # SNR = trace(Cov(Y))/n*sigma^2_eps
  snr = grid_sim[s, "snr"]
  sigma_eps_test = sqrt(sum(diag(cov(t(signal_test))))/n_t/snr)
  
  # simulate Y test
  for(i in 1:n_test){
    Y_test[,i] = rowSums(Z_test[,1:p,i]*Bphi) + intercept + rnorm(n_t, 0, sigma_eps_test)
  }
  
  for (rep in 1:sim_nrep) {
    
    Y = array(dim = c(n_t, n))
    
    # simulate covariates 
    X = array(dim = c(n_t, p, n))
    count = 1
    for(j in 1:n_comp){
      Sigma_x = diag(1, p_vec[j]) # compound symmetry correlation structure X 
      Sigma_x[Sigma_x != diag(Sigma_x)] = corr_x
      sigma2_x = 9
      Sigma_sim = sigma2_x*kronecker(Sigma_t, Sigma_x)
      L_Sigmasim = t(chol(Sigma_sim)) # Sigma_sim = L*t(L)
      for(i in 1:n){
        # w_i = (w_i(t_1), wi(t_2),..., wi(t_nt)) ~ N(0, Sigma_sim)
        # w_i(t_l) R^{p_vec[j]}
        wi = L_Sigmasim%*%matrix(rnorm(dim(Sigma_sim)[1]), dim(Sigma_sim)[1], 1)
        wij = matrix(wi, nrow = n_t, ncol = p_vec[j], byrow = T)
        Xij = apply(wij, 1, function(x) exp(x)/sum(exp(x)))
        X[,count:(count+p_vec[j]-1),i] = t(Xij)
      }
      count = count + p_vec[j]
    }
    Z = log(X)
    
    signal = array(dim = c(n_t, n))
    for(i in 1:n){
      signal[,i] = rowSums(Z[,1:p,i]*Bphi) + intercept
    }
    
    # SNR = trace(Cov(Y))/n*sigma^2_eps
    snr = grid_sim[s, "snr"]
    sigma_eps = sqrt(sum(diag(cov(t(signal))))/n_t/snr)
    
    # simulate Y
    for(i in 1:n){
      Y[,i] = rowSums(Z[,1:p,i]*Bphi) + intercept + rnorm(n_t, 0, sigma_eps)
    }
    
    #####
    # CGL
    #####
    res_cv_cgl = cv_cgl(Y = Y, Z = Z, L = L, t = t, nfold = 10, k_range = 4:6, ncores = ncores, eps = 1e-6)
    obj_res_cv_cgl = plot_cv(res_cv_cgl)
    k_cgl = obj_res_cv_cgl$k_min1se
    res[rep, 9, 1, s] = k_cgl
    lambda_cgl = exp(obj_res_cv_cgl$min1se$log_lambda)
    mat_cgl = mat_comp(Z, Y, t, k_cgl)
    res_cgl = alm_cgl_path(mat_cgl, L, lambda_cgl)
    idx_null_cgl = which(res_cgl$index_null == 1)
    B_cgl_est = matrix(res_cgl$beta, ncol = k_cgl, byrow = T)
    
    res[rep, 1, 1, s] = error_comp(res_cgl$beta, res_cgl$betac, Y, Z, mat_cgl$Phi, t)$APMSE
    
    est_cgl_fdobj = fd(coef = t(B_cgl_est),
                       basisobj = create.bspline.basis(rangeval = range, nbasis = k_cgl))
    err_cgl_fdobj = est_cgl_fdobj-Beta_fdobj
    res[rep, 2, 1, s] = mean(sqrt(diag(inprod(err_cgl_fdobj, err_cgl_fdobj))))
    
    res[rep, 3, 1, s] = (1-mean(idx_null%in%idx_null_cgl))*100
    res[rep, 4, 1, s] = (mean((1:p)[-idx_null]%in%idx_null_cgl))*100
    
    # selection via threshold
    # (int_{0,1} abs(beta_est(t))^2 dt)^(0.5)
    energy = sqrt(diag(inprod(est_cgl_fdobj, est_cgl_fdobj)))
    idx_null_cgl = which(energy/sum(energy) <= 1/p)
    res[rep, 5, 1, s] = (1-mean(idx_null%in%idx_null_cgl))*100
    res[rep, 6, 1, s] =  (mean((1:p)[-idx_null]%in%idx_null_cgl))*100
    
    # est error after selection via threshold
    B_cgl_est[idx_null_cgl,] = 0
    est_cgl_fdobj = fd(coef = t(B_cgl_est),
                       basisobj = create.bspline.basis(rangeval = range, nbasis = k_cgl))
    err_cgl_fdobj = est_cgl_fdobj-Beta_fdobj
    res[rep, 7, 1, s] = mean(sqrt(diag(inprod(err_cgl_fdobj, err_cgl_fdobj))))
    
    # prediction test set
    res[rep, 8, 1, s] =  error_comp(res_cgl$beta, res_cgl$betac, Y_test, Z_test, mat_cgl$Phi, t)$APMSE
    
    #####
    # RGL
    #####
    idx_start = c(1, cumsum(head(p_vec, -1)) + 1)
    idx_end = cumsum(p_vec)
    r = sapply(1:n_comp, function(j) sample(idx_start[j]:idx_end[j], 1)) # random reference category
    idx_Zr = (1:p)[-r]
    idx_null_Zr = idx_null[idx_null%in%idx_Zr]
    idx_temp = unlist(lapply(1:n_comp, function(i) rep(r[i], p_vec[i])))
    Zr = array(dim = c(n_t, p-n_comp, n))
    for (i in 1:n_t) {
      for (l in 1:n) {
        zril = log(X[i,,l])-log(X[i,idx_temp,l])
        Zr[i,,l] = zril[zril!=0]
      }
    }
    Zr_test = array(dim = c(n_t, p-n_comp, n_test))
    for (i in 1:n_t) {
      for (l in 1:n_test) {
        zril = log(X_test[i,,l])-log(X_test[i,idx_temp,l])
        Zr_test[i,,l] = zril[zril!=0]
      }
    }
    
    res_cv_bgl = cv_bgl(Y = Y, Z = Zr, t = t, nfold = 10, k_range = 4:6, ncores = ncores)
    obj_res_cv_bgl = plot_cv(res_cv_bgl)
    k_bgl = obj_res_cv_bgl$k_min1se
    res[rep, 9, 2, s] = k_bgl
    lambda_bgl = exp(obj_res_cv_bgl$min1se$log_lambda)
    mat_bgl = mat_comp(Zr, Y, t, k_bgl)
    res_bgl = bgl_path(mat_bgl, lambda_bgl)
    idx_null_bgl = which(res_bgl$index_null == 1)
    B_bgl_est = matrix(0, ncol = k_bgl, nrow = p)
    B_bgl_est[idx_Zr,] = matrix(res_bgl$beta, ncol = k_bgl, byrow = T)
    B_bgl_est[r,] = -t(sapply(1:n_comp, function(j) colSums(B_bgl_est[(idx_start[j]:idx_end[j]),])))
    
    res[rep, 1, 2, s] = error_comp(res_bgl$beta , res_bgl$beta0, Y, Zr, mat_bgl$Phi, t)$APMSE
    
    est_bgl_fdobj = fd(coef = t(B_bgl_est),
                       basisobj = create.bspline.basis(rangeval = range, nbasis = k_bgl))
    err_bgl_fdobj = est_bgl_fdobj-Beta_fdobj
    res[rep, 2, 2, s] = mean(sqrt(diag(inprod(err_bgl_fdobj, err_bgl_fdobj))))
    
    res[rep, 3, 2, s] = (1-mean(idx_null%in%idx_null_bgl))*100
    res[rep, 4, 2, s] =  (mean((1:p)[-idx_null]%in%idx_null_bgl))*100
    
    # selection via threshold
    # (int_{0,1} abs(beta_est(t))^2 dt)^(0.5)
    energy = sqrt(diag(inprod(est_bgl_fdobj, est_bgl_fdobj)))
    idx_null_bgl = which(energy/sum(energy) <= 1/p)
    res[rep, 5, 2, s] = (1-mean(idx_null%in%idx_null_bgl))*100
    res[rep, 6, 2, s] =  (mean((1:p)[-idx_null]%in%idx_null_bgl))*100
    
    # est error after selection via threshold
    B_bgl_est[idx_null_bgl,] = 0
    est_bgl_fdobj = fd(coef = t(B_bgl_est),
                       basisobj = create.bspline.basis(rangeval = range, nbasis = k_bgl))
    err_bgl_fdobj = est_bgl_fdobj-Beta_fdobj
    res[rep, 7, 2, s] = mean(sqrt(diag(inprod(err_bgl_fdobj, err_bgl_fdobj))))
    
    # prediction test set
    res[rep, 8, 2, s] =  error_comp(res_bgl$beta , res_bgl$beta0, Y_test, Zr_test, mat_bgl$Phi, t)$APMSE
    
    
    #######
    # naive
    #######
    Beta_naive = array(dim = c(n_t, p))
    beta0_naive = rep(0, n_t)
    cv_full = cv_clt_full(Z = Z, Y = Y, L = L, t = t, nfold = 10, ncores = ncores)
    for(i in 1:n_t) {
      temp = alm_clt_path(Z = t(Z[i,,]), Y = Y[i,], L = L,
                     lambda = exp(cv_full$log_lambda_min1se[i]))
      Beta_naive[i,] =  temp$beta
      beta0_naive[i] = temp$beta0
    }
    # smoothing splines
    Beta_naive_smooth = apply(Beta_naive, 2,
                              function(x) smooth.spline(x, all.knots = T)$y)
    beta0_naive_smooth = smooth.spline(beta0_naive, all.knots = T)$y
    Beta_naive_smooth_coef = apply(Beta_naive, 2, function(x) smooth.spline(x, all.knots = T)$fit$coef)
    # smoothing splines basis to create fdobj
    splines_basis = create.bspline.basis(rangeval = c(0, 1), nbasis = n_t+2)
    beta_smooth_fdobj = fd(coef = Beta_naive_smooth_coef, basisobj = splines_basis)
    # (int_{0,1} abs(beta_est(t))^2 dt)^(0.5)
    energy = sqrt(diag(inprod(beta_smooth_fdobj, beta_smooth_fdobj)))
    idx_null_naive = which(energy/sum(energy) <= 1/p)
    beta_smooth_fdobj = fd(coef = Beta_naive_smooth_coef, basisobj = splines_basis)
    
    Y_est_naive = apply(Z, 3, function(x) rowSums(x*Beta_naive_smooth)) + beta0_naive
    res[rep, 1, 3, s] = mean((Y_est_naive-Y)^2)
    
    err_naive_fdobj = beta_smooth_fdobj-Beta_fdobj
    res[rep, 2, 3, s] = mean(sqrt(diag(inprod(err_naive_fdobj, err_naive_fdobj))))
    
    res[rep, 5, 3, s] = (1-mean(idx_null%in%idx_null_naive))*100
    res[rep, 6, 3, s] =  (mean((1:p)[-idx_null]%in%idx_null_naive))*100
    
    # est error after selection via threshold
    # set estimated null curves to zero
    Beta_naive_smooth[, idx_null_naive] = 0
    Beta_naive_smooth_coef[, idx_null_naive] = 0
    beta_smooth_fdobj = fd(coef = Beta_naive_smooth_coef, basisobj = splines_basis)
    err_naive_fdobj = beta_smooth_fdobj-Beta_fdobj
    res[rep, 7, 3, s] = mean(sqrt(diag(inprod(err_naive_fdobj, err_naive_fdobj))))
    
    Y_est_test_naive = apply(Z_test, 3, function(x) rowSums(x*Beta_naive_smooth)) + beta0_naive
    res[rep, 8, 3, s] = mean((Y_est_test_naive-Y_test)^2)
    
    print(paste("iter", rep, "setup", s))
    
  }
  
  # save partial result setup s
  temp = res[,,,s]
  attr(temp, "info") = grid_sim[s,]
  saveRDS(res[,,,s], paste0("sim_setup", s, ".rds"))
  # print summary setup s
  mean_res = apply(res[,1:8,,s], c(2, 3), mean)
  se_res = apply(res[,1:8,,s], c(2, 3), sd)/sqrt(sim_nrep)
  print("k_mode")
  print(apply(res[,9,,s], 2, function(x) names(sort(-table(x)))[1]))
  print("k_max")
  print(apply(res[,9,,s], 2, max))
  
  print(attr(temp, "info"))
  print(round(mean_res, 3))
  print(round(se_res, 3))
  
}

saveRDS(res, "sim_full_4comp.rds")

