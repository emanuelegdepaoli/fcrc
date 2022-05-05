# Alternative Direction Method of Multipliers (ADMM) for Group Lasso
# inspired by https://stanford.edu/~boyd/papers/admm/lasso/lasso.html
# and related paper Boyd, Stephen, Neal Parikh, and Eric Chu. 
# Distributed optimization and statistical learning via the alternating direction 
# method of multipliers. 
# Now Publishers Inc, 2011.
# scaled augmented lagragian implementation
# rho can vary for max. 50 iterations and then at least 50 more iterations
# are done
# loss function: 1/2*t(theta)*xtx*theta - t(theta)*xty + lambda*sum(||theta_j||_2)
# s.t theta_j = beta_j

# last edit: 19/01/2022
# added max(abs(beta[index == j])) <= 1e-10 as criterion for checking if beta_j = 0 in KKT conditions
# (on theta, the thresholded solution of ADMM)

################################################################################
# auxiliary functions
################################################################################

# soft-thresholding operator
soft_thre = function(x, par) x*max(0, 1-par/norm(x, '2'))

# objective function group lasso
obj_fun_grplasso = function(xty, xtx, beta, index, lambda){
  grps = unique(index)
  beta_norms = sapply(grps, function(j) norm(beta[index == j], '2'))
  0.5*t(beta)%*%xtx%*%beta-sum(beta*xty)+lambda*sum(beta_norms)
}

# function to check KKT conditions for group lasso
KKTcheck_grplasso = function(xty, xtx, index, beta, lambda, 
                             tol = 1e-2, print_warnings = T){
  grps = unique(index)
  idx_null = NULL
  warnings = rep(F, length(grps))
  
  for(j in grps){
    norm = norm(beta[index == j], "2")
    pj = length(beta[index == j])
    # beta_j=0
    # KKT beta_j=0 ||-x_j(y-Xbeta)||_2 <= lambda
    # beta_j!=0
    # KKT beta_j!=0 ||-x_j(y-Xbeta)||_2 = lambda
    if(max(abs(beta[index == j])) <= 1e-10){ 
      idx_null = c(idx_null, j)
      z = norm(-xty[index == j]+xtx[index == j,]%*%beta, '2')-lambda
      if(z > sqrt(pj)*tol) warnings[j] = T
    } else{
      z = mean(abs(-xty[index == j]+xtx[index == j,]%*%beta+
                   lambda*beta[index == j]/norm))
      if(z > tol) warnings[j] = T 
    }
    if(warnings[j] & print_warnings) warning(paste('KKT condition not satisfied', 'idx =', j, 
                                  'norm =', round(norm, 3), 'z = ', round(z, 3)))
  }
  
  return(list(warnings = warnings, index_null = idx_null))
}

################################################################################
# ADMM for group Lasso - cpp
################################################################################

admm_grplasso_cpp = function(xty, xtx, index, lambda, beta_start = NULL, u_start = NULL,
                         varying_rho = T, rho = 1, abs_tol = 1e-4, rel_tol = 1e-3,
                         KKT_tol = 1e-2, max_iter = 1e3, print_warnings = T){
  # initialization optional arguments
  p = ncol(xtx)
  if(is.null(beta_start)) beta_start = matrix(rep(0, p), ncol = 1) 
  if(is.null(u_start)) u_start = matrix(rep(0, p), ncol = 1) 
  
  res = admm_grplasso_int(xty, xtx, index, lambda, beta_start, u_start, varying_rho,
                          rho, abs_tol, rel_tol, max_iter)
  #index_null = which(sapply(unique(index), function(j) mean(abs(res$beta[index == j])) < 1e-6))
  index_null = KKTcheck_grplasso(xty, xtx, index, res$beta, lambda, tol = KKT_tol,
                                 print_warnings = print_warnings)$index_null
  res$index_null = index_null
  
  return(res)
}
