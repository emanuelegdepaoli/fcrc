################################################################################
# auxiliary functions
# last edit: 04/05/2022
# NULL instead of missing for optional arguments 
# added functions to plot coef by composition 
# fixed min1se CV 
################################################################################

# integral approximation: trapezoidal method
trapz = function(x, y){
  l = length(x)
  xs = x[-1]-x[-l]
  ys = 0.5*(y[-l]+y[-1])
  return(sum(xs*ys))
}

# constraint matrix
comp_L = function(p_vec){
  p = sum(p_vec)
  n_constr = length(p_vec)
  L = array(0, dim = c(n_constr, p))
  count = 1
  for(j in 1:n_constr){
    L[j, count:(count+p_vec[j]-1)] = 1
    count = count + p_vec[j]
  }
  return(L)
}

# matrices for computations
mat_comp = function(Z, Y, t, k, Zc = NULL){
  beta_basis = create.bspline.basis(rangeval = c(min(t), max(t)), nbasis = k, 
                                    norder = 4) # cubic-splines 
  Phi = eval.basis(evalarg = t, basisobj = beta_basis)
  
  # initialization
  p = dim(Z)[2]
  n = dim(Z)[3]
  n_t = dim(Z)[1]
  # if Zc is missing, use just a vector of one for the intercept
  if (is.null(Zc)) Zc = array(1, dim = c(n_t, 1, n))
  pc = dim(Zc)[2] # number of control variables, including the intercept
  J = array(dim = c(p*k, 1))
  K = array(dim = c(p*k, p*k))
  M = array(dim = c(pc*k, pc*k))
  R = array(dim = c(k, k))
  Q = array(dim = c(pc*k, p*k))
  P = array(dim = c(pc*k, 1))
  
  # computation J, K
  count = 1
  for(i in 1:p){
    # z_i(t)*y(t)
    zyt = rowSums(Z[,i,]*Y)
    for(l in 1:k){
      # int phi_l(t)*z_i(t)*y(t)
      J[count] = trapz(t, zyt*Phi[,l])
      count1 = 1
      for(j in 1:p){
        # z_i(t)*z_j(t)
        zizjt = rowSums(Z[,i,]*Z[,j,])
        for(m in 1:k){
          # int phi_l(t)*z_i(t)*z_j(t)*phi_m(t)
          K[count, count1] = trapz(t, Phi[,l]*zizjt*Phi[,m]) 
          count1 = count1+1
        }
      }
      count = count+1
    }
  }
  
  # computation M, P
  count = 1
  for(i in 1:pc){
    # z_i(t)*y(t)
    zyt = rowSums(Zc[,i,]*Y)
    for(l in 1:k){
      # int phi_l(t)*zc_i(t)*y(t)
      P[count] = trapz(t, zyt*Phi[,l])
      count1 = 1
      for(j in 1:pc){
        # zc_i(t)*z_j(t)
        zizjt = rowSums(Zc[,i,]*Zc[,j,])
        for(m in 1:k){
          # int phi_l(t)*zc_i(t)*zc_j(t)*phi_m(t)
          M[count, count1] = trapz(t, Phi[,l]*zizjt*Phi[,m]) 
          count1 = count1+1
        }
      }
      count = count+1
    }
  }
  
  # computation Q
  count = 1
  for(i in 1:pc){
    # z_i(t)*y(t)
    zyt = rowSums(Zc[,i,]*Y)
    for(l in 1:k){
      count1 = 1
      for(j in 1:p){
        # zc_i(t)*z_j(t)
        zizjt = rowSums(Zc[,i,]*Z[,j,])
        for(m in 1:k){
          # int phi_l(t)*z_i(t)*z_j(t)*phi_m(t)
          Q[count, count1] = trapz(t, Phi[,l]*zizjt*Phi[,m]) 
          count1 = count1+1
        }
      }
      count = count+1
    }
  }
  
  # computation R int phi(t)*phi(t)^t dt 
  for (i in 1:k) {
    for (j in 1:k) {
      R[i,j] = trapz(t, Phi[,i]*Phi[,j]) 
    }
  }
  
  # int d2 phi(t) * d2 phi(t)
  S = inprod(beta_basis, beta_basis, Lfdobj1 = 2, Lfdobj2 = 2)
  
  return(list(J = J, K = K, M = M, R = R, S = S, P = P, Q = Q, 
              p = p, pc = pc, n = n, k = k, Phi = Phi, beta_basis = beta_basis))
}

# MSE and prediction error
error_comp = function(beta, betac, Y, Z, Phi, t, Zc = NULL){
  n = dim(Z)[3]
  p = dim(Z)[2]
  k = ncol(Phi)
  n_t = length(t)
  if(is.null(Zc)) Zc = array(1, dim = c(n_t, 1, n))
  B_est = matrix(beta, ncol = k, byrow = T)
  Bc_est = matrix(betac, ncol = k, byrow = T)
  Y_est = apply(Z, 3, function(x) rowSums(x*Phi%*%t(B_est))) + 
    apply(Zc, 3, function(x) rowSums(x*Phi%*%t(Bc_est)))
  sse = sum(sapply(1:n, function(i) trapz(t, (Y[,i]-Y_est[,i])^2)))
  sumnorms = sum(sapply(1:n, function(i) trapz(t, Y[,i]^2)))
  apmse = mean((Y-Y_est)^2)
  
  return(list(MSE = sse/sumnorms, APMSE = apmse))
}

plot_cv = function(results){
  
  lambda_k = results$lambda_k
  k_range = results$k_range
  n_k = length(k_range)
  n = length(results$reord)
  mean = NULL
  se = NULL
  npred = NULL
  apmse_min = array(0, dim = c(n_k, 2))
  colnames(apmse_min) = c("log_lambda", "val")
  rownames(apmse_min) = k_range
  min = Inf
  idx_min = Inf
  idx = rep(k_range, each = 100)
  n_fold = max(results$fold)
  apsse = results$apsse
  apmse = array(sapply(1:n_fold, function(i) apsse[,,i]/sum(results$folds == i)),
                dim = c(100, n_k, n_fold))
  for(j in 1:n_k){
    mean_j = rowSums(apsse[,j,])/n
    mean = c(mean, mean_j)
    se = c(se, apply(apmse[,j,], 1, sd)/sqrt(n_fold))
    npred_j = results$res_data[[j]]$npred
    npred = c(npred, npred_j)
    idx_min_j = which(mean_j == min(mean_j))
    min_j = mean_j[idx_min_j]
    apmse_min[j,] = c(log(lambda_k[idx_min_j, j]), min_j)
    if(min_j < min){
      min = min_j
      npred_min = npred_j[idx_min_j]
    }
  }
  idx_min1se = which((mean < min+se[which(mean == min)]) & (npred <= npred_min))[1]
  k_min1se = idx[idx_min1se]
  k_min = k_range[apmse_min[,2] == min(apmse_min[,2])]
  min1se = list(log_lambda = log(as.vector(lambda_k)[idx_min1se]), val = mean[idx_min1se])
  min = list(log_lambda = apmse_min[apmse_min[,2] == min(apmse_min[,2]), 1], 
             val = apmse_min[apmse_min[,2] == min(apmse_min[,2]), 2])
  
  df = data.frame(mean = mean, lower = mean-se,
                  upper = mean+se, 
                  log_lambda = as.vector(log(lambda_k)),
                  label = idx)
  df2 = data.frame(log_lambda = apmse_min[,1], text_pos = max(mean),
                   label = k_range, 
                   text = "min")
  plot = ggplot(df, aes(x = log_lambda, y = mean, group = label)) + 
    xlab("log(lambda)") + ylab("APMSE") +
    facet_wrap(~label) +
    geom_errorbar(aes(ymin = lower, ymax = upper), col = "gray60") + 
    geom_point(col = "red") +
    geom_vline(data = df2, mapping = aes(xintercept = log_lambda), linetype = "dashed") +
    geom_text(data = df2, mapping = aes(x = log_lambda, y = text_pos, label = text),
              angle = 90, vjust = 1) +
    geom_vline(data.frame(x = min1se$log_lambda, label = k_min1se), 
               mapping = aes(xintercept = x), linetype = "dashed") +
    geom_text(data.frame(x = min1se$log_lambda, label = k_min1se, text = "min1se"), 
              mapping = aes(x = min1se$log_lambda, y = max(mean), label = text),
              angle = 90, vjust = 1) 
  
  return(list(apmse_min = apmse_min, k_min1se = k_min1se, k_min = k_min,
              min1se = min1se, min = min, plot = plot))
}

