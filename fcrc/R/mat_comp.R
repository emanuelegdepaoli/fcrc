#' Matrix computation for model estimation
#' 
#' Compute the matrices required for estimating the functional model, using B-splines
#' of order 4 as basis functions. 
#'
#' @param Z A \eqn{n_t*p*n} array, containing \eqn{n} samples of \eqn{p} functional
#'          predictors observed at \eqn{n_t} points. 
#' @param Y A \eqn{n_t*n} array, containing \eqn{n} samples of the functional 
#'          response at \eqn{n_t} points. 
#' @param t A vector containing the values at which functions are observed.
#' @param k The number of basis function.
#' @param Zc optional. An array of \eqn{p_c} control variables, with same structure of Z. The
#'           first predictor must be a vector of ones. The default value is an array of ones,
#'           to estimate the intercept. 
#' @return A list containing the matrices \eqn{J}, \eqn{K}, \eqn{M}, \eqn{P} and 
#'         other quantities used for model estimation. 
#' @importFrom fda create.bspline.basis eval.basis 
#' @export
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
  
  # # computation R int phi(t)*phi(t)^t dt 
  # for (i in 1:k) {
  #   for (j in 1:k) {
  #     R[i,j] = trapz(t, Phi[,i]*Phi[,j]) 
  #   }
  # }
  # 
  # # int d2 phi(t) * d2 phi(t)
  # S = inprod(beta_basis, beta_basis, Lfdobj1 = 2, Lfdobj2 = 2)
  
  return(list(J = J, K = K, M = M, P = P, Q = Q, p = p, pc = pc, n = n, k = k, 
              Phi = Phi, beta_basis = beta_basis))
}