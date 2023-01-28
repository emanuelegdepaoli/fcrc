#' Fitting linear regression models with linear constraints 
#' 
#' Compute the least square solution for estimating a linear regression model
#' \eqn{y=Z*beta} under the constraint \eqn{L*beta=c}.
#' 
#' The associated Lagrangian is \eqn{L(beta, u) = -t(y)*Z*beta + 0.5*t(beta)*t(Z)*Z*beta 
#' + t(u)*(L*beta-c)}.
#'
#' @param y A \eqn{n}-dimensional response vector.
#' @param Z A \eqn{n*p} matrix of predictors. 
#' @param L A \eqn{q*p} matrix.
#' @param c A \eqn{q}-dimensional vector.
#' @param intercept logical. Should an intercept be included in the model? 
#'                  default is TRUE.
#' @return A list containing \eqn{beta}, the intercept \eqn{beta_0} and \eqn{u}.
#' @export
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
  
  return(list(beta = beta, beta0 = beta0, u = u))
} 