#' Compute linear regression with constraints 
#' 
#' Compute the least square solution from linear regression with constraint 
#' \eqn{L*beta=c}.
#'
#' @param y tba
#' @param Z tba
#' @param L tba
#' @param c tba
#' @param intercept tba 
#' @return A list 
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