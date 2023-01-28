#' Compute functional concurrent regression with constraints 
#' 
#' Compute the least square solution from functional concurrent regression with constraint 
#' \eqn{L*beta=c}.
#' 
#' The Lagrangian is \eqn{L(beta, u) = -t(y)*x*beta+t(beta)*t(x)*x*beta+t(u)*(Ltilda*beta-c)}.
#'
#' @param matrices tba.
#' @param Ltilda tba.
#' @param c tba.
#' @return A list.
#' @export
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