#' Error computation
#' 
#' Compute the mean square error (MSE) and the average prediction MSE (APMSE). 
#'
#' @param beta to be added
#' @param betac to be added
#' @param Y to be added
#' @param Z to be added
#' @param Phi to de added 
#' @param t to be added
#' @param Zc to be added
#' @return A list containing the MSE and the APMSE. 
#' @export
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