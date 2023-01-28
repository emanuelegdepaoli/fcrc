#' Error computation
#' 
#' Compute the mean square error (MSE) and the average prediction MSE (APMSE). 
#'
#' @param beta The estimated \eqn{beta} (for a single \eqn{lambda})
#' @param betac The estimated \eqn{betac} (for a single \eqn{lambda}).
#' @param Z A \eqn{n_t*p*n} array, containing \eqn{n} samples of \eqn{p} functional
#'          predictors observed at \eqn{n_t} points. 
#' @param Y A \eqn{n_t*n} array, containing \eqn{n} samples of the functional 
#'          response at \eqn{n_t} points. 
#' @param Phi The matrix \eqn{Phi} from \eqn{mat_comp}. 
#' @param t A vector containing the values at which functions are observed.
#' @param Zc optional. An array of \eqn{p_c} control variables, with same structure of Z. The
#'           first predictor must be a vector of ones. The default value is an array of ones,
#'           to estimate the intercept.
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