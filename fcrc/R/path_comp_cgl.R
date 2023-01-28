#' Compute regularization path 
#' 
#' Compute the regularization path to be used by \link{alm_cgl_path}.
#'
#' @param matrices A list from \link{mat_comp}.
#' @return A vector with the regularization path. 
#' @export
# Regularization path
path_comp_cgl = function(matrices){
  p = matrices$p
  k = matrices$k
  ncoef = p*k
  idx = rep(1:p, each = k)
  grps = unique(idx)
  n = matrices$n
  # scale K and J
  qtminvp = t(matrices$Q)%*%solve(matrices$M)%*%matrices$P
  xtystar = matrices$J - qtminvp
  lambda_max = max(sapply(grps, function(j) norm(xtystar[idx == j], "2")))
  # grid as in Regularization Paths for Generalized Linear Models
  # Friedman et al, 2010
  lambda_path = exp(seq(log(lambda_max), log(lambda_max*0.001),
                        length.out = 100))
  
  return(lambda_path)
}