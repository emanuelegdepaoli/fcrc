#' Constraint matrix computation 
#' 
#' Compute the matrix representing the set of linear constraints. 
#'
#' @param p_vec A vector of length \eqn{q} whose elements are the number of 
#'              components of each composition. 
#' @return The matrix \eqn{L}.
#' @export
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