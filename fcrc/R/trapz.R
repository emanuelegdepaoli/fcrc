#' Trapezoidal Integration
#' 
#' Compute the area of a function with values \eqn{y} at points \eqn{x}, using the
#' trapezoidal rule.
#'
#' @param x \eqn{x}-coordinates of points.
#' @param y \eqn{y}-coordinates of points.
#' @return Approximated integral of the function. 

trapz = function(x, y){
  l = length(x)
  xs = x[-1]-x[-l]
  ys = 0.5*(y[-l]+y[-1])
  return(sum(xs*ys))
}