#' Compute the gradient of the Gaussian loss function
#' with respect to the inverse covariance matrix
#' @param Sigma A p x p covariance matrix
#' @param X A n x p data matrix
#' @return A p x p matrix representing the gradient of the Gaussian loss function
#' @export
gaussian_loss_gradient <- function(Sigma, X){
  p <- nrow(Sigma)
  n <- nrow(X)
  out <- 0
  for(h in 1:n){
    XXT <- crossprod(t(X[h,]))
    G <- fastmatrix::dupl.prod(p, x=fastmatrix::vec(-Sigma+XXT)/2, side="left", transposed=T)
    out <- out + tcrossprod(G)/n
  }
  out
}

#' Compute the Hessian of the Gaussian loss function
#' with respect to the inverse covariance matrix
#' @param Sigma A p x p covariance matrix
#' @return A (p^2) x (p^2) matrix representing the Hessian of the Gaussian loss function
#' @export
gaussian_loss_hessian <- function(Sigma){
  p <- nrow(Sigma)
  fastmatrix::dupl.cross(p,x=fastmatrix::kronecker.prod(Sigma))/2
}
