loss_gradient <- function(Sigma, X){
  p <- nrow(Sigma)
  Sn <- cov(X)
  fastmatrix::dupl.prod(p, x=fastmatrix::vec(Sn - Sigma)/2, side="left", transposed=T)
}

loss_gradient_variance <- function(Sigma, X){
  p <- nrow(Sigma)
  n <- nrow(X)
  out <- 0
  for(h in 1:n){
    XXT <- crossprod(t(X[h,]))
    G <- fastmatrix::dupl.prod(p, x=fastmatrix::vec(Sigma - XXT)/2, side="left", transposed=T)
    out <- out + tcrossprod(G)/n
  }
  out
}

loss_hessian <- function(Sigma){
  p <- nrow(Sigma)
  fastmatrix::dupl.cross(p,x=fastmatrix::kronecker.prod(Sigma))/2
}
