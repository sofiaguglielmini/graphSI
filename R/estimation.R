#' Graph estimation within the selected graph
#' @description Estimate the precision matrix and variance within the selected graph structure.
#' @param X Data matrix (n x p)
#' @param selected Output from graphSelect function
#' @param sandwich.variance Logical indicating whether to use sandwich variance estimator (default: FALSE)
#' @return A list containing the sample covariance matrix Sn, one-step
graph_estimate <- function(X, selected, sandwich.variance){
  n <- nrow(X)
  p <- ncol(X)
  Sn <- cov(X)

  E <- selected$E
  loss <- selected$loss

  if(loss == "Gaussian"){

    # Fit within the selected graph
    refitting_step <- GGMncv::constrained(Sn, selected$adjacency.matrix)
    Theta_bar <- refitting_step$Theta
    Sigma_bar <- refitting_step$Sigma
    theta_bar <- fastmatrix::vech(Theta_bar)

    # Variance estimator
    Hn <- loss_hessian(Sigma_bar)
    Gn <- loss_gradient(Sigma_bar, X)
    Hn_EE_inv <- solve(Hn[E,E])
    if(sandwich.variance){
      Jn <- loss_gradient_variance(Sigma_bar, X)
      Sigma_E <- Hn_EE_inv %*% Jn[E,E] %*% Hn_EE_inv
    } else {
      Jn <- Hn
      Sigma_E <- Hn_EE_inv
    }
  }
  list(Sn = Sn,
       theta_bar = theta_bar,
       Sigma_E = Sigma_E,
       Jn = Jn,
       Hn = Hn,
       Gn = Gn,
       sandwich.variance = sandwich.variance)
}
