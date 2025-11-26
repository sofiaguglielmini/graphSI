#' Data Splitting Inference for Graphical Models
#' @description Perform data-splitting inference for a selected graphical model.
#' @param X Data matrix (n x p)
#' @param j Index of the edge to perform inference on (in vech form)
#' @param nullvalue Null value for the hypothesis test
#' @param selected Output from graphSelect function
#' @param sandwich.variance Logical indicating whether to use sandwich variance estimator (default: FALSE)
#' @param alpha Significance level for confidence intervals (default: 0.05)
#' @return A list containing the p-value, lower and upper bounds of the confidence interval
graphInference_datasplitting <- function(X, j, nullvalue, selected,
                                sandwich.variance,
                                alpha){
  n2 <- nrow(X)
  E <- selected$E
  estimated <- graph_estimate(X = X, selected = selected, sandwich.variance = sandwich.variance)
  inference <- inference_Gaussian(estimated$theta_bar[E][j], sqrt(estimated$Sigma_E[j,j]/n2), nullvalue, alpha)
  inference
}

