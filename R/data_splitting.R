#' Data Splitting Inference for Graphical Models
#' @description Perform data-splitting inference for a selected graphical model.
#' @param X Data matrix (n x p)
#' @param j Indices of the edges to perform inference on (in vech form)
#' @param nullvalue Null value for the hypothesis test
#' @param selected Output from graphSelect function
#' @param sandwich.variance Logical indicating whether to use sandwich variance estimator (default: FALSE)
#' @param alpha Significance level for confidence intervals (default: 0.05)
#' @return A list containing the p-value, lower and upper bounds of the confidence interval
graphInference_datasplitting <- function(X, j, nullvalue, selected,
                                         sandwich.variance = FALSE,
                                         alpha = 0.05){
  n2 <- nrow(X)
  E <- selected$E
  estimated <- graph_estimate(X = X, selected = selected, sandwich.variance = sandwich.variance)
  if(is.character(j) && j == "all"){
    p <- ncol(X)
    diags <- cumsum(p:1) - (p - 1:p)
    j <- 1:length(E)
    j <- j[!E %in% diags]
  }
  inference <- lapply(j, function(idx) {
    inference_Gaussian(estimated$theta_bar[E][idx], sqrt(estimated$Sigma_E[idx, idx]/n2), nullvalue, alpha)
  })
  inference <- do.call(rbind, inference)
  inference <- cbind(selected$selected.edges[E,][j,], inference)
  return(list(inference = inference, estimated.graph = estimated$Theta_bar))
}

