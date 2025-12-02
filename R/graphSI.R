#' @importFrom stats cov pnorm qnorm
#'
#' @title Regularized undirected graph selection
#' @description Select an undirected graphical model using regularized estimation methods.
#' @param data Data matrix (n x p)
#' @param lambda Tuning parameter for the penalty (default: sqrt(2*log(p)/n)). For valid selective inference, lambda should be chosen independently of the data used for selection.
#' @param gamma Additional tuning parameter for certain penalties (default: 0.5 for elastic net, 3.7 for SCAD, 2.0 for MCP)
#' @param data.splitting Logical indicating whether to use data splitting for inference (default: FALSE)
#' @param split.proportion Proportion of data to use for selection if data.splitting is TRUE (default: 0.5)
#' @param loss Loss function to use (default: "Gaussian")
#' @param penalty Penalty to use for graph selection (default: "lasso"). Options: "lasso", "elastic net", "scad", "mcp".
#' @param penalize.diagonal Logical indicating whether to penalize diagonal elements of the precision matrix (default: FALSE)
#' @param seed Random seed for data splitting (default: NULL)
#' @return An S3 object of class 'graphSelect' with user-facing elements: adjacency matrix, selected indices, data splitting info. Internal elements also stored but hidden from print().
#' @export
graphSelect <- function(data, lambda = NULL, gamma = NULL,
                        data.splitting = FALSE,
                        split.proportion = NULL,
                        loss = c("Gaussian"),
                        penalty = c("lasso", "elastic net", "scad", "mcp"),
                        penalize.diagonal = FALSE,
                        seed = NULL){

  penalty <- match.arg(penalty)
  loss <- match.arg(loss)

  if(data.splitting){
    if(!is.null(seed)) set.seed(seed)
    if(is.null(split.proportion)) split.proportion <- 0.5
    n <- nrow(data)
    n1 <- floor(n * split.proportion)
    X <- data[1:n1, ]
  } else{
    X <- data
    n1 <- NULL
  }

  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X)
  Sn <- cov(X)

  if(is.null(lambda)) lambda <- sqrt(2*log(p)/n)
  if(is.null(gamma)){
    gamma <- switch(penalty,
                    "elastic net" = 0.5,
                    "scad" = 3.7,
                    "mcp" = 2.0,
                    NA)
  }

  if(loss == "Gaussian"){
    if(penalty == "lasso"){
      Lambda <- matrix(lambda, nrow=p, ncol=p)
      if(!penalize.diagonal) diag(Lambda) <- 0
      selection_step <- glassoFast::glassoFast(Sn, Lambda)
      Theta_hat <- selection_step$wi
      Sigma_hat <- selection_step$w
    } else if(penalty == "elastic net"){
      selection_step <- GLassoElnetFast::gelnet(Sn, lambda, gamma, penalize.diagonal=penalize.diagonal, Target=0*diag(p))
      Theta_hat <- selection_step$Theta
      Sigma_hat <- selection_step$W
    } else if(penalty %in% c("scad", "mcp")){
      selection_step <- suppressMessages(GGMncv::ggmncv(Sn, n, penalty=penalty, lambda=lambda, gamma=gamma,
                                       penalize_diagonal=penalize.diagonal, LLA = TRUE, maxit=1e05, thr=1e-06, initial=glassoFast::glassoFast(Sn, lambda)$wi))
      Theta_hat <- selection_step$Theta
      Sigma_hat <- selection_step$Sigma
    }
  } else stop("Loss function not recognized.")

  theta_hat <- fastmatrix::vech(Theta_hat)
  E <- which(fastmatrix::vech(1 * (Theta_hat != 0)) != 0)
  nE <- which(fastmatrix::vech(1 * (Theta_hat != 0)) == 0)

  selected.indices <- which(Theta_hat != 0 & row(Theta_hat) >= col(Theta_hat), arr.ind = TRUE)
  selected.edges <- which(Theta_hat != 0 & row(Theta_hat) > col(Theta_hat), arr.ind = TRUE)

  if (!penalize.diagonal) {
    selected.indices <- selected.edges
  }

  # create S3 object with both user-facing and internal elements
  out <- structure(
    list(

      # user-facing elements
      adjacency.matrix = Theta_hat != 0,
      selected.edges = selected.edges,
      data.splitting = data.splitting,
      split.proportion = split.proportion,

      # internal elements
      selected.indices = selected.indices,
      Theta_hat = Theta_hat,
      Sigma_hat = Sigma_hat,
      theta_hat = theta_hat,
      E = E,
      nE = nE,
      penalty = penalty,
      loss = loss,
      lambda = lambda,
      gamma = gamma,
      penalize.diagonal = penalize.diagonal,
      n1 = n1
    ),
    class = "graphSelect"
  )
  return(out)
}

# custom print method for graphSelect objects
# only shows user-facing elements
#' @export
print.graphSelect <- function(x, ...) {
  cat("Selected Graph:\n")
  print(x$adjacency.matrix)
  cat("\nSelected edges:\n")
  print(x$selected.edges)
  if(x$data.splitting) {
    cat("\nData splitting used, proportion:", x$split.proportion, "\n")
  }
}

#' @title Inference for edges in the selected graph
#' @description Perform inference for a selected edge in the graphical model using either polyhedral or data-splitting methods.
#' @param data Data matrix (n x p)
#' @param selected Output from graphSelect function
#' @param j Index of the edge to perform inference on (in vech form)
#' @param nullvalue Null value for the hypothesis test
#' @param sandwich.variance Logical indicating whether to use sandwich variance estimator (default: FALSE)
#' @param alpha Significance level for confidence intervals (default: 0.05)
#' @param seed Random seed for data splitting (default: NULL)
#' @return A list containing the p-value, lower and upper bounds of the confidence interval
#' @export
graphInference <- function(data, selected, j, nullvalue,
                    sandwich.variance = FALSE,
                    alpha = 0.05, seed = NULL){
  if(selected$data.splitting){
    X <- data[(selected$n1+1):nrow(data), ]
    method <- "Data splitting"
    message("Using data splitting for inference.")
  } else{
    X <- data
    method <- "Polyhedral"
    message("Using polyhedral method for inference.")
  }
  X <- scale(X)
  if(method == "Polyhedral"){
    out <- graphInference_polyhedral(X = X, j, nullvalue, selected = selected,
                              sandwich.variance = sandwich.variance,
                              alpha = alpha)
  } else if(method == "Data splitting"){
    out <- graphInference_datasplitting(X = X, j, nullvalue, selected = selected,
                                 sandwich.variance = sandwich.variance,
                                 alpha = alpha)
  }
  out
}

