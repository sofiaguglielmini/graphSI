#' Get Gaussian p-value and confidence interval
#' @param theta Estimated parameter
#' @param se Standard error of the estimated parameter
#' @param nullvalue Null hypothesis value for the parameter
#' @param alpha Significance level for confidence intervals
#' @return A list containing the p-value, lower and upper bounds of the confidence interval
inference_Gaussian <- function(theta, se, nullvalue, alpha){
  z_score <- (theta - nullvalue)/se
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  z_alpha <- qnorm(1 - alpha/2)
  ci_lower <- theta - z_alpha * se
  ci_upper <- theta + z_alpha * se
  list(p_value = p_value,
       ci_lower = ci_lower,
       ci_upper = ci_upper)
}

#' Get truncated Gaussian p-value and confidence interval
#' @param theta Estimated parameter
#' @param j Index of the parameter in the full parameter vector
#' @param Var Covariance matrix of the estimated parameters
#' @param nullvalue Null hypothesis value for the parameter
#' @param A Affine constraint matrix
#' @param b Affine constraint vector
#' @param alpha Significance level for confidence intervals
#' @return A list containing the p-value, lower and upper bounds of the confidence interval
inference_truncatedGaussian <- function(theta, j, Var, nullvalue, A, b, alpha){
  truncgauss <- selectiveInference::TG.pvalue(Z = theta, A = A, b =  b,
                                             eta = c(rep(0, j-1), 1, rep(0, length(theta)-j)),
                                             null_value = nullvalue,
                                             Sigma = Var)
  p_value <- 2 * min(truncgauss$pv, 1 - truncgauss$pv)
  ci_lower <- truncgauss$int[1]
  ci_upper <- truncgauss$int[2]
  list(p_value = p_value,
       ci_lower = ci_lower,
       ci_upper = ci_upper)
}

