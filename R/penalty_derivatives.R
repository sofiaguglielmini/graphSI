#' Derivative of the Lasso Penalty
#'
#' This function computes the derivative of the Lasso penalty for a given input.
#' The Lasso penalty is defined as:
#' \deqn{P(x) = \lambda |x|}
#' where \eqn{\lambda} is the regularization parameter.
#' @param x Numeric value or vector for which to compute the derivative.
#' @param lambda Regularization parameter (non-negative scalar).
#' @return Numeric value or vector representing the derivative of the Lasso penalty at \code{x}.
#' @examples
#' lasso_derivative(2, lambda = 0.5)
#' lasso_derivative(c(-1, 0, 1), lambda = 1)
#' @export
lasso_derivative <- function(x, lambda){
  lambda*sign(x)
}

#' Derivative of the Elastic Net Penalty
#'
#' This function computes the derivative of the elastic net penalty for a given input.
#' The elastic net penalty is defined as:
#' \deqn{P(x) = \lambda \left( \gamma |x| + \frac{1 - \gamma}{2} x^2 \right)}
#' where \eqn{\lambda} is the regularization parameter and \eqn{\gamma} is the mixing parameter.
#'
#' @param x Numeric value or vector for which to compute the derivative.
#' @param lambda Regularization parameter (non-negative scalar).
#' @param gamma Mixing parameter between L1 and L2 penalties (between 0 and 1).
#' @return Numeric value or vector representing the derivative of the elastic net penalty at \code{x}.
#' @examples
#' elnet_derivative(2, lambda = 0.5, gamma = 0.3)
#' elnet_derivative(c(-1, 0, 1), lambda = 1, gamma = 0.5)
#' @export
elnet_derivative <- function(x, lambda, gamma){
  lambda*(gamma+(1-gamma)*abs(x))
}

#' Derivative of the SCAD Penalty
#'
#' This function computes the derivative of the Smoothly Clipped Absolute Deviation (SCAD) penalty for a given input.
#' The SCAD penalty is defined as:
#' \deqn{P(x) = \begin{cases}
#' \lambda |x| & \text{if } |x| \leq \lambda \\
#' \frac{2 \gamma \lambda |x| - x^2 - \lambda^2}{2(\gamma - 1)} & \text{if } \lambda < |x| \leq \gamma \lambda \\
#' \frac{(\gamma + 1) \lambda^2}{2} & \text{if } |x| > \gamma \lambda
#' \end{cases}}
#' where \eqn{\lambda} is the regularization parameter and \eqn{\gamma} is a parameter greater than 2.
#' @param x Numeric value or vector for which to compute the derivative.
#' @param lambda Regularization parameter (non-negative scalar).
#' @param gamma Parameter greater than 2 that controls the concavity of the penalty.
#' @return Numeric value or vector representing the derivative of the SCAD penalty at \code{x}.
#' @examples
#' scad_derivative(2, lambda = 0.5, gamma = 3)
#' scad_derivative(c(-1, 0, 1), lambda = 1, gamma = 4)
#' @export
scad_derivative <- function(x, lambda, gamma) {
  absx <- abs(x)
  ifelse(absx <= lambda, lambda,
                ifelse(absx > gamma * lambda, 0, (gamma * lambda - absx) / (gamma - 1)))
}

#' Derivative of the MCP Penalty
#'
#' This function computes the derivative of the Minimax Concave Penalty (MCP) for a given input.
#' The MCP penalty is defined as:
#' \deqn{P(x) = \begin{cases}
#' \lambda |x| - \frac{x^2}{2\gamma} & \text{if } |x| \leq \gamma \lambda \\
#' \frac{\gamma \lambda^2}{2} & \text{if } |x| > \gamma \lambda
#' \end{cases}}
#' where \eqn{\lambda} is the regularization parameter and \eqn{\gamma} is a parameter greater than 1.
#' @param x Numeric value or vector for which to compute the derivative.
#' @param lambda Regularization parameter (non-negative scalar).
#' @param gamma Parameter greater than 1 that controls the concavity of the penalty.
#' @return Numeric value or vector representing the derivative of the MCP penalty at \code{x}.
#' @examples
#' mcp_derivative(2, lambda = 0.5, gamma = 3)
#' mcp_derivative(c(-1, 0, 1), lambda = 1, gamma = 4)
#' @export
mcp_derivative <- function(x, lambda, gamma){
  absx <- abs(x)
  ifelse(absx <= gamma*lambda, lambda-absx/gamma, 0)
}
