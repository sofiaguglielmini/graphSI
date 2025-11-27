#' @title Polyhedral inference for graphical models
#' @description Perform hypothesis testing and construct confidence intervals for edges in a selected graphical model using the polyhedral method.
#' @param X Data matrix (n x p)
#' @param j Index of the edge to perform inference on (in vech form)
#' @param nullvalue Null value for the hypothesis test
#' @param selected Output from graphSelect function
#' @param sandwich.variance Logical indicating whether to use sandwich variance estimator (default: FALSE)
#' @param alpha Significance level for confidence intervals (default: 0.05)
#' @return A list containing the p-value, lower and upper bounds of the confidence interval
graphInference_polyhedral <- function(X, j, nullvalue, selected,
                             sandwich.variance,
                             alpha){
  E <- selected$E
  estimated <- graph_estimate(X = X,
                              selected = selected,
                              sandwich.variance = sandwich.variance)

  conditional <- graph_polyhedral_conditioning(X = X,
                                               selected = selected,
                                               estimated = estimated)

  inference <- inference_truncatedGaussian(conditional$theta_onestep[E], j, conditional$Sigma_E/nrow(X),
                                           nullvalue, conditional$A, conditional$b, alpha)
  inference
}

#' @title Polyhedral conditioning for graphical models
#' @description Compute the affine constraints and one-step estimator for polyhedral inference in graphical models.
#' @param X Data matrix (n x p)
#' @param selected Output from graphSelect function
#' @param estimated Output from graph_estimate function
#' @return A list containing the affine constraint matrix A, vector b, one-step estimator theta_onestep, and covariance matrix Sigma_E
graph_polyhedral_conditioning <- function(X, selected, estimated){
  n <- nrow(X)
  p <- ncol(X)

  Hn <- estimated$Hn
  Jn <- estimated$Jn
  Gn <- estimated$Gn
  Sigma_E <- estimated$Sigma_E
  Sn <- estimated$Sn
  sandwich.variance <- estimated$sandwich.variance

  E <- selected$E
  nE <- selected$nE
  diagsE <- which(selected$selected.indices[,1]==selected$selected.indices[,2])
  theta_hat <- selected$theta_hat
  Sigma_hat <- selected$Sigma_hat
  lambda <- selected$lambda
  gamma <- selected$gamma
  penalty <- selected$penalty
  loss <- selected$loss
  penalize.diagonal <- selected$penalize.diagonal

  sigma_hat <- fastmatrix::vech(Sigma_hat)
  sn <- fastmatrix::vech(Sn)
  H_hat_EE_inv <- solve(loss_hessian(Sigma_hat)[E,E])
  Jn_EE_inv <- solve(Jn[E,E, drop=F])
  G_hat <- loss_gradient(Sigma_hat, X)
  signE <- sign(theta_hat[E])

  # # Compute subgradient
  # subgrad <- rep(NA, length(theta_hat))
  # # subgrad[E] <- sign(theta_hat[E])
  # if(loss == "Gaussian"){
  #   for(edge in 1:length(theta_hat)){
  #     subgrad[edge] <- (sigma_hat[edge] - sn[edge])/penalty_derivative(theta_hat[edge], penalty=penalty, lambda, gamma)
  #   }
  # }
  # if(!penalize.diagonal){
  #   diag_indices <- which(fastmatrix::vech_diag_indices(p))
  #   subgrad[diag_indices] <- 0
  # }
  # print(subgrad[E])
  # print(subgrad[nE])

  # One-step estimator
  theta_onestep <- rep(0, length(theta_hat))
  theta_onestep[E] <- theta_hat[E] - H_hat_EE_inv %*% G_hat[E]

  # check lasso kkt
  # diagonal elements of E

  # print(cbind(G_hat[E], -lambda * signE, selected$selected.indices))

  # Nuisance parameter
  A_perp <- Hn[nE,E, drop=F] - Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E, drop=F]
  theta_perp <- rep(0, length(theta_hat))
  # theta_perp[nE] <- Gn[nE] + A_perp %*% theta_onestep[E]
  theta_perp[nE] <- G_hat[nE] + lambda * Hn[nE,E, drop=F]%*% H_hat_EE_inv %*% signE + A_perp %*% theta_onestep[E]

  V_hat_E <- abs(theta_hat[E])

  # Affine constraints
  if(penalty == "lasso"){
    A1 <- - diag(signE, nrow = length(E))
    b1 <- - lambda * diag(signE, nrow = length(E)) %*% H_hat_EE_inv %*% signE

    b_term <- lambda * Hn[nE,E, drop=F] %*% H_hat_EE_inv %*% signE
    A00 <- - A_perp
    b00 <- lambda - theta_perp[nE] - b_term
    A01 <- A_perp
    b01 <- lambda + theta_perp[nE] + b_term

  } else if(penalty == "elastic net"){
    K_inv <- solve(diag(length(E)) + lambda * (1-gamma) * H_hat_EE_inv)

    A1 <- - diag(signE, nrow = length(E)) %*% K_inv
    b1 <- - lambda * gamma * diag(signE, nrow = length(E)) %*% K_inv %*% H_hat_EE_inv %*% signE

    A_term1 <- Hn[nE,E, drop=F] %*% K_inv
    A_term2 <- Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E, drop=F]
    b_term <- lambda * gamma * Hn[nE,E, drop=F] %*% K_inv %*% H_hat_EE_inv %*% signE
    A00 <- - A_term1 + A_term2
    b00 <- lambda * gamma - theta_perp[nE] - b_term
    A01 <- - A00
    b01 <- lambda * gamma + theta_perp[nE] + b_term

  } else if(penalty=="scad"){

    I1 <- 1 * diag(V_hat_E > lambda & V_hat_E <= gamma * lambda)
    if(length(I1) == 0) I1 <- 0

    I2 <- 1 * diag(V_hat_E <= lambda)
    if(length(I2) == 0) I2 <- 0

    K_inv <- round(solve(diag(length(E)) - 1/(gamma-1) * H_hat_EE_inv %*% I1), 5)

    A1 <- - diag(signE, nrow = length(E)) %*% K_inv
    b1 <- - lambda * diag(signE, nrow = length(E)) %*% K_inv %*% H_hat_EE_inv %*% (I2 + gamma/(gamma-1) * I1) %*% signE

    b_term <- lambda * Hn[nE,E, drop=F] %*% K_inv %*% H_hat_EE_inv %*% (I2 + gamma/(gamma-1) * I1) %*% signE
    A00 <- - Hn[nE,E, drop=F] %*% K_inv + Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E, drop=F]
    b00 <- lambda - theta_perp[nE] - b_term
    A01 <- - A00
    b01 <- lambda + theta_perp[nE] + b_term

    E0 <- which(V_hat_E > lambda*gamma)
    E2 <- which(V_hat_E <= lambda)
    E1 <- which(V_hat_E > lambda & V_hat_E <= lambda * gamma)

    M1 <- diag(length(E1)) - H_hat_EE_inv[E1,E1, drop=F] / (gamma - 1)

    if(length(E1)!=0){
      M1_inv <- round(solve(M1), 5)
    } else{ M1_inv <- matrix(NA, length(E1), length(E1)) }

    K0 <- cbind(diag(length(E0)), H_hat_EE_inv[E0,E1, drop=F] %*% M1_inv / (gamma - 1), matrix(0, length(E0), length(E2)))
    k0 <- (lambda * gamma / (gamma - 1) * H_hat_EE_inv[E0,E1, drop=F] %*% signE[E1]
           + lambda * H_hat_EE_inv[E0,E2, drop=F] %*% signE[E2]
           + lambda * gamma / (gamma - 1)^2 * H_hat_EE_inv[E0,E1, drop=F] %*% M1_inv %*% H_hat_EE_inv[E1,E1, drop=F] %*% signE[E1]
           + lambda / (gamma - 1) * H_hat_EE_inv[E0,E1, drop=F] %*% M1_inv %*% H_hat_EE_inv[E1,E2, drop=F] %*% signE[E2])

    if(length(E1) != 0){
      K1 <- cbind(matrix(0, length(E1), length(E0)), M1_inv, matrix(0, length(E1), length(E2)))
      k1 <- (lambda * gamma / (gamma - 1) * M1_inv %*% H_hat_EE_inv[E1,E1, drop=F] %*% signE[E1]
             + lambda * M1_inv %*% H_hat_EE_inv[E1,E2, drop=F] %*% signE[E2])
    } else{
      K1 <- 0
      k1 <- 0
    }
    K2 <- cbind(matrix(0, length(E2), length(E0)), H_hat_EE_inv[E2,E1, drop=F] %*% M1_inv / (gamma - 1), diag(length(E2)))
    k2 <- (lambda * gamma / (gamma - 1) * H_hat_EE_inv[E2,E1, drop=F] %*% signE[E1]
           + lambda * H_hat_EE_inv[E2,E2, drop=F] %*% signE[E2]
           + lambda * gamma / (gamma - 1)^2 * H_hat_EE_inv[E2,E1, drop=F] %*% M1_inv %*% H_hat_EE_inv[E1,E1, drop=F] %*% signE[E1]
           + lambda / (gamma - 1) * H_hat_EE_inv[E2,E1, drop=F] %*% M1_inv %*% H_hat_EE_inv[E1,E2, drop=F] %*% signE[E2])

    A00_0 <- - diag(signE[E0], nrow = length(E0)) %*% K0
    b00_0 <- - gamma * lambda - diag(signE[E0], nrow = length(E0)) %*% k0

    A00_10 <- diag(signE[E1], nrow = length(E1)) %*% K1
    b00_10 <- gamma * lambda + diag(signE[E1], nrow = length(E1)) %*% k1

    A00_11 <- - diag(signE[E1], nrow = length(E1)) %*% K1
    b00_11 <- - lambda - diag(signE[E1], nrow = length(E1)) %*% k1
    if(length(E1)==0){
      A00_11 <- NULL
      b00_11 <- NULL
      A00_10 <- NULL
      b00_10 <- NULL
    }
    A00_2 <- diag(signE[E2], nrow = length(E2)) %*% K2
    b00_2 <- lambda + diag(signE[E2], nrow = length(E2)) %*% k2

    og_order <- order(c(E0, E1, E2))

    A00_cond <- rbind(rbind(A00_0, A00_10, A00_2)[og_order,og_order], rbind(A00_0, A00_11, A00_2)[og_order,og_order])
    b00_cond <- c(c(b00_0, b00_10, b00_2)[og_order], c(b00_0, b00_11, b00_2)[og_order])

    A00 <- rbind(A00, A00_cond)
    b00 <- c(b00, b00_cond)

  } else if(penalty=="mcp"){

    I1 <- 1 * diag(V_hat_E<=lambda*gamma)
    if(length(I1)==0) I1 <- 0

    K_inv <- solve(diag(length(E)) - (1 / gamma) * H_hat_EE_inv %*% I1)

    A1 <- - diag(signE, nrow = length(E)) %*% K_inv
    b1 <- -lambda * diag(signE, nrow = length(E)) %*% K_inv %*% H_hat_EE_inv %*% I1 %*% signE

    A_term1 <- Hn[nE,E, drop=F] %*% K_inv
    A_term2 <- Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E]
    b_term <- lambda * Hn[nE,E, drop=F] %*% K_inv %*% H_hat_EE_inv %*% I1 %*% signE

    A00 <- - A_term1 + A_term2
    b00 <- lambda - theta_perp[nE] - b_term

    A01 <- -A00
    b01 <- lambda + theta_perp[nE] + b_term

    E1 <- which(V_hat_E <= lambda * gamma)
    E0 <- which(V_hat_E > lambda * gamma)

    M1 <- diag(length(E1)) - H_hat_EE_inv[E1,E1] / gamma
    M1_inv <- solve(M1)

    K1 <- cbind(matrix(0, length(E1), length(E0)), M1_inv)
    k1 <- lambda * M1_inv %*% H_hat_EE_inv[E1,E1, drop=F] %*% signE[E1]

    K0 <- cbind(diag(length(E0)), H_hat_EE_inv[E0,E1, drop=F] %*% M1_inv / gamma)
    k0 <- (lambda * H_hat_EE_inv[E0,E1, drop=F] %*% signE[E1]
           + (lambda/gamma) * H_hat_EE_inv[E0,E1, drop=F] %*% M1_inv %*% H_hat_EE_inv[E1,E1, drop=F] %*% signE[E1])

    A00_0 <- - diag(signE[E0], nrow=length(E0)) %*% K0
    b00_0 <- - gamma * lambda - diag(signE[E0], nrow=length(E0)) %*% k0

    A00_1 <- diag(signE[E1], nrow=length(E1)) %*% K1
    b00_1 <- gamma * lambda + diag(signE[E1], nrow=length(E1)) %*% k1

    og_order <- order(c(E0, E1))

    A00_cond <- rbind(rbind(A00_0, A00_1)[og_order,og_order])
    b00_cond <- c(c(b00_0, b00_1)[og_order])

    A00 <- rbind(A00, A00_cond)
    b00 <- c(b00, b00_cond)
  }
  A <- rbind(A1, A00, A01)
  b <- c(b1, b00, b01)
  print(cbind(A %*% theta_onestep[E], b, A %*% theta_onestep[E] - b))
  if(!all(A %*% theta_onestep[E] <= b)) print("Affine constraint not satisfied.")

  list(A = A, b = b, theta_onestep = theta_onestep, Sigma_E = Sigma_E)
}
