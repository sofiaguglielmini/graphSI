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

  inference <- inference_truncatedGaussian(conditional$theta_onestepE, j, conditional$Sigma_E/nrow(X),
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
  diags <- cumsum(p:1) - (p - 1:p)
  offdiags <- setdiff(1:(p*(p+1)/2), diags)

  E <- selected$E
  nE <- selected$nE
  theta_hat <- selected$theta_hat
  Sigma_hat <- selected$Sigma_hat
  lambda <- selected$lambda
  gamma <- selected$gamma
  penalty <- selected$penalty
  loss <- selected$loss
  penalize.diagonal <- selected$penalize.diagonal

  Sigma_E <- estimated$Sigma_E
  Sn <- estimated$Sn
  sandwich.variance <- estimated$sandwich.variance

  # If diagonal is not penalized, we don't need to adjust it
  Ep <- E
  if(!penalize.diagonal) Ep <- setdiff(E, diags)

  H_hat_EE_inv <- solve(loss_hessian(Sigma_hat)[E,E])
  H_hat_EpEp_inv <- H_hat_EE_inv
  if(!penalize.diagonal) H_hat_EpEp_inv <- H_hat_EE_inv[which(E %in% offdiags), which(E %in% offdiags), drop=F]
  G_hat <- loss_gradient(Sigma_hat, X)
  Sigma_Ep <- Sigma_E
  if(!penalize.diagonal) Sigma_Ep <- Sigma_E[which(E %in% offdiags), which(E %in% offdiags), drop=F]

  p_prime <- penalty_derivative(theta_hat, penalty=penalty, lambda, gamma)
  if(!penalize.diagonal) p_prime[diags] <- 0

  # Compute subgradient
  subgrad <- rep(0, length(theta_hat))
  if(loss == "Gaussian"){
    for(edge in 1:length(theta_hat)){

      # if the penalty is 0, we don't need to compute the subgradient
      if(p_prime[edge] == 0) next

      subgrad[edge] <- fastmatrix::vech(Sigma_hat - Sn)[edge]/p_prime[edge]
    }
  }
  # Subgrad <- rockchalk::vech2mat(subgrad)
  # subgrad <- fastmatrix::dupl.prod(p, x=fastmatrix::vec(Subgrad)/2, side="left", transposed=T)

  # One-step estimator
  theta_onestep <- theta_hat
  theta_onestep[Ep] <- theta_hat[Ep] - H_hat_EpEp_inv %*% G_hat[Ep]

  Hn <- loss_hessian(Sigma_hat)
  Jn <- if(sandwich.variance) loss_gradient_variance(Sigma_hat, X) else Hn

  # Nuisance parameter
  Jn_EE_inv <- solve(Jn[E,E, drop=F])
  Jn_EpEp_inv <- Jn_EE_inv
  if(!penalize.diagonal) Jn_EpEp_inv <- Jn_EE_inv[which(E %in% offdiags), which(E %in% offdiags), drop=F]
  A_perp <- Hn[nE,Ep, drop=F] - Jn[nE,Ep, drop=F] %*% Jn_EpEp_inv %*% Hn[Ep,Ep, drop=F]
  theta_perp <- rep(0, length(theta_hat))

  V_hat_E <- abs(theta_hat[Ep])

  # Affine constraints
  if(penalty == "lasso"){
    theta_perp[nE] <- lambda * subgrad[nE] - lambda*Hn[nE,Ep, drop=F]%*%H_hat_EpEp_inv%*%subgrad[Ep] + A_perp %*% theta_onestep[Ep]

    A1 <- - diag(sign(theta_hat[Ep]), nrow = length(Ep))
    b1 <- - lambda * diag(sign(theta_hat[Ep]), nrow = length(Ep)) %*% H_hat_EpEp_inv %*% subgrad[Ep]

    b_term <- lambda * Hn[nE,Ep, drop=F] %*% H_hat_EpEp_inv %*% subgrad[Ep]
    A00 <- - A_perp
    b00 <- lambda - theta_perp[nE] - b_term
    A01 <- A_perp
    b01 <- lambda + theta_perp[nE] + b_term

  } else if(penalty == "elastic net"){
    K_E <- diag(length(E)) + lambda * (1 - gamma) *  H_hat_EE_inv
    K_Ep <- K_E[which(E %in% Ep), which(E %in% Ep), drop=F]
    K_inv <- solve(K_Ep)

    theta_perp[nE] <- lambda * gamma * subgrad[nE] + Hn[nE,Ep, drop=F] %*% (K_inv %*% theta_onestep[Ep] - lambda * gamma * K_inv %*% H_hat_EpEp_inv %*% subgrad[Ep]) - Jn[nE,Ep, drop=F] %*% Jn_EpEp_inv %*% Hn[Ep,Ep] %*% theta_onestep[Ep]

    A1 <- - diag(sign(theta_hat[Ep]), nrow = length(Ep)) %*% K_inv
    b1 <- - lambda * gamma * diag(sign(theta_hat[Ep]), nrow = length(Ep)) %*% K_inv %*% H_hat_EpEp_inv %*% subgrad[Ep]

    b_term <- lambda * gamma * Hn[nE,Ep, drop=F] %*% K_inv %*% H_hat_EpEp_inv %*% subgrad[Ep]
    A00 <- - Hn[nE,Ep, drop=F] %*% K_inv + Jn[nE,Ep, drop=F] %*% Jn_EpEp_inv %*% Hn[Ep,Ep, drop=F]
    b00 <- lambda * gamma - theta_perp[nE] - b_term
    A01 <- - A00
    b01 <- lambda * gamma + theta_perp[nE] + b_term

  } else if(penalty=="scad"){

    theta_perp[nE] <- lambda * subgrad[nE] + Hn[nE,Ep, drop=F] %*% theta_hat[Ep] - Jn[nE,Ep, drop=F] %*% Jn_EpEp_inv %*% Hn[Ep,Ep, drop=F] %*% theta_onestep[Ep]

    I1 <- 1 * diag(V_hat_E > lambda & V_hat_E <= gamma * lambda)
    if(length(I1) == 0) I1 <- 0

    I2 <- 1 * diag(V_hat_E <= lambda)
    if(length(I2) == 0) I2 <- 0

    K_Ep <- diag(length(Ep)) - 1/(gamma-1) * H_hat_EpEp_inv %*% I1
    K_inv <- solve(K_Ep)

    A1 <- - diag(sign(theta_hat[Ep]), nrow = length(Ep)) %*% K_inv
    b1 <- - lambda * diag(sign(theta_hat[Ep]), nrow = length(Ep)) %*% K_inv %*% H_hat_EpEp_inv %*% (I2 + gamma/(gamma-1) * I1) %*% subgrad[Ep]

    b_term <- lambda * Hn[nE,Ep, drop=F] %*% K_inv %*% H_hat_EpEp_inv %*% (I2 + gamma/(gamma-1) * I1) %*% subgrad[Ep]
    A00 <- - Hn[nE,Ep, drop=F] %*% K_inv + Jn[nE,Ep, drop=F] %*% Jn_EpEp_inv %*% Hn[Ep,Ep, drop=F]
    b00 <- lambda - theta_perp[nE] - b_term
    A01 <- - A00
    b01 <- lambda + theta_perp[nE] + b_term

    E0 <- which(V_hat_E > lambda*gamma)
    E2 <- which(V_hat_E <= lambda)
    E1 <- which(V_hat_E > lambda & V_hat_E <= lambda * gamma)

    M1 <- diag(length(E1)) - H_hat_EpEp_inv[E1,E1, drop=F] / (gamma - 1)

    if(length(E1)!=0){
      M1_inv <- solve(M1)
    } else{ M1_inv <- matrix(NA, length(E1), length(E1)) }

    K0 <- cbind(diag(length(E0)), H_hat_EpEp_inv[E0,E1, drop=F] %*% M1_inv / (gamma - 1), matrix(0, length(E0), length(E2)))
    k0 <- (lambda * gamma / (gamma - 1) * H_hat_EpEp_inv[E0,E1, drop=F] %*% subgrad[Ep][E1]
           + lambda * H_hat_EpEp_inv[E0,E2, drop=F] %*% subgrad[Ep][E2]
           + lambda * gamma / (gamma - 1)^2 * H_hat_EpEp_inv[E0,E1, drop=F] %*% M1_inv %*% H_hat_EpEp_inv[E1,E1, drop=F] %*% subgrad[Ep][E1]
           + lambda / (gamma - 1) * H_hat_EpEp_inv[E0,E1, drop=F] %*% M1_inv %*% H_hat_EpEp_inv[E1,E2, drop=F] %*% subgrad[Ep][E2])

    if(length(E1) != 0){
      K1 <- cbind(matrix(0, length(E1), length(E0)), M1_inv, matrix(0, length(E1), length(E2)))
      k1 <- (lambda * gamma / (gamma - 1) * M1_inv %*% H_hat_EpEp_inv[E1,E1, drop=F] %*% subgrad[Ep][E1]
             + lambda * M1_inv %*% H_hat_EpEp_inv[E1,E2, drop=F] %*% subgrad[Ep][E2])
    } else{
      K1 <- 0
      k1 <- 0
    }
    K2 <- cbind(matrix(0, length(E2), length(E0)), H_hat_EpEp_inv[E2,E1, drop=F] %*% M1_inv / (gamma - 1), diag(length(E2)))
    k2 <- (lambda * gamma / (gamma - 1) * H_hat_EpEp_inv[E2,E1, drop=F] %*% subgrad[Ep][E1]
           + lambda * H_hat_EpEp_inv[E2,E2, drop=F] %*% subgrad[Ep][E2]
           + lambda * gamma / (gamma - 1)^2 * H_hat_EpEp_inv[E2,E1, drop=F] %*% M1_inv %*% H_hat_EpEp_inv[E1,E1, drop=F] %*% subgrad[Ep][E1]
           + lambda / (gamma - 1) * H_hat_EpEp_inv[E2,E1, drop=F] %*% M1_inv %*% H_hat_EpEp_inv[E1,E2, drop=F] %*% subgrad[Ep][E2])

    A00_0 <- - diag(sign(theta_hat[Ep])[E0], nrow = length(E0)) %*% K0
    b00_0 <- - gamma * lambda - diag(sign(theta_hat[Ep])[E0], nrow = length(E0)) %*% k0

    A00_10 <- diag(sign(theta_hat[Ep])[E1], nrow = length(E1)) %*% K1
    b00_10 <- gamma * lambda + diag(sign(theta_hat[Ep])[E1], nrow = length(E1)) %*% k1

    A00_11 <- - diag(sign(theta_hat[Ep])[E1], nrow = length(E1)) %*% K1
    b00_11 <- - lambda - diag(sign(theta_hat[Ep])[E1], nrow = length(E1)) %*% k1
    if(length(E1)==0){
      A00_11 <- NULL
      b00_11 <- NULL
      A00_10 <- NULL
      b00_10 <- NULL
    }
    A00_2 <- diag(sign(theta_hat[Ep])[E2], nrow = length(E2)) %*% K2
    b00_2 <- lambda + diag(sign(theta_hat[Ep])[E2], nrow = length(E2)) %*% k2

    og_order <- order(c(E0, E1, E2))

    A00_cond <- rbind(rbind(A00_0, A00_10, A00_2)[og_order,og_order], rbind(A00_0, A00_11, A00_2)[og_order,og_order])
    b00_cond <- c(c(b00_0, b00_10, b00_2)[og_order], c(b00_0, b00_11, b00_2)[og_order])

    A00 <- rbind(A00, A00_cond)
    b00 <- c(b00, b00_cond)

  } else if(penalty=="mcp"){

    theta_perp[nE] <- lambda*subgrad[nE] + Hn[nE,Ep, drop=F] %*% theta_hat[Ep] - Jn[nE,Ep, drop=F] %*% Jn_EpEp_inv %*% Hn[Ep,Ep] %*% theta_onestep[Ep]

    I1 <- 1 * diag(V_hat_E<=lambda*gamma)
    if(length(I1)==0) I1 <- 0

    K_inv <- solve(diag(length(Ep)) - (1 / gamma) * H_hat_EpEp_inv %*% I1)

    A1 <- - diag(sign(theta_hat[Ep]), nrow = length(Ep)) %*% K_inv
    b1 <- -lambda * diag(subgrad[Ep], nrow = length(Ep)) %*% K_inv %*% H_hat_EpEp_inv %*% I1 %*% subgrad[Ep]

    A_term1 <- Hn[nE,Ep, drop=F] %*% K_inv
    A_term2 <- Jn[nE,Ep, drop=F] %*% Jn_EpEp_inv %*% Hn[Ep,Ep]
    b_term <- lambda * Hn[nE,Ep, drop=F] %*% K_inv %*% H_hat_EpEp_inv %*% I1 %*% subgrad[Ep]

    A00 <- - A_term1 + A_term2
    b00 <- lambda - theta_perp[nE] - b_term

    A01 <- -A00
    b01 <- lambda + theta_perp[nE] + b_term

    E1 <- which(V_hat_E <= lambda * gamma)
    E0 <- which(V_hat_E > lambda * gamma)

    M1 <- diag(length(E1)) - H_hat_EpEp_inv[E1,E1] / gamma
    M1_inv <- solve(M1)

    K1 <- cbind(matrix(0, length(E1), length(E0)), M1_inv)
    k1 <- lambda * M1_inv %*% H_hat_EpEp_inv[E1,E1, drop=F] %*% subgrad[Ep][E1]

    K0 <- cbind(diag(length(E0)), H_hat_EpEp_inv[E0,E1, drop=F] %*% M1_inv / gamma)
    k0 <- (lambda * H_hat_EpEp_inv[E0,E1, drop=F] %*% subgrad[Ep][E1]
           + (lambda/gamma) * H_hat_EpEp_inv[E0,E1, drop=F] %*% M1_inv %*% H_hat_EpEp_inv[E1,E1, drop=F] %*% subgrad[Ep][E1])

    A00_0 <- - diag(sign(theta_hat[Ep])[E0], nrow=length(E0)) %*% K0
    b00_0 <- - gamma * lambda - diag(sign(theta_hat[Ep])[E0], nrow=length(E0)) %*% k0

    A00_1 <- diag(sign(theta_hat[Ep])[E1], nrow=length(E1)) %*% K1
    b00_1 <- gamma * lambda + diag(sign(theta_hat[Ep])[E1], nrow=length(E1)) %*% k1

    og_order <- order(c(E0, E1))

    A00_cond <- rbind(rbind(A00_0, A00_1)[og_order,og_order])
    b00_cond <- c(c(b00_0, b00_1)[og_order])

    A00 <- rbind(A00, A00_cond)
    b00 <- c(b00, b00_cond)
  }
  A <- rbind(A1, A00, A01)
  b <- c(b1, b00, b01)

  if(!all(A %*% theta_onestep[Ep] <= b)){
    print(subgrad[E])
    print(subgrad[nE])
    stop("Affine constraint not satisfied.")
  }
  list(A = A, b = b, theta_onestepE = theta_onestep[Ep], Sigma_E = Sigma_Ep)
}
