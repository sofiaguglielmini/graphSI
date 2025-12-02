# test that graphSelect and graphInference work
set.seed(1)
n <- 100
p <- 5
library(MASS)
Sigma <- diag(p)
# dense Sigma
for(i in 1:(p-1)){
  for(j in (i+1):p){
    Sigma[i,j] <- 0.9^(abs(i-j)+p)
    Sigma[j,i] <- Sigma[i,j]
  }
}
Theta <- round(solve(Sigma),3)
Theta
pval <- NA
devtools::load_all()
for(jsim in 1:500){
  data <- mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
  selected <- graphSelect(data, penalty="lasso", lambda=NULL, data.splitting=F, penalize.diagonal = F)
  j <- 1:nrow(selected$selected.indices)
  true.value <- Theta[selected$selected.indices[j,1],selected$selected.indices[j,2]]
  if(true.value!=0) next

  inference <- graphInference(data, selected, j, nullvalue=0, sandwich.variance=FALSE, alpha=0.05, seed=1)
  inference$inference
  pval[jsim] <- inference[[1]]$p_value
}
hist(pval)
mean(pval<0.05, na.rm=TRUE)







# print(subgradnE == subgrad[nE])
# print(cbind(abs(subgrad[nE]) <= 1))
# print(cbind(abs((theta_perp[nE]-(Hn[nE,E, drop=F] %*% theta_hat[E] - Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E]))/(lambda*gamma)) <= 1))
# print((theta_perp[nE]-(Hn[nE,E, drop=F] %*% theta_hat[E] - Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E]))/(lambda*gamma) <= 1)
# print(-(theta_perp[nE]-(Hn[nE,E, drop=F] %*% theta_hat[E] - Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E]))/(lambda*gamma) <= 1)

# print(theta_perp[nE] - Hn[nE,E, drop=F] %*% theta_hat[E] + Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E] <= lambda*gamma)
# print(-theta_perp[nE] + Hn[nE,E, drop=F] %*% theta_hat[E] - Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E] <= lambda*gamma)

# print(theta_perp[nE] - Hn[nE,E, drop=F] %*% (theta_onestep[E] + H_hat_EE_inv%*%(-DsubgradE%*%rhoprimeE)) + Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E] <= lambda*gamma)
# print(-theta_perp[nE] + Hn[nE,E, drop=F] %*% (theta_onestep[E] + H_hat_EE_inv%*%(-DsubgradE%*%rhoprimeE)) - Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E] <= lambda*gamma)

# print(theta_perp[nE] - Hn[nE,E, drop=F] %*% theta_onestep[E] - Hn[nE,E, drop=F] %*%H_hat_EE_inv%*%(-DsubgradE%*%rhoprimeE) + Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E] <= lambda*gamma)
# print(-theta_perp[nE] + Hn[nE,E, drop=F] %*% theta_onestep[E] + Hn[nE,E, drop=F] %*%H_hat_EE_inv%*%(-DsubgradE%*%rhoprimeE) - Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E] %*% theta_onestep[E] <= lambda*gamma)

# print(theta_perp[nE] + (- Hn[nE,E, drop=F]+ Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E]) %*% theta_onestep[E] - Hn[nE,E, drop=F] %*%H_hat_EE_inv%*%(-DsubgradE%*%rhoprimeE) <= lambda*gamma)
# print(-theta_perp[nE] + (Hn[nE,E, drop=F]- Jn[nE,E, drop=F] %*% Jn_EE_inv %*% Hn[E,E]) %*% theta_onestep[E] + Hn[nE,E, drop=F] %*%H_hat_EE_inv%*%(-DsubgradE%*%rhoprimeE) <= lambda*gamma)

# print( - A_perp %*% theta_onestep[E] <= -theta_perp[nE] + lambda*gamma + Hn[nE,E, drop=F] %*%H_hat_EE_inv%*%(-DsubgradE%*%rhoprimeE))
# print( A_perp %*% theta_onestep[E] <= theta_perp[nE] + lambda*gamma - Hn[nE,E, drop=F] %*%H_hat_EE_inv%*%(-DsubgradE%*%rhoprimeE))f
