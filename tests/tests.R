# test that graphSelect and graphInference work
set.seed(1)
n <- 50
p <- 50
library(MASS)
Sigma <- diag(p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    if(sample(c(0,0,0,0,1),1)==0) next
    Sigma[i,j] <- 0.5^(abs(i-j))
    Sigma[j,i] <- Sigma[i,j]
  }
}
Theta <- round(solve(Sigma),3)
pval <- NA
for(jsim in 1:1000){
  data <- mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
  selected <- graphSelect(data, penalty="lasso", lambda=NULL, data.splitting=F)
  j <- sample(1:nrow(selected$selected.indices), 3)
  true.value <- Theta[selected$selected.indices[j,1],selected$selected.indices[j,2]]
  # if(true.value!=0) next
  inference <- graphInference(data, selected, j, nullvalue=0, sandwich.variance=FALSE, alpha=0.05, seed=1)
  pval[jsim] <- inference$p_value
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
