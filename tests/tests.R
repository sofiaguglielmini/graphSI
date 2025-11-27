# test that graphSelect and graphInference work
set.seed(1)
n <- 100
p <- 10
library(MASS)
Sigma <- diag(p)
for(i in 1:(p-1)){
  if(sample(c(TRUE,FALSE),1)) next
  Sigma[i, i+1] <- 0.5
  Sigma[i+1, i] <- 0.5
}
Theta <- solve(Sigma)
data <- mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
selected <- graphSelect(data, penalty="lasso", lambda=NULL, data.splitting=F, seed=1, penalize.diagonal = F)
j <- 4 # j-th selected edge in selected$selected.indices
Theta[selected$selected.indices[j,1],selected$selected.indices[j,2]]
inference <- graphInference(data, selected, j, nullvalue=0, sandwich.variance=FALSE, alpha=0.05, seed=1)
print(inference)

