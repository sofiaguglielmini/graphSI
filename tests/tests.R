# test that graphSelect and graphInference work
set.seed(1)
n <- 50
p <- 100
library(MASS)
Sigma <- diag(p)
# for(i in 1:(p-1)){
#   Sigma[i, i+1] <- 0.5
#   Sigma[i+1, i] <- 0.5
# }
Theta <- solve(Sigma)
data <- mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
selected <- graphSelect(data, penalty="elastic net", lambda=NULL, data.splitting=FALSE, seed=1, penalize.diagonal = T)
selected
j <- 4 # j-th selected edge in selected$selected.indices
inference <- graphInference(data, selected, j, nullvalue=0, sandwich.variance=FALSE, alpha=0.05, seed=1)
print(inference)
