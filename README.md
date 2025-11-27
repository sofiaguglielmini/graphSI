# graphSI
GraphSI provides tools for selection and post-selection inference in undirected graphical models. 

## Features
- **Graph selection**: Estimate the structure of an undirected graphical model from a data matrix using Lasso or Elastic Net
- **Data-splitting inference**: Split the data to separate selection and inference, overcoming selection bias.
- **Polyhedral selective inference**: Compute asymptotic post-selection inference for the selected edges, using the entire data for selection and inference and conditioning on selection

## Installation
```r
# Install from GitHub
remotes::install_github("yourusername/GraphSI")
```

## Usage
```r
library(MASS)
set.seed(1)
n <- 100
p <- 10
Sigma <- diag(p)
data <- mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
# Selection
selected <- graphSelect(data, penalty="lasso", data.splitting=FALSE, seed=1)
selected
# Inference on the j-th selected edge with polyhedral selective inference:
inference <- graphInference(data, selected, j, nullvalue=0, sandwich.variance=FALSE, alpha=0.05, seed=1)
inference
```

## References

### Polyhedral selective inference
- Guglielmini, Sofia, and Gerda Claeskens. "Asymptotic post-selection inference for regularized graphical models." Statistics and Computing 35.2 (2025): 36.
- Lee, Jason D., et al. "Exact post-selection inference, with application to the lasso." (2016): 907-927.

### Graphical lasso selection
- Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. *Biostatistics*, 9(3), 432–441.
- Yuan, Ming, and Yi Lin. "Model selection and estimation in the Gaussian graphical model." Biometrika 94.1 (2007): 19-35.

### Graphical elastic net selection
- Kovács, Solt, et al. "Graphical elastic net and target matrices: Fast algorithms and software for sparse precision matrix estimation." arXiv preprint arXiv:2101.02148 (2021).

