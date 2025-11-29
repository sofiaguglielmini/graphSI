# graphSI: Selection and Inference in graphical models
GraphSI provides tools for selection and post-selection inference in undirected graphical models. 

## Features
- **Graph selection**: Estimate the structure of an undirected graphical model from a data matrix using the graphical lasso (Friedman et al. 2008), graphical elastic net (Kovács et al. 2021), or the SCAD (Fan et al. 2009) or MCP (Zhang et al. 2010) penalties. The selection with the non-convex penalties SCAD and MCP uses local linear approximations (GGMncv; Williams, 2020).
- **Data-splitting inference**: Split the data to separate selection and inference, overcoming selection bias.
- **Polyhedral selective inference**: Compute asymptotic post-selection inference for the selected edges, using the entire data for selection and inference and conditioning on selection (Guglielmini and Claeskens, 2025), using the polyhedral lemma (Lee et al. 2016).

## Installation
```r
remotes::install_github("sofiaguglielmini/graphSI")
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
- Fan, Jianqing, Yang Feng, and Yichao Wu. "Network exploration via the adaptive LASSO and SCAD penalties." The annals of applied statistics 3.2 (2009): 521.
- Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. "Sparse inverse covariance estimation with the graphical lasso." Biostatistics 9.3 (2008): 432-441.
- Guglielmini, Sofia, and Gerda Claeskens. "Asymptotic post-selection inference for regularized graphical models." Statistics and Computing 35.2 (2025): 36.
- Kovács, Solt, et al. "Graphical elastic net and target matrices: Fast algorithms and software for sparse precision matrix estimation." arXiv preprint arXiv:2101.02148 (2021).
- Lee, Jason D., et al. "Exact post-selection inference, with application to the lasso." (2016): 907-927.
- Williams, Donald R. "Beyond lasso: A survey of nonconvex regularization in Gaussian graphical models." (2020).
- Yuan, Ming, and Yi Lin. "Model selection and estimation in the Gaussian graphical model." Biometrika 94.1 (2007): 19-35.
- Zhang, Cun-Hui. "Nearly unbiased variable selection under minimax concave penalty." (2010): 894-942.

