# graphSI: Selection and Inference in Graphical Models
The aim of the graphSI package is to provide tools for selection and post-selection inference in Gaussian graphical models. 

Once a model has been selected, a question of crucial importance is often that of carrying out statistical inference and making statements on the variability of the estimators. In his vote of thanks for Tibshirani (1996), Bühlmann (2010) 
> Suggest[s] that we interpret the second 's' in lasso as 'screening' rather than 'selection'. Once we have the screening property, the task is to remove the false positive selections. [...] The issue of assigning uncertainty and variability in high dimensional statistical inference deserves further research. For example, questions about power are largely unanswered.

In an analogous way, in the context of regularized Gaussian graphical models, Williams (2020) states how
> Researchers always want to do more than detect nonzero relations in GGMs. For example, to determine which edges are the strongest or to rule relations out of the 'network' (i.e., conditional independence), each of which requires more than merely mining data. In other words, statistical inference still requires a p-value or confidence interval, neither of which is straightforward to obtain after data-driven model selection. [...] [non-convex regularization] has its place, for example, to gain the first glimpse into a dependence structure or to formulate hypotheses to then test with inferential statistics.

## Installation
```r
remotes::install_github("sofiaguglielmini/graphSI")
```

## Features
- **Graph selection**: Estimate the structure of an undirected graphical model from a data matrix using the graphical lasso (Friedman et al. 2008), graphical elastic net (Kovács et al. 2021), or the SCAD (Fan et al. 2009) or MCP (Zhang et al. 2010) penalties. The selection with the non-convex penalties SCAD and MCP uses local linear approximations (GGMncv; Williams, 2020).
- - **Polyhedral selective inference**: Compute asymptotic post-selection inference for the selected edges, using the entire data for selection and inference and conditioning on selection (Guglielmini and Claeskens, 2025), using the polyhedral lemma (Lee et al. 2016).
- **Data-splitting inference**: Split the data to separate selection and inference.

## Functions

`graphSelect()`
Implements model selection for Gaussian graphical models using regularized estimation.  Supports lasso, elastic net, SCAD, and MCP penalties, with optional diagonal penalization and optional data splitting. Returns an object containing the adjacency matrix and selected edge indices.
For valid selective inference, `\lambda` must be chosen independently of the data used for selection.  

`graphInference()`
Performs inference for a selected edge using either the polyhedral method or data splitting. The method is determined by the `graphSelect` object passed through the `selected` argument. 
Returns a p–value and confidence interval for the chosen elements of the precision matrix.

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

# View the selected edges
selected

# Inference on the j-th selected edge with polyhedral selective inference:
inference <- graphInference(data, selected, j, nullvalue=0, sandwich.variance=FALSE, alpha=0.05, seed=1)
inference
```

## References
Bühlmann, Peter. "Proposing the vote of thanks: Regression shrinkage and selection via the Lasso: a retrospective by Robert Tibshirani." 2010.

Fan, Jianqing, Yang Feng, and Yichao Wu. "Network exploration via the adaptive LASSO and SCAD penalties." _The Annals of Applied Statistics_ 3.2 (2009): 521.

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. "Sparse inverse covariance estimation with the graphical lasso." _Biostatistics_ 9.3 (2008): 432-441.

Guglielmini, Sofia, and Gerda Claeskens. "Asymptotic post-selection inference for regularized graphical models." _Statistics and Computing_ 35.2 (2025): 36.

Kovács, Solt, et al. "Graphical elastic net and target matrices: Fast algorithms and software for sparse precision matrix estimation." _arXiv preprint_ arXiv:2101.02148 (2021).

Lee, Jason D., et al. "Exact post-selection inference, with application to the lasso." Annals of Statistics 44.3 (2016): 907-927.

Tibshirani, Robert. "Regression shrinkage and selection via the lasso." _Journal of the Royal Statistical Society Series B: Statistical Methodology_ 58.1 (1996): 267-288.

Williams, Donald R. 2020. "Beyond Lasso: A Survey of Nonconvex Regularization in Gaussian Graphical Models." _PsyArXiv_.

Yuan, Ming, and Yi Lin. "Model selection and estimation in the Gaussian graphical model." _Biometrika_ 94.1 (2007): 19-35.

Zhang, Cun Hui. "Nearly unbiased variable selection under minimax concave penalty." _Annals of Statistics_ 38.2 (2010): 894-942.

