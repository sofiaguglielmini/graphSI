lasso_derivative <- function(x, lambda){
  lambda*sign(x)
}

elnet_derivative <- function(x, lambda, gamma){
  lambda*(gamma+(1-gamma)*abs(x))
}
scad_derivative <- function(x, lambda, gamma) {
  absx <- abs(x)
  ifelse(absx <= lambda, lambda,
                ifelse(absx > gamma * lambda, 0, (gamma * lambda - absx) / (gamma - 1)))
}

mcp_derivative <- function(x, lambda, gamma){
  absx <- abs(x)
  ifelse(absx <= gamma*lambda, lambda-absx/gamma, 0)
}

penalty_derivative <- function(x, penalty, lambda, gamma = NULL){
  if(penalty == "lasso"){
    return(lasso_derivative(x, lambda))
  } else if(penalty == "elastic net"){
    return(elnet_derivative(x, lambda, gamma))
  } else if(penalty == "SCAD"){
    return(scad_derivative(x, lambda, gamma))
  } else if(penalty == "MCP"){
    return(mcp_derivative(x, lambda, gamma))
  } else{
    stop("Penalty not recognized.")
  }
}
