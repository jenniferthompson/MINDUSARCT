################################################################################
## Functions in this script are written by Qi Liu in conjunction with Liu et al:
## Liu Q, Shepherd BE, Li C, Harrell FE. Modeling continuous response variables
##   using ordinal regression. Statistics in Medicine 2017; 36: 4316-4335.
##
## At the time of writing, source code can be found at
##   http://biostat.mc.vanderbilt.edu/wiki/pub/Main/ArchivedAnalyses/orm_code.r
## My only changes are to formatting/white space.
################################################################################

##-- Estimate conditional mean and its standard error for orm model ------------
mean.orm <- function(mod, new.data, se=TRUE){
  if(is.null(mod$yunique)) {
    stop("Need to set x=TRUE and y=TRUE for orm") 
  } else{
    order.y <- mod$yunique
    n.alpha <- length(order.y)-1
    xb <- as.matrix(new.data)%*%matrix(coef(mod)[colnames(new.data)])
    m.alpha <- mod$coef[1:n.alpha]
    lb <- t(outer(m.alpha, xb, "+")[,,1])
    m.s <- mod$trans$cumprob(lb)
    m.f <-
      t(apply(m.s, 1, FUN = function(x) c(1,x[1:n.alpha]) - c(x[1:n.alpha], 0)))
    m.mean <- apply(m.f, 1, FUN = function(x) sum(x*order.y))
    
    if(se){
      if(mod$family == "logistic"){
        mod$trans$deriv <- function(x) exp(-x)/(1+exp(-x))^2
      }
      
      dmean.dalpha <-
        t(apply(mod$trans$deriv(lb),
                1,
                FUN = function(x){
                  x*(order.y[2:length(order.y)] - order.y[1:n.alpha])
                }
        ))
      dmean.dbeta <- apply(dmean.dalpha, 1, sum)*as.matrix(new.data)
      dmean.dtheta <- cbind(dmean.dalpha, dmean.dbeta)   
      mean.var <-diag(dmean.dtheta%*%solve(mod$info.matrix)%*%t(dmean.dtheta))
      mean.se <- sqrt(mean.var)   
      result <- cbind(m.mean, mean.se)
      ci <- t(apply(result,
                    1,
                    FUN = function(x){
                      c(x[1]- qnorm(0.975)*x[2], x[1]+ qnorm(0.975)*x[2])
                    }))
      result <- cbind(result, ci)
      colnames(result) <- c("est", "se", "lb", "ub")
    } else{
      result <- matrix(m.mean)
      colnames(result) <- c("est")
    }
    
    return(result)
  } 
}

## -- Estimate conditional quantiles and confidence intervals for orm model ----
quantile.orm <- function(mod, new.data, probs = 0.5, se = TRUE){
  
  quantile <- matrix(NA, nrow = dim(new.data)[1], ncol = length(probs))
  order.y <- mod$yunique
  #n.alpha <- length(order.y)-1
  xb <- as.matrix(new.data) %*% matrix(coef(mod)[colnames(new.data)])
  alpha <- mod$coef[1:(length(unique(order.y)) - 1)]
  lb <- t(outer(alpha, xb, "+")[, , 1])
  m.cdf <- 1 - mod$trans$cumprob(lb)
  m.cdf <- cbind(0, m.cdf, 1)
  for(i in 1:length(probs)){
    try({
      index.1 <-
        apply(m.cdf, 1, FUN = function(x){ max(which(x <= probs[i]))[1]} )
      index.2 <-
        apply(m.cdf, 1, FUN = function(x){ min(which(x >= probs[i]))[1]} )
      
      index.y1 <- ifelse(index.1 > length(order.y), Inf, order.y[index.1])
      index.y2 <- ifelse(index.2 > length(order.y), Inf, order.y[index.2])
      
      index.y1.cdf <-
        ifelse(index.1 == 0, 0, m.cdf[cbind(1:dim(new.data)[1], index.1)])
      
      index.y2.cdf <- ifelse(index.2 > length(order.y),
                             1,
                             m.cdf[cbind(1:dim(new.data)[1], index.2)])
      
      quantile[, i] <-
        ifelse(
          index.1 == index.2,
          index.y1,
          (index.y2 - index.y1)/(index.y2.cdf - index.y1.cdf)*(probs[i] - index.y1.cdf) + index.y1)
      quantile[, i] <-
        ifelse(is.infinite(quantile[,i]), max(order.y), quantile[, i])
    })
  }
  result <- quantile
  
  if(se){
    if(mod$family == "logistic"){
      mod$trans$deriv <- function(x) exp(-x)/(1 + exp(-x))^2
    }
    
    quantile.lb <- quantile.ub <-
      matrix(NA, nrow = dim(new.data)[1], ncol = length(probs))
    
    lb.se <- matrix(NA, ncol = dim(lb)[2], nrow = dim(new.data)[1])
    var <- as.matrix(solve(mod$info.matrix))
    
    for(i in 1:dim(lb)[2]){
      var.i <- var[c(i, which(names(coef(mod)) %in% colnames(new.data))), 
                   c(i, which(names(coef(mod)) %in% colnames(new.data)))]
      
      dcdf.dtheta <- cbind(-mod$trans$deriv(lb[,i]),  
                           -mod$trans$deriv(lb[,i]) * as.matrix(new.data))
      dlb.dtheta <- as.matrix(cbind(1, new.data))
      lb.se[, i] <- sqrt(diag(dlb.dtheta %*% var.i %*% t(dlb.dtheta)))
    }
    
    ci.lb <- sapply(
      1:dim(lb)[2],
      FUN = function(i){
        1- mod$trans$cumprob(lb[, i] + qnorm(0.975) * lb.se[, i])
      }
    )
    ci.ub <- sapply(
      1:dim(lb)[2],
      FUN = function(i){
        1- mod$trans$cumprob(lb[, i] -qnorm(0.975) * lb.se[, i])
      }
    )
    ci.lb <- matrix(ci.lb, nrow = dim(new.data)[1])
    ci.ub <- matrix(ci.ub, nrow = dim(new.data)[1])
    
    ci.lb <- cbind(0, ci.lb, 1)
    ci.ub <- cbind(0, ci.ub, 1)
    
    for(i in 1: length(probs)){
      try({
        index.1 <- apply(ci.lb, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
        index.2 <- apply(ci.lb, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
        
        index.y1 <- ifelse(index.1>length(order.y), Inf, order.y[index.1])
        index.y2 <- ifelse(index.2>length(order.y), Inf, order.y[index.2])
        
        index.y1.cdf <-
          ifelse(index.1==0, 0, ci.lb[cbind(1:dim(new.data)[1], index.1)])
        
        index.y2.cdf <- ifelse(
          index.2 > length(order.y), 1, ci.lb[cbind(1:dim(new.data)[1], index.2)]
        )
        
        quantile.lb[,i] <- ifelse(
          index.1 == index.2,
          index.y1,
          (index.y2 - index.y1)/(index.y2.cdf - index.y1.cdf) * (probs[i] - index.y1.cdf) + index.y1
        )
        quantile.lb[, i] <-
          ifelse(is.infinite(quantile.lb[, i]), max(order.y), quantile.lb[, i])
        
        index.1 <-
          apply(ci.ub, 1, FUN=function(x){ max(which(x <= probs[i]))[1]} )
        index.2 <-
          apply(ci.ub, 1, FUN=function(x){ min(which(x >= probs[i]))[1]} )
        
        index.y1 <- ifelse(index.1 > length(order.y), Inf, order.y[index.1])
        index.y2 <- ifelse(index.2 > length(order.y), Inf, order.y[index.2])
        
        index.y1.cdf <- ifelse(
          index.1 == 0, 0, ci.ub[cbind(1:dim(new.data)[1], index.1)]
        )
        
        index.y2.cdf <- ifelse(
          index.2 > length(order.y), 1, ci.ub[cbind(1:dim(new.data)[1], index.2)]
        )
        
        quantile.ub[,i] <- ifelse(
          index.1 == index.2,
          index.y1, 
          (index.y2 - index.y1)/(index.y2.cdf - index.y1.cdf) * (probs[i] - index.y1.cdf) + index.y1
        ) 
        quantile.ub[, i] <-
          ifelse(is.infinite(quantile.ub[, i]), max(order.y), quantile.ub[, i])
      })
    }
    
    result <- list(quantile = quantile,
                   lb = quantile.ub,
                   ub = quantile.lb)
  }
  
  return(result)
}

## Wrapper to put these in a data.frame
quantile_orm_df <- function(
  ## arguments passed to quantile.orm()
  mod, new.data = pred_df,
  ## levels of treatment to append to final df
  trt_levels = levels(model_df$trt),
  ...
){
  quantile.orm(mod = mod, new.data = new.data, ...) %>%
    ## Each quantity initially stored as a matrix, in case we want multiple
    ##   quantiles (we don't); need to be numeric
    map_dfc(as.numeric) %>%
    bind_cols(data.frame(trt = trt_levels), .)
}

## -- Estimate conditional CDF and its standard error for orm models -----------
cdf.orm <- function(mod, new.data, at.y = 0, se = TRUE){
  if(is.null(mod$yunique)) {
    stop("Need to set x=TRUE and y=TRUE for orm") 
  } else{
    order.y <- mod$yunique
    xb <- as.matrix(new.data)%*%matrix(coef(mod)[colnames(new.data)])
    
    index <- sapply(
      at.y,
      FUN = function(x){
        if(x < min(order.y)[1]) result <- Inf 
        else if(x == min(order.y)[1]) result <- 1
        else if(x >= max(order.y)[1]) result <- -Inf
        else which(order.y >= x)[1]-1
      }
    )
    
    m.alpha <- mod$coef[index]
    m.alpha <- ifelse(is.infinite(index), index, m.alpha)
    if(length(at.y) == 1){
      lb <- as.matrix(outer(m.alpha, xb, "+")[,,1])
    } else lb <- t(outer(m.alpha, xb, "+")[,,1])
    m.cdf <- 1- mod$trans$cumprob(lb)
    
    if(se){
      if(mod$family == "logistic") mod$trans$deriv <- function(x) exp(-x)/(1+exp(-x))^2
      cdf.se <- matrix(NA, ncol = length(at.y), nrow = dim(new.data)[1])
      lb.se <- matrix(NA, ncol = length(at.y), nrow = dim(new.data)[1])
      
      var <- as.matrix(solve(mod$info.matrix))
      
      for(i in 1:length(at.y)) {
        var.i <- var[c(index[i], which(names(coef(mod)) %in% colnames(new.data))), 
                     c(index[i], which(names(coef(mod)) %in% colnames(new.data)))]
        dcdf.dtheta <- cbind(-mod$trans$deriv(lb[,i]),  
                             -mod$trans$deriv(lb[,i])*as.matrix(new.data) )
        dlb.dtheta <- as.matrix(cbind(1, new.data))
        cdf.se[,i] <- sqrt(diag(dcdf.dtheta %*% var.i%*% t(dcdf.dtheta)))
        lb.se[, i] <- sqrt(diag(dlb.dtheta%*%var.i%*% t(dlb.dtheta)))
      }
      
      ci.lb <- sapply(
        1:length(at.y),
        FUN=function(i) { 1- mod$trans$cumprob(lb[, i] +qnorm(0.975)*lb.se[, i])}
      )
      ci.ub <- sapply(
        1:length(at.y),
        FUN=function(i) { 1- mod$trans$cumprob(lb[, i] -qnorm(0.975)*lb.se[, i])}
      )
      
      result <- list(est = m.cdf,
                     se = cdf.se,
                     lb = ci.lb,
                     ub = ci.ub)
    } else{
      result <- list(est = m.cdf)
    }
    
    return(result)
  } 
}
