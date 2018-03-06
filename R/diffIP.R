#' @title diffIP
#' @param x The RNADMethyl data list
#' @param Covariates Matrix of covariates to include in the PoissonGamma test. Default is NULL.
#' @param p.adjBy The method used to adjust for the p value. Default is 'fdr'. Can be "bonferroni","BY","fdr".
#' @param exclude Names of samples to exclude (as outlier) in the test. Default is NULL
#' @param maxPsi The cutoff value for random effect parameter Psi estimation
#' @export
diffIP <- function(
  x, 
  Covariates = NULL, # covariates
  plotPvalue = TRUE,
  p.adjBy = "fdr", # the method to calculate adjusted p_value
  exclude = NULL,
  maxPsi = 100
){
  #### start run PoissonGamma test
  if( is.null(Covariates) ){ # run simple PoissonGamma
    
    if(is.null(exclude) ){ ## running default
      allY <- x$ip_adjExpr_filtered
      X <- x$X
    }else if( all( exclude %in% x$samplenames) ){
      exc.id <- match(exclude,x$samplenames)
      allY <- x$ip_adjExpr_filtered[,-c(exc.id)]
      X <- x$X[-c(exc.id)]
      cat(paste0("Sample ",exclude, " has been removed from test...\n"))
    }else{
      stop("Samples to exclude is not in samplenames...")
    }
    
    tmp <- rafalib::as.fumeric(X) - 1
    names(tmp) <- X
    cat("The categorical variable has been converted:\n")
    print(tmp)
    X <- rafalib::as.fumeric(X) - 1 # convert categorical variable into numerical variable.
    
    
    psi <- 10
    
    print("running PoissonGamma test at single beta mode...")
    pb <- txtProgressBar(min = 1, max = nrow(allY), style = 3) ##creat a progress bar to track loop progress
    all.est <- NULL
    all.id <- NULL
    for(kk in 1:nrow(allY)){
      Y <- unlist(allY[kk, ])
      model1 <- glm(Y ~ X, family = poisson(link = 'log'))
      coef <- model1$coefficients
      mu2 <- coef[1]
      beta <- coef[2]
      est <- try(unlist(PoissonGamma::PoissionGamma(Y, X, beta, psi, mu2, gamma = 0.75, steps = 50, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) != "try-error"){
        all.est <- rbind(all.est, est)
        all.id <- c(all.id, kk)
      }
      setTxtProgressBar(pb, kk) # update progress bar
    }
    
    rownames(all.est) <- rownames(allY)[all.id]
    if(plotPvalue){ # plot the distribution of p values
      hist(all.est[,"p_value"],main = "P Values", xlab = "p_value")
    }
    
    ##calcualte adjusted p value
    padj <- p.adjust(all.est[,"p_value"],method = p.adjBy )
    all.est <- cbind( all.est, padj )
    
    x <- c(x,list('all.est'=all.est))
    return(x)
    
  } else if( !is.null(Covariates) & dim( cbind(X,Covariates) )[1] == length(X)  ){ # include covariates
    
    if(is.null(exclude) ){ ## running default
      allY <- x$ip_adjExpr_filtered
      X <- x$X
    }else if( all( exclude %in% x$samplenames) ){
      exc.id <- match(exclude,x$samplenames)
      allY <- x$ip_adjExpr_filtered[,-c(exc.id)]
      X <- x$X[-c(exc.id)]
      Covariates <- Covariates[-c(exc.id),]
      cat(paste0("Sample ",exclude, " has been removed from test...\n"))
    }else{
      stop("Samples to exclude is not in samplenames...")
    }
    
    tmp <- rafalib::as.fumeric(X) - 1
    names(tmp) <- X
    print("The categorical variable has been converted:")
    print(tmp)
    X1 <- rafalib::as.fumeric(X) - 1 # convert categorical variable into numerical variable.
    
    ## organize covariates
    if(!is.numeric(Covariates)){stop("Please convert covariates into numerical variables")}

    X.all <- cbind(X1,Covariates) # new design matrix
    colnames(X.all) <- paste("X",1:ncol(X.all),sep = "")
    
    #design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste("X.all[", rownames(X.all), sep = "",collapse = " + ")) )
    design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste("X.all[,", 1:ncol(X.all),"]", sep = "",collapse = " + ")) )
    ## Run multi-beta PoissonGamma
    psi <- 10
    
    cat("running PoissonGamma test at multi-beta mode...\n")
    pb <- txtProgressBar(min = 1, max = nrow(allY), style = 3) ##creat a progress bar to track loop progress
    all1 <- NULL
    all.id <- NULL
    for(kk in 1:nrow(allY)){
      Y <- unlist(allY[kk, ] )
      aa <- unlist(summary( lm( design.multiBeta ) )$coefficients[, 1])
      mu2 <- aa[1]
      beta <- aa[2:(ncol(X.all)+1 )]
      est <- try(unlist(PoissonGamma::PoissionGamma_multiple_beta(Y, X.all, beta, psi, mu2, gamma = 0.25, steps = 10, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) != "try-error"){
        all1 <- rbind(all1, est)
        all.id <- c(all.id, kk)
      }
      setTxtProgressBar(pb, kk) # update progress bar
    }
    
    rownames(all1) <- rownames(allY)[all.id] ## assign window names to test statistics
    if(plotPvalue){ # plot the distribution of p values
      hist(all1[,"p_value3"],main = "P Values", xlab = "p_value")
    }
    ##calcualte adjusted p value
    padj <- p.adjust(all1[,"p_value3"],method = p.adjBy )
    all1 <- cbind( all1, padj )
    
    x <- c(x,list('all.est'=all1))
    return(x)
    
  } else{ ## if Covariate is not null and not the same dimension of X, throw a error
    stop("Covariates need to be the same dimension as X")
  }
  
}

#' @title diffIP_parallel
#' @param x The RNADMethyl data list
#' @param Covariates Matrix of covariates to include in the PoissonGamma test. Default is NULL.
#' @param p.adjBy The method used to adjust for the p value. Default is 'fdr'. Can be "bonferroni","BY","fdr".
#' @param exclude Names of samples to exclude (as outlier) in the test. Default is NULL
#' @param maxPsi The cutoff value for random effect parameter Psi estimation
#' @param thread The number threads to run in parallel
#' @export
diffIP_parallel <- function(x, 
                            Covariates = NULL, # covariates
                            plotPvalue = TRUE,
                            p.adjBy = "fdr", # the method to calculate adjusted p_value
                            exclude = NULL,
                            maxPsi = 100,
                            thread = 4){
  
  #### start run PoissonGamma test
  if( is.null(Covariates) ){ # run simple PoissonGamma
    
    if(is.null(exclude) ){ ## running default
      allY <- x$ip_adjExpr_filtered
      X <- x$X
    }else if( all( exclude %in% x$samplenames) ){
      exc.id <- match(exclude,x$samplenames)
      allY <- x$ip_adjExpr_filtered[,-c(exc.id)]
      X <- x$X[-c(exc.id)]
      cat(paste0("Sample ",exclude, " has been removed from test...\n"))
    }else{
      stop("Samples to exclude is not in samplenames...")
    }
    
    tmp <- rafalib::as.fumeric(X) - 1
    names(tmp) <- X
    cat("The categorical variable has been converted:\n")
    print(tmp)
    X <- rafalib::as.fumeric(X) - 1 # convert categorical variable into numerical variable.
    
    
    psi <- 10
    
    cat("running PoissonGamma test at single beta mode...\n")
    start_time <- Sys.time()
    ## register cluster for hyperthread computing
    doParallel::registerDoParallel(cores=thread)
    cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
    cat(paste("Using",getDoParWorkers(),"thread(s) to run PoissonGamma test...\n"))
    
    error.id <- NULL
    all.est <- foreach(kk = 1:nrow(allY), .combine = rbind, .errorhandling = "remove") %dopar% {
      Y <- unlist(allY[kk, ])
      model1 <- glm(Y ~ X, family = poisson(link = 'log'))
      coef <- model1$coefficients
      mu2 <- coef[1]
      beta <- coef[2]
      est <- try(unlist(PoissonGamma::PoissionGamma(Y, X, beta, psi, mu2, gamma = 0.75, steps = 50, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) == "try-error"){
        error.id <- c(error.id, kk)
      }
      est
    }
    rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
    end_time <- Sys.time()
    cat(paste("Time used to run PoissonGamma test:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    
    all.id <- which(! 1:nrow(allY) %in% error.id )
    rownames(all.est) <- rownames(allY)[all.id]
    
    if(plotPvalue){ # plot the distribution of p values
      hist(all.est[,"p_value"],main = "P Values", xlab = "p_value")
    }
    
    ##calcualte adjusted p value
    padj <- p.adjust(all.est[,"p_value"],method = p.adjBy )
    all.est <- cbind( all.est, padj )
    
    x <- c(x,list('all.est'=all.est))
    return(x)
    
  } else if( !is.null(Covariates) & dim( cbind(X,Covariates) )[1] == length(X)  ){ # include covariates
    
    if(is.null(exclude) ){ ## running default
      allY <- x$ip_adjExpr_filtered
      X <- x$X
    }else if( all( exclude %in% x$samplenames) ){
      exc.id <- match(exclude,x$samplenames)
      allY <- x$ip_adjExpr_filtered[,-c(exc.id)]
      X <- x$X[-c(exc.id)]
      Covariates <- Covariates[-c(exc.id),]
      cat(paste0("Sample ",exclude, " has been removed from test...\n"))
    }else{
      stop("Samples to exclude is not in samplenames...")
    }
    
    tmp <- rafalib::as.fumeric(X) - 1
    names(tmp) <- X
    print("The categorical variable has been converted:")
    print(tmp)
    X1 <- rafalib::as.fumeric(X) - 1 # convert categorical variable into numerical variable.
    
    ## organize covariates
    if(!is.numeric(Covariates)){stop("Please convert covariates into numerical variables")}
    
    X.all <- cbind(X1,Covariates) # new design matrix
    colnames(X.all) <- paste("X",1:ncol(X.all),sep = "")
    
    #design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste("X.all[", rownames(X.all), sep = "",collapse = " + ")) )
    #design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste("X.all[,", 1:ncol(X.all),"]", sep = "",collapse = " + ")) )
    ## Run multi-beta PoissonGamma
    psi <- 10
    
    cat("running PoissonGamma test at multi-beta mode...\n")
    error.id <- NULL
    start_time <- Sys.time()
    ## register cluster for hyperthread computing
    doParallel::registerDoParallel(cores=thread)
    cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
    cat(paste("Using",getDoParWorkers(),"thread(s) to run multi-beta PoissonGamma test...\n"))
    all1 <- foreach(kk = 1:nrow(allY),.combine = rbind, .errorhandling = "remove") %dopar% {
      
      Y <- unlist(allY[kk, ] )
      y <- log(Y+1)
      lmY <- as.data.frame(cbind(y,X.all) )
      aa <- unlist(summary( lm( y~.,data=lmY ) )$coefficients[, 1])
      mu2 <- aa[1]
      beta <- aa[2:(ncol(X.all)+1 )]
      est <- try(unlist(PoissonGamma::PoissionGamma_multiple_beta(Y, X.all, beta, psi, mu2, gamma = 0.25, steps = 10, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) == "try-error"){
        error.id <- c(error.id, kk)
      }
      est
    }
    rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
    end_time <- Sys.time()
    cat(paste("Time used to run multi-beta PoissonGamma test:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    all.id <- which(! 1:nrow(all1) %in% error.id )
    rownames(all1) <- rownames(allY)[all.id] ## assign window names to test statistics
    
    if(plotPvalue){ # plot the distribution of p values
      hist(all1[,"p_value3"],main = "P Values", xlab = "p_value")
    }
    ##calcualte adjusted p value
    padj <- p.adjust(all1[,"p_value3"],method = p.adjBy )
    all1 <- cbind( all1, padj )
    
    x <- c(x,list('all.est'=all1))
    return(x)
    
  } else{ ## if Covariate is not null and not the same dimension of X, throw a error
    stop("Covariates need to be the same dimension as X")
  }
  
}

