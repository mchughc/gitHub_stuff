estVarComp <- function(scanAnnot, covMatList, outcome, covar.vec = NULL, scan.include = NULL, start = NULL, family = gaussian, group.var = NULL, AIREML.tol = 1e-6, maxIter = 100, dropZeros = TRUE, verbose = TRUE){
    
    # check family
    if(is.character(family)){
        family <- get(family)
    }
    if(is.function(family)){
        family <- family()
    }
    if(is.null(family$family)){
        stop("'family' not recognized")
    }
    
    # create design matrices and outcome vector
    if(verbose) message("Reading in Phenotype and Covariate Data...")
    dat <- .createDesignMatrix(scanAnnot = scanAnnot, outcome = outcome, covar.vec = covar.vec, ivar.vec = NULL, scan.include = scan.include)
    # outcome vector
    Y <- dat$Y
    # design matrix
    W <- dat$W
    k <- ncol(W)
    rm(dat)
    # scans to include
    scan.include <- rownames(W)

    # sample size
    n <- length(scan.include)

    # determine unique groups and create index
    if(!is.null(group.var)){
        if(family$family != "gaussian"){
            stop("heterogeneous group residual variances can only be used when family = gaussian")
        }
        # read in group variable and subset to included samples
        group <- getVariable(scanAnnot, group.var)[getScanID(scanAnnot) %in% scan.include]
        # determine unique groups
        group.names <- as.character(unique(group))
        g <- length(group.names)
        group.idx <- list()
        for(i in 1:g){
            group.idx[[i]] <- which(group == group.names[i])
        }
    }else{
        group.names <- "E"
        g <- 1
    }
    # flag for heterogeneous residual variances
    hetResid <- (g > 1)
    
    # if covMatList is a matrix, convert to a list
    if(class(covMatList) == "matrix"){
        covMatList <- list(A = covMatList)
    }
    
    # number of covariance structure matrices
    m <- length(covMatList)
    
    # if covMatList doesn't have names, assign them
    if(is.null(names(covMatList))){
        names(covMatList) <- paste("A",1:m,sep="")
    }
    
    # check for starting values
    if(!is.null(start)){
        if(family$family == "gaussian"){
            if(length(start) != (m+g)){
                stop("length of start must equal the length of covMatList + number of groups for gaussian")
            }
        }else{
            if(length(start) != m){
                stop("length of start must equal the length of covMatList")
            }
        }
    }
    
    # subest covariance structure matrices
    for(i in 1:m){
        if(!all(scan.include %in% colnames(covMatList[[i]]))){
            stop(paste("All of the included Samples must be in matrix", i, "of covMatList"))
        }
        # subset matrix
        keepMat <- colnames(covMatList[[i]]) %in% scan.include
        covMatList[[i]] <- covMatList[[i]][keepMat,keepMat]
        
        # check that names match
        if(!all(colnames(covMatList[[i]]) == scan.include)){
            stop("Column and Row names of matrix", i, "of covMatList must match the scanIDs of scanAnnot")
        }
    }
    
    if(verbose) message("Using AIREML Procedure...")
    if(verbose) message(paste("Sample Size: ", n))
    
    if(family$family == "gaussian"){
        # estimate variance components
        if(verbose) message("Computing Variance Component Estimates...")
        if(verbose) message(paste(paste("Sigma^2_",c(names(covMatList),group.names),sep="", collapse="     "), "log-lik", "RSS", sep="     "))
        
        # estimate variance components
        out <- .runAIREMLgaussian(Y=Y, W=W, start=start, m=m, g=g, n=n, k=k, covMatList=covMatList, group.idx=group.idx, AIREML.tol=AIREML.tol, dropZeros=dropZeros, maxIter=maxIter, verbose=verbose)
        
        mu <- family$linkinv(out$eta)

    }else{
        # initial fit
        fit0 <- glm(Y ~ -1 + W, family=family)
        #beta <- coef(fit0)
        eta <- fit0$linear.predictors  # W %*% beta
        mu <- family$linkinv(eta) # exp(eta)/(1 + exp(eta)) for binomial
        # weights
        vmu <- family$variance(mu) # mu(1-mu) for binomial
        # inverse of g'(mu)
        gmuinv <- family$mu.eta(eta) # = vmu for canonical link
        # working vector
        Y <- eta + (fit0$y - mu)/gmuinv
        
        # set starting values
        newstart <- start
        
        Yreps <- 0
        repeat({
            Yreps <- Yreps + 1
            
            if(verbose) message("Computing Variance Component Estimates...")
            if(verbose) message(paste(paste("Sigma^2_",c(names(covMatList)),sep="", collapse="     "), "log-lik", "RSS", sep="     "))
            
            # estimate variance components
            out <- .runAIREMLother(Y=Y, W=W, start=newstart, m=m, n=n, k=k, covMatList=covMatList, AIREML.tol=AIREML.tol, dropZeros=dropZeros, maxIter=maxIter, verbose=verbose, vmu=vmu, gmuinv=gmuinv)
            
            # update parameters
            if(verbose) message("Updating WorkingY Vector...")
            mu <- family$linkinv(out$eta) # exp(eta)/(1 + exp(eta)) for binomial
            # weights
            vmu <- family$variance(mu) # mu(1-mu) for binomial
            # inverse of g'(mu)
            gmuinv <- family$mu.eta(out$eta) # = vmu for canonical link
            # working vector
            Y <- out$eta + (fit0$y - mu)/gmuinv
            
            # current variance component estimate
            newstart <- out$sigma2.k
            newstart[out$zeroFLAG] <- 2*AIREML.tol
            
            # test for convergence
            stat <- sqrt(sum((out$eta - eta)^2))
            if(verbose) message(paste("Checking for Convergence...", stat, sep = "\t"))
            eta <- out$eta
            if(stat < AIREML.tol){ break() }
            if(Yreps == maxIter){
                out$converged <- FALSE
                warning("Maximum number of iterations for workingY reached without convergence!")
                break()
            }
        })
        
    }
    
    # get estimates and covariance of estimates
    varComp <- out$sigma2.k
    
    if(family$family == "gaussian"){
        names(varComp) <- paste("V_",c(names(covMatList),group.names),sep="")
        varCompCov <- matrix(NA, nrow=(m+g), ncol=(m+g))
        colnames(varCompCov) <- paste("V_",c(names(covMatList),group.names),sep="")
        rownames(varCompCov) <- paste("V_",c(names(covMatList),group.names),sep="")
    }else{
        names(varComp) <- paste("V_",c(names(covMatList)),sep="")
        varCompCov <- matrix(NA, nrow=m, ncol=m)
        colnames(varCompCov) <- paste("V_",c(names(covMatList)),sep="")
        rownames(varCompCov) <- paste("V_",c(names(covMatList)),sep="")
    }
    
    if(dropZeros){
        varCompCov[!out$zeroFLAG, !out$zeroFLAG] <- solve(out$AI)
    }else{
        varCompCov <- solve(out$AI)
    }
    
    # compute Cholesky decomposition of Covariance matrix
    if(verbose) message("Computing Cholesky Decomposition of Inverse Covariance Matrix of Phenotype...")
    # Covariance Matrix
    if(family$family == "gaussian"){
        Sigma <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
        if(g == 1){
            diag(Sigma) <- diag(Sigma) + varComp[m+1]
        }else{
            diagV <- rep(0,n)
            for(i in 1:g){
                diagV[group.idx[[i]]] <- varComp[m+i]
            }
            diag(Sigma) <- diag(Sigma) + diagV
        }
    }else{
        Sigma <- matrix(0, nrow=n, ncol=n)
        for(i in 1:m){
            Sigma <- Sigma + covMatList[[i]]*varComp[i]
        }
        Sigma <- Sigma + diag(as.vector(vmu)/as.vector(gmuinv)^2)
    }
    
    # Inverse
    SigmaInv <- chol2inv(chol(Sigma))
    # Cholesky Decomposition
    cholSigmaInv <- t(chol(SigmaInv))
    colnames(cholSigmaInv) <- colnames(covMatList[[1]])
    rownames(cholSigmaInv) <- rownames(covMatList[[1]])
    
    return(list(varComp = varComp, varCompCov = varCompCov, cholSigmaInv = cholSigmaInv, beta = as.vector(out$beta), workingY=as.vector(Y), eta=as.vector(out$eta), mu=as.vector(mu), converged = out$converged, zeroFLAG = out$zeroFLAG, hetResid = hetResid, logLikR = out$logLikR, logLik = out$logLik, dispersion = out$RSS,Sigma=Sigma,Sinv=SigmaInv))
    
    
}



.runAIREMLgaussian <- function(Y, W, start, m, g, n, k, covMatList, group.idx, AIREML.tol, dropZeros, maxIter, verbose){
    
    # initial values
    sigma2.p <- var(Y)
    AIREML.tol <- AIREML.tol*sigma2.p  # set convergence tolerance dependent on trait
    if(is.null(start)){
        sigma2.k <- rep((1/(m+1))*sigma2.p, (m+g))
    }else{
        sigma2.k <- as.vector(start)
    }
    
    reps <- 0
    repeat({
        reps <- reps+1
        
        zeroFLAG <- sigma2.k < AIREML.tol # which elements have converged to "0"
        sigma2.k[zeroFLAG] <- 0 # set these to 0
        
        # phenotype covariance matrix
        Vre <- Reduce("+", mapply("*", covMatList, sigma2.k[1:m], SIMPLIFY=FALSE))
        
        V <- Vre
        if(g == 1){
            diag(V) <- diag(V) + sigma2.k[m+1]
        }else{
            diagV <- rep(0,n)
            for(i in 1:g){
                diagV[group.idx[[i]]] <- sigma2.k[m+i]
            }
            diag(V) <- diag(V) + diagV
        }
        
        # cholesky decomposition
        cholV <- chol(V)
        # inverse
        Vinv <- chol2inv(cholV)
        VinvW <- crossprod(Vinv,W)
        cholWtVinvW <- chol(crossprod(W, VinvW))
        WtVinvWInv <- chol2inv(cholWtVinvW)
        beta <- crossprod(WtVinvWInv, crossprod(VinvW,Y))
        fits <- tcrossprod(W, t(beta))
        residM <- as.vector(Y - fits)
        VinvR <- crossprod(Vinv, residM)
        RVinvR <- crossprod(residM, VinvR)
        # residual sum of squares
        RSS <- as.numeric(RVinvR/(n-k))
        # log likelihood
        logLik <- as.numeric( -0.5*n*log(2*pi*RSS) - sum(log(diag(cholV))) - 0.5*RVinvR/RSS )
        logLikR <- as.numeric( logLik + 0.5*k*log(2*pi*RSS) - sum(log(diag(cholWtVinvW))) )
        
        # print current estimates
        if(verbose) print(c(sigma2.k, logLikR, RSS))
        
        # projection matrix
        P <- Vinv - tcrossprod(tcrossprod(VinvW,WtVinvWInv),VinvW)
        
        # vector for later use
        PY <- crossprod(P,Y)
        
        if(reps > 1){
            # Average Information and Scores
            AI <- matrix(NA, nrow=(m+g), ncol=(m+g))
            score <- rep(NA,(m+g))
            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                score[i] <- -0.5*(sum(P*covMatList[[i]]) - crossprod(Y, PAPY)) # tr(PA) - YPAPY
                AI[i,i] <- 0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY)) # YPAPAPY
                if((i+1) <= m){
                    for(j in (i+1):m){
                        AI[i,j] <- 0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY)) # YPDPAPY
                        AI[j,i] <- AI[i,j]
                    }
                }
                if(g == 1){
                    AI[i,(m+1)] <- 0.5*crossprod(PY, PAPY) # YPIPAPY
                    AI[(m+1),i] <- AI[i,(m+1)]
                }else{
                    for(j in 1:g){
                        AI[i,m+j] <- 0.5*crossprod(PY[group.idx[[j]]], PAPY[group.idx[[j]]]) # YP(I_group)PAPY
                        AI[m+j,i] <- AI[i,m+j]
                    }
                }
            }
            if(g == 1){
                score[m+1] <- -0.5*(sum(diag(P)) - crossprod(PY)) # tr(P) - YPIPY
                AI[(m+1),(m+1)] <- 0.5*crossprod(PY,crossprod(P,PY)) # YPIPIPY
            }else{
                for(i in 1:g){
                    PIPY <- crossprod(P[group.idx[[i]], ],PY[group.idx[[i]]])
                    score[m+i] <- -0.5*(sum(diag(P)[group.idx[[i]]]) - crossprod(PY[group.idx[[i]]])) # tr(P(I_group)) - YP(I_group)PY
                    AI[m+i,m+i] <- 0.5*crossprod(PY[group.idx[[i]]], PIPY[group.idx[[i]]]) # YP(I_group)P(I_group)PY
                    if((i+1) <= g){
                        for(j in (i+1):g){
                            AI[m+i,m+j] <- 0.5*crossprod(PY[group.idx[[j]]], PIPY[group.idx[[j]]]) # YP(I_group2)P(I_group)PY
                            AI[m+j,m+i] <- AI[m+i,m+j]
                        }
                    }
                }
            }
            
            if(dropZeros){
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            
            if(dropZeros){
                sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + AIinvScore
                sigma2.kplus1[zeroFLAG] <- 0
            }else{
                sigma2.kplus1 <- sigma2.k + AIinvScore
                sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
            }
            
            # step-halving if step too far
            tau <- 1
            while(!all(sigma2.kplus1 >= 0)){
                tau <- 0.5*tau
                if(dropZeros){
                    sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG] <- 0
                }else{
                    sigma2.kplus1 <- sigma2.k + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
                }
            }
            
            # test for convergence
            stat <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            # update estimates
            sigma2.k <- sigma2.kplus1
            if(stat < AIREML.tol){
                converged <- TRUE
                break()
            }
            if(reps == maxIter){
                converged <- FALSE
                warning("Maximum number of iterations reached without convergence!")
                break()
            }
            
        }else{
            # EM step
            sigma2.kplus1 <- rep(NA,(m+g))
            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                sigma2.kplus1[i] <- (1/n)*(sigma2.k[i]^2*crossprod(Y,PAPY) + n*sigma2.k[i] - sigma2.k[i]^2*sum(P*covMatList[[i]]))
            }
            if(g == 1){
                sigma2.kplus1[m+1] <- (1/n)*(sigma2.k[m+1]^2*crossprod(PY) + n*sigma2.k[m+1] - sigma2.k[m+1]^2*sum(diag(P)))
            }else{
                for(i in 1:g){
                    sigma2.kplus1[m+i] <- (1/n)*(sigma2.k[m+i]^2*crossprod(PY[group.idx[[i]]]) + n*sigma2.k[m+i] - sigma2.k[m+i]^2*sum(diag(P)[group.idx[[i]]]))
                }
            }
            sigma2.k <- sigma2.kplus1
        }
        
    })
    
    # linear predictor
    eta <- fits + crossprod(Vre, VinvR) # X\beta + Zb
    
    return(list(sigma2.k = sigma2.k, AI = AI, converged = converged, zeroFLAG = zeroFLAG, beta=beta, eta=eta, logLikR=logLikR, logLik=logLik, RSS=RSS))
    
}



.runAIREMLother <- function(Y, W, start, m, n, k, covMatList, AIREML.tol, dropZeros, maxIter, verbose, vmu, gmuinv){
    
    # initial values for variance components
    if(is.null(start)){
        sigma2.k <- rep(10*sqrt(AIREML.tol), m)
    }else{
        sigma2.k <- as.vector(start)
    }
    
    reps <- 0
    repeat({
        reps <- reps+1
        
        zeroFLAG <- sigma2.k < AIREML.tol # which elements have converged to "0"
        sigma2.k[zeroFLAG] <- 0 # set these to 0
        
        # variance matrix
        Vre <- Reduce("+", mapply("*", covMatList, sigma2.k[1:m], SIMPLIFY=FALSE))
        V <- Vre + diag(as.vector(vmu)/as.vector(gmuinv)^2)

        # cholesky decomposition
        cholV <- chol(V)
        # inverse
        Vinv <- chol2inv(cholV)
        VinvW <- crossprod(Vinv,W)
        cholWtVinvW <- chol(crossprod(W, VinvW))
        WtVinvWInv <- chol2inv(cholWtVinvW)
        beta <- crossprod(WtVinvWInv, crossprod(VinvW,Y))
        fits <- tcrossprod(W, t(beta))
        residM <- as.vector(Y - fits)
        VinvR <- crossprod(Vinv, residM)
        RVinvR <- crossprod(residM, VinvR)
        # residual sum of squares
        RSS <- as.numeric(RVinvR/(n-k))
        # log likelihood
        logLik <- as.numeric( -0.5*n*log(2*pi*RSS) - sum(log(diag(cholV))) - 0.5*RVinvR/RSS )
        logLikR <- as.numeric( logLik + 0.5*k*log(2*pi*RSS) - sum(log(diag(cholWtVinvW))) )
 
        # print current estimates
        if(verbose) print(c(sigma2.k, logLikR, RSS))
        
        # projection matrix
        P <- Vinv - tcrossprod(tcrossprod(VinvW,WtVinvWInv),VinvW)
        
        # matrices for later use
        PY <- crossprod(P,Y)
        
        if(reps > 1){
            # Average Information and Scores
            AI <- matrix(NA, nrow=m, ncol=m)
            score <- rep(NA,m)
            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                score[i] <- -0.5*(sum(P*covMatList[[i]]) - crossprod(Y, PAPY))
                AI[i,i] <- 0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY)) # YPAPAPY
                if((i+1) <= m){
                    for(j in (i+1):m){
                        AI[i,j] <- 0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY)) # YPDPAPY
                        AI[j,i] <- AI[i,j]
                    }
                }
            }
            
            if(dropZeros){
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            
            if(dropZeros){
                sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + AIinvScore
                sigma2.kplus1[zeroFLAG] <- 0
            }else{
                sigma2.kplus1 <- sigma2.k + AIinvScore
                sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
            }
            
            # step-halving if step too far
            tau <- 1
            while(!all(sigma2.kplus1 >= 0)){
                tau <- 0.5*tau
                if(dropZeros){
                    sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG] <- 0
                }else{
                    sigma2.kplus1 <- sigma2.k + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
                }
            }
            
            # test for convergence
            stat <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            sigma2.k <- sigma2.kplus1
            if(stat < AIREML.tol){
                converged <- TRUE
                break()
            }
            if(reps == maxIter){
                converged <- FALSE
                warning("Maximum number of iterations reached without convergence!")
                break()
            }
            
        }else{
            # EM step
            sigma2.kplus1 <- rep(NA,m)
            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                sigma2.kplus1[i] <- (1/n)*((sigma2.k[i])^2*crossprod(Y,PAPY) + (n*sigma2.k[i] - (sigma2.k[i])^2*sum(P*covMatList[[i]])))
            }
            sigma2.k <- sigma2.kplus1
        }
        
    })
    
    # linear predictor
    eta <- fits + crossprod(Vre, VinvR) # X\beta + Zb
    
    return(list(sigma2.k = sigma2.k, AI = AI, converged = converged, zeroFLAG = zeroFLAG, beta = beta, eta = eta, logLikR=logLikR, logLik=logLik, RSS=RSS))

}



# x is the output from estVarComp
estVarCompCI <- function(x, prop=TRUE){
    if(prop){
        if(x$hetResid){ 
            stop("Estimates of proportional variance are not supported with heterogeneous group residual variances")
        }
        ci <- matrix(NA, nrow=length(x$varComp), ncol=2)
        est <- x$varComp/sum(x$varComp)
        varCompCov <- x$varCompCov
        varCompCov[is.na(varCompCov)] <- 0
        for(i in 1:length(est)){
            deltaH <- rep(-x$varComp[i]/(sum(x$varComp)^2),length(x$varComp))
            deltaH[i] <- deltaH[i] + sum(x$varComp)/(sum(x$varComp)^2)
            varH <- crossprod(deltaH, crossprod(varCompCov, deltaH))
            ci[i,] <- est[i] + sqrt(varH)*qnorm(c(0.025,0.975))
        }
        ci[x$zeroFLAG,] <- NA
        res <- as.data.frame(cbind(est, ci))
        names(res) <- c("Proportion", "Lower 95", "Upper 95")
        
    }else{
        ci <- x$varComp + sqrt(diag(x$varCompCov)) %o% qnorm(c(0.025,0.975))
        res <- as.data.frame(cbind(x$varComp, ci))
        names(res) <- c("Est", "Lower 95", "Upper 95")
    }
    
    print(res)
}


.createDesignMatrix <- function(scanAnnot, outcome, covar.vec, ivar.vec, scan.include){
    # set which samples to keep
    scanID <- getScanID(scanAnnot)
      # samples to include for the analysis
      if(!is.null(scan.include)){
            if(!all(scan.include %in% scanID)){
                    stop("Not all of the scan IDs in scan.include are in scanAnnot")
                  }
                keep <- scanID %in% scan.include
          }else{
                keep <- rep(TRUE, length(scanID))
              }

      if(!is.null(covar.vec)){
            cvnames <- unique(unlist(strsplit(covar.vec,"[*:]")))
                # read in data
                dat <- as.data.frame(getVariable(scanAnnot, c(outcome,cvnames)))
                # identify samples with any missing data
                keep <- keep & apply(dat,1,function(x){ all(!is.na(x)) })

                # read in interaction variable data
                if(!is.null(ivar.vec)){
                        ivnames <- unique(unlist(strsplit(ivar.vec,"[*:]")))
                              # read in data
                              idat <- as.data.frame(getVariable(scanAnnot, c(outcome,ivnames)))
                              # identify samples with any missing data
                              keep <- keep & apply(idat,1,function(x){ all(!is.na(x)) })
                              # remove samples with any missing data
                              idat <- as.data.frame(idat[keep,])
                              # create design matrix for interaction variables
                              iformula <- as.formula(paste(paste(outcome,"~"), paste(ivar.vec,collapse="+")))
                              V <- model.matrix(iformula, data=idat)
                              rownames(V) <- scanID[keep]
                      }

                # remove samples with any missing data
                dat <- as.data.frame(dat[keep,])
                # outcome vector
                Y <- dat[,outcome]
                # create design matrix
                model.formula <- as.formula(paste(paste(outcome,"~"), paste(covar.vec,collapse="+")))
                W <- model.matrix(model.formula, data=dat)
                rownames(W) <- scanID[keep]

          }else{
                # read in data
                dat <- getVariable(scanAnnot,outcome)
                    # identify samples with any missing data
                    keep <- keep & !is.na(dat)

                    # read in interaction variable data
                    if(!is.null(ivar.vec)){
                            ivnames <- unique(unlist(strsplit(ivar.vec,"[*:]")))
                                  # read in data
                                  idat <- as.data.frame(getVariable(scanAnnot, c(outcome,ivnames)))
                                  # identify samples with any missing data
                                  keep <- keep & apply(idat,1,function(x){ all(!is.na(x)) })
                                  # remove samples with any missing data
                                  idat <- as.data.frame(idat[keep,])
                                  # create matrix to store interaction variables
                                  iformula <- as.formula(paste(paste(outcome,"~"), paste(ivar.vec,collapse="+")))
                                  V <- model.matrix(iformula, data=idat)
                                  rownames(V) <- scanID[keep]
                          }

                    # outcome
                    Y <- dat[keep]
                    # design matrix
                    W <- matrix(1,nrow=length(Y),ncol=1)
                    rownames(W) <- scanID[keep]
              }

      # return output
      if(!is.null(ivar.vec)){
            return(list(Y = Y, W = W, V = V))
          }else{
                return(list(Y = Y, W = W))
              }

  }
