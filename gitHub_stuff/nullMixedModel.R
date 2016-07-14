nullMixedModel <- function (scanAnnot, VC, outcome, covar.vec = NULL, scan.include = NULL, 
          group.var = NULL, verbose = TRUE) 
{
  cholSigmaInv <- VC$cholSigmaInv
  varComp <- VC$varComp
  if (verbose) 
    message("Reading in Phenotype and Covariate Data...")
  dat <- .createDesignMatrix(scanAnnot = scanAnnot, outcome = outcome, 
                             covar.vec = covar.vec, ivar.vec = NULL, scan.include = scan.include)
  Y <- dat$Y
  W <- dat$W
  k <- ncol(W)
  rm(dat)
  scan.include <- rownames(W)
  if (!all(scan.include %in% colnames(cholSigmaInv))) {
    stop("All of the included Samples must be in the cholSigmaInv matrix")
  }
  chol.idx <- which(!(colnames(cholSigmaInv) %in% scan.include))
  cholSigmaInv <- .subsetCholSigmaInv(cholSigmaInv, chol.idx)
  n <- length(scan.include)
  if (verbose) 
    message("Fitting Model with ", n, " Samples")
  CW <- crossprod(cholSigmaInv, W)
  CY <- crossprod(cholSigmaInv, Y)
  XtSigInvX <- crossprod(CW)
  XtSigInvXInv <- chol2inv(chol(XtSigInvX))
  beta <- crossprod(XtSigInvXInv, crossprod(CW, CY))
  fits <- tcrossprod(W, t(beta))
  residM <- as.vector(Y - fits)
  m <- length(varComp)
  residtmp <- crossprod(residM, cholSigmaInv)
  if (VC$hetResid) {
    if (is.null(group.var)) {
      stop("group.var must be specified when heterogeneous residual variances were used to estimate variance components")
    }
    sigma2epsilon <- rep(NA, n)
    group <- getVariable(scanAnnot, group.var)[getScanID(scanAnnot) %in% 
                                                 scan.include]
    group.names <- as.character(unique(group))
    g <- length(group.names)
    for (i in 1:g) {
      sigma2epsilon[group == group.names[i]] <- varComp[m - 
                                                          g + i]
    }
    residC <- as.vector(sigma2epsilon * tcrossprod(cholSigmaInv, 
                                                   residtmp))
  }
  else {
    residC <- as.vector(varComp[m] * tcrossprod(cholSigmaInv, 
                                                residtmp))
  }
  RSS <- sum(residtmp^2)/(n - k)
  Vbeta <- RSS * XtSigInvXInv
  SE <- sqrt(diag(Vbeta))
  Stat <- (beta/SE)^2
  pval <- pchisq(Stat, df = 1, lower.tail = FALSE)
  logLik <- as.numeric(-0.5 * n * log(2 * pi * RSS) + sum(log(diag(cholSigmaInv))) - 
                         0.5 * tcrossprod(residtmp)/RSS)
  AIC <- 2 * (k + m) - 2 * logLik
  logLikR <- as.numeric(logLik + 0.5 * k * log(2 * pi * RSS) - 
                          0.5 * log(det(XtSigInvX)))
  dimnames(Vbeta) <- list(colnames(W), colnames(W))
  fixef <- data.frame(Est = beta, SE = SE, Stat = Stat, pval = pval)
  rownames(fixef) <- colnames(W)
  return(list(fixef = fixef, varComp = varComp, resid.marginal = residM, 
              resid.conditional = residC, logLik = logLik, AIC = AIC, 
              logLikR = logLikR, model.matrix = W, Vbeta = Vbeta, RSS = RSS))
}
