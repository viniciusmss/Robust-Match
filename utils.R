# Load imports
library(rbounds)
library(rgenoud)
library(Matching)
library(sensitivitymv)
library(foreign)
library(segmented)
library(MatchingFrontier)

# Load Data
lalonde <- read.dta("http://www.nber.org/~rdehejia/data/nsw_dw.dta")[,-1]
dw.treat <- lalonde[lalonde$treat == 1,] 
cps1 <- read.dta("http://www.nber.org/~rdehejia/data/cps_controls.dta")[,-1]
lalonde.cps1 <- rbind(dw.treat, cps1)
cps2 <- read.dta("http://www.nber.org/~rdehejia/data/cps_controls2.dta")[,-1]
lalonde.cps2 <- rbind(dw.treat, cps2)
cps3 <- read.dta("http://www.nber.org/~rdehejia/data/cps_controls3.dta")[,-1]
lalonde.cps3 <- rbind(dw.treat, cps3)

# Create matrix of matched set outcomes
make_ymat<-function(matches, y){
  
  # Get indices of treated and control units
  treated <- unique(matches[, 1]) # Remove repeated indices
  controls <- matches[, 2]        # Keep repeated indices
  
  # Get outcomes of control units
  y_ctrls <- y[controls]
  
  # Create a table to check how many matches for each control  
  trt_table <- table(treated)
  
  # Get number of sets and size of largest set
  n_sets <- length(trt_table)
  max_ctrls <- max(trt_table)
  
  # Smallest table necessary is number of sets vs. size of largest set
  y_mat <- matrix(NA, n_sets, max_ctrls + 1)
  
  m <- 0 # Auxiliary indexer that will run linearly through the matches sets
  
  # For each set
  for (i in 1:n_sets){
    
    y_mat[i, 1] <- y[treated[i]] # Get treated outcome
    
    # Get controls outcomes
    y_mat[i, 2:(1+trt_table[i])] <- y_ctrls[(m+1):(m+trt_table[i])]
    
    # Advance the indexer
    m <- m + trt_table[i]
  }
  y_mat
}

# Sensitivity optim function
robust.fitfunc <- function(matches, BM, gamma=2) {
  
  # robust.fitfunc requires that the last column of the BM matrix contain the outcomes
  # hence, it is VERY important not to run standard pval balance genetic matching
  # with the same BM!
  
  # Get outcomes
  outcomes <- BM[ , ncol(BM)]
  
  # Construct matrix of outcomes
  ymat <- make_ymat(matches, outcomes)
  
  # Compute sensitivity
  pval <- senmv(y=ymat, gamma=gamma, inner=0, trim=Inf, lambda = 1/2, TonT=TRUE)$pval
  
  return(pval)
}

# Get cut point
getCutpoint_ <- function(dataset, base.form, cov, median = TRUE){
  if(!median){
    mod.form <- as.formula(paste(as.character(base.form[2]),
                                 as.character(base.form[1]),
                                 cov))
    
    if(length(unique(dataset[[as.character(base.form[2])]])) == 2){
      base.mod <- glm(mod.form, data = dataset, family = 'binomial')
    } else{
      base.mod <- lm(mod.form, data = dataset)
    }
    
    ### CHANGE : When the median is zero, this leads segmented() to break
    psi0 <-median(dataset[[cov]])
    if (psi0 == 0) psi0 <- mean(dataset[[cov]])
    ###
    
    seg.reg <- segmented(base.mod, seg.Z=mod.form[c(1,3)], 
                         psi = psi0, control = seg.control(it.max = 10000))
    cutpoint <- seg.reg$psi[2]
    
  }
  else{
    cutpoint <- median(dataset[[cov]])
    
    ### CHANGE : When the median is zero, this leads the lm.fit() in 
    # modelDependence() to break
    ###
    if (cutpoint == 0) cutpoint <- 0.001
  }
  return(cutpoint)
}

trim <- function(x){
  gsub("^\\s+|\\s+$", "", x)
}

# Get model dependence
modelDependence_ <- 
  function(dataset, treatment, base.form, verbose = TRUE, seed = 12345, cutpoints = NA, median = TRUE){
    set.seed(seed)
    
    base.form <- as.formula(base.form)
    
    covs <- strsplit(as.character(base.form[3]), '\\+')
    covs <- unlist(lapply(covs, trim))
    
    base.theta <- lm(base.form, data = dataset)$coefficients[[treatment]]
    
    if(verbose){
      cat(paste('Estimate from base model:', round(base.theta, 2), '\n'))
    }
    
    N <- nrow(dataset)
    # estimate theta_p
    
    theta.Ps <- c()
    
    for(cov in covs){
      if(cov == treatment){next}
      
      # Formula for this iteration
      this.form <- paste(as.character(base.form[2]),
                         as.character(base.form[1]),
                         paste(covs[!(covs %in% cov)], collapse = ' + '))
      
      
      base.mod <- lm(base.form, data = dataset)
      
      # Split data
      if(length(unique(dataset[[cov]])) == 2){
        split.inds <- dataset[[cov]] == unique(dataset[[cov]])[1]            
        dat1 <- dataset[split.inds,]
        dat2 <- dataset[!split.inds,]
      }else{
        if(cov %in% names(cutpoints)){
          cutpoint <- cutpoints[names(cutpoints) == cov]
        }else{
          cutpoint <- getCutpoint_(dataset, base.form, cov, median)
        }
        split.inds <- dataset[[cov]] < cutpoint
        dat1 <- dataset[split.inds,]
        dat2 <- dataset[!split.inds,]
      }
      
      # Get theta_ps
      dat1.est <- lm(this.form, data = dat1)$coefficients[[treatment]]
      dat2.est <- lm(this.form, data = dat2)$coefficients[[treatment]]
      
      this.theta.p <- dat1.est * (nrow(dat1) / N) + dat2.est * (nrow(dat2) / N)        
      
      if(verbose){
        cat(paste('Estimate from', cov, 'partition:', round(this.theta.p, 2), '\n'))
      }
      theta.Ps <- c(theta.Ps, this.theta.p)      
    }
    
    covs <- covs[!(covs %in% treatment)]
    failed.covs <-covs[is.na(theta.Ps)]
    
    theta.Ps <- theta.Ps[!is.na(theta.Ps)]
    
    sigma.hat.theta <- sqrt(sum((theta.Ps - base.theta) ^ 2) / length(theta.Ps))
    
    return(sigma.hat.theta)
  }

robust.fitfunc.plus <- function(matches, BM, gamma=2) {
  ### Requirements
  # treatment column should be named 'treat'
  # Balance matrix should have the last column as the outcome
  # Balance matrix should also contain the treatment vector
  # Balance matrix should have appropriate column names
  
  # Get auxiliary objects
  outcomes <- BM[ , ncol(BM)]
  covars <- BM[ , -ncol(BM)]
  vars <- colnames(BM) 
  nvars <- length(vars)
  if (is.null(vars)) stop("Balance Matrix should have appropriate column names.") 
  base.form <- paste(
    vars[nvars],
    "~",
    paste(vars[-nvars], collapse="+")
  )
  
  # Construct matrix of outcomes
  ymat <- make_ymat(matches, outcomes)
  
  # Compute sensitivity
  pval <- senmv(y=ymat, gamma=gamma, inner=0, trim=Inf, lambda = 1/2, TonT=TRUE)$pval
  
  # Compute robustness
  df <- as.data.frame(rbind(BM[matches[,1],], BM[matches[,2],]), stringsAsFactors = FALSE)
  md <- modelDependence_(dataset = df, treatment = 'treat', 
                         verbose=FALSE, base.form = base.form, median=TRUE)
  
  # If robustness p-value falls under significance threshold, 
  # coerce it to zero to focus on robustness to model mispecification
  if (pval < 0.05) pval <- 0
  
  return(c(pval, md))
}

robust.fitfunc.wbal <- function(matches, BM, gamma=2) {
  ### Requirements
  # treatment column should be named 'treat'
  # Balance matrix should have the last column as the outcome
  # Balance matrix should also contain the treatment vector
  # Balance matrix should have appropriate column names
  
  # Get auxiliary objects
  outcomes <- BM[ , ncol(BM)]
  covars <- BM[ , -ncol(BM)]
  vars <- colnames(BM) 
  covar_names <- colnames(covars)
  nvars <- length(vars)
  if (is.null(vars)) stop("Balance Matrix should have appropriate column names.") 
  base.form <- paste(
    vars[nvars],
    "~",
    paste(vars[-nvars], collapse="+")
  )
  
  # Construct matrix of outcomes
  ymat <- make_ymat(matches, outcomes)
  
  # Compute sensitivity
  pval <- senmv(y=ymat, gamma=gamma, inner=0, trim=Inf, lambda = 1/2, TonT=TRUE)$pval
  
  # Compute robustness
  df <- as.data.frame(rbind(BM[matches[,1],], BM[matches[,2],]), stringsAsFactors = FALSE)
  md <- modelDependence_(dataset = df, treatment = 'treat', 
                         verbose=FALSE, base.form = base.form, median=TRUE)
  
  # If robustness p-value falls under significance threshold, 
  # coerce it to zero to focus on robustness to model mispecification
  if (pval < 0.05) pval <- 0
  
  # Compute balance pval
  
  # MatchBalance formula
  bal.form <- as.formula(paste("treat ~",
                               paste(covar_names[covar_names != "treat"], 
                                     collapse="+")))
  # Get the matched dataset
  matcheddf <- covars[c(matches[,1], matches[,2]),]
  mbout <- MatchBalance(bal.form, data=matcheddf, print.level = 0)
  min.pval <- mbout$BMsmallest.p.value
  
  return(c(pval, 1-min.pval, md))
}
