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
cps3 <- read.dta("http://www.nber.org/~rdehejia/data/cps_controls3.dta")[,-1]
lalonde.cps3 <- rbind(dw.treat, cps3)

# Create data objects
attach(lalonde)
Y <- re78
Tr <- treat
X <- cbind(age, education, black, hispanic, married, nodegree, re74, re75)
BalanceMat <- cbind(age, I(age^2), education, I(education^2), black,
                    hispanic, married, nodegree, re74 , I(re74^2), re75, I(re75^2),
                    I(re74*re75), I(age*nodegree), I(education*re74), I(education*re75))

set.seed(12345)
# Genetic matching 
gen1 <- GenMatch(Tr=Tr, X=X, BalanceMat=BalanceMat, pop.size=50,
                 data.type.int=FALSE, print=1, replace=FALSE)

# Matching
mgen1 <- Match(Y=Y, Tr=Tr, X=X, Weight.matrix=gen1, replace=FALSE)
summary(mgen1)

# Balance
mbgen1 <- MatchBalance(Tr ~ age + I(age^2)+ education + I(education^2) + black +
                         hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
                           I(re74*re75) + I(age*nodegree) + I(education*re74) + I(education*re75),
                       data=lalonde, match.out = mgen1, print.level=0)

# Sensitivity Analysis
p <- psens(mgen1, Gamma=1.5, GammaInc=.1)
p
1 - p$bounds$`Upper bound`
hlsens(mgen1, Gamma=1.5, GammaInc=.1)

# Is it because of the stochastic behavior of ties = FALSE?

gen2 <- GenMatch(Tr=Tr, X=X, BalanceMat=BalanceMat, pop.size=50,
                 data.type.int=FALSE, print=0)

mgen2 <- Match(Y=Y, Tr=Tr, X=X, Weight.matrix=gen2)

mbgen2 <- MatchBalance(Tr ~ age + I(age^2)+ education + I(education^2) + black +
                         hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
                        I(re74*re75) + I(age*nodegree) + I(education*re74) + I(education*re75),
                       data=lalonde, match.out = mgen2, print=0)

stopifnot(gen1$value[1] == mbgen1$AMsmallest.p.value)

m = sum(gen2$matches[, 2] != mgen2$index.control)
n = nrow(gen2$matches)
print(sprintf("%d / %d are different matches", m, n))

print("yes :)")

####################################
# WITH ROBUSTNESS TO UNOBSERVABLES #
####################################

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

ymat <- make_ymat(gen2$matches, re78)

senmv.out <- senmv(y=ymat, gamma=1.5, inner=0, trim=Inf, lambda = 1/2, TonT=TRUE)
senmv.out$pval  ## Note that this is a one-tailed test of positive treatment effect


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

# Demonstration of optimization function
BM <- cbind(BalanceMat, re78)
df <- data.frame(trt = mgen2$index.treated, ctrl = mgen2$index.control)
robust.fitfunc(df, BM)           # With matches from Match()
robust.fitfunc(gen2$matches, BM) # With matches from GenMatch()

# Genetic Optimization
genout3 <- GenMatch(Tr=Tr, X=X, BalanceMat=BM, pop.size=100,
                    print=1, ties=TRUE, wait.generations = 5, 
                    fit.func = robust.fitfunc)

mout3 <- Match(Y=Y, Tr=treat, X=X, ties=TRUE, Weight.matrix=genout3)
summary(mout3)

robust.fitfunc(genout3$matches, BM, gamma=2)
robust.fitfunc(cbind(mout3$index.treated, mout3$index.control), BM, gamma=2)

MatchBalance(Tr ~ age + I(age^2)+ education + I(education^2) + black +
                   hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
                     I(re74*re75) + I(age*nodegree) + I(education*re74) + I(education*re75),
                   data=lalonde, match.out = mout3, print.level=1)


######################################
# WITH ROBUSTNESS TO MISPECIFICATION #
######################################

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


base.form <- as.formula('re78 ~ treat + age + education + black + hispanic + married + nodegree + re74 + re75')
md <- modelDependence_(dataset = lalonde.cps3, treatment = 'treat', verbose=TRUE, base.form = base.form, median=TRUE)
print(md)

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
  
  return(c(pval, md))
}

# Genetic Optimization
genout3 <- GenMatch(Tr=lalonde.cps3$treat, X=lalonde.cps3[,-which(names(lalonde.cps3) == "re78")], 
                    BalanceMatrix = lalonde.cps3, pop.size=100,
                    print=1, ties=TRUE, wait.generations = 5, 
                    fit.func = robust.fitfunc.plus)

mout3 <- Match(Y=lalonde.cps3$re78, Tr=lalonde.cps3$treat, X=lalonde.cps3[,-which(names(lalonde.cps3) == "re78")],
               ties=TRUE, Weight.matrix=genout3)
summary(mout3)

robust.fitfunc.plus(genout3$matches, lalonde.cps3, gamma=2)
robust.fitfunc.plus(cbind(mout3$index.treated, mout3$index.control), lalonde.cps3, gamma=2)

MatchBalance(treat ~ age + I(age^2)+ education + I(education^2) + black +
               hispanic + married + nodegree + re74 + I(re74^2) + re75 + I(re75^2) +
               I(re74*re75) + I(age*nodegree) + I(education*re74) + I(education*re75),
             data=lalonde.cps3, match.out = mout3, print.level=1)
