# Load imports
library(rbounds)
library(rgenoud)

# Load data
data(lalonde)
attach(lalonde)

# Create data objects
Y <- re78
Tr <- treat
X <- cbind(age, educ, black, hisp, married, nodegr, re74, re75, u74, u75)
BalanceMat <- cbind(age, I(age^2), educ, I(educ^2), black,
                    hisp, married, nodegr, re74 , I(re74^2), re75, I(re75^2),
                    u74, u75, I(re74*re75), I(age*nodegr), I(educ*re74), I(educ*75))

set.seed(12345)
# Genetic matching 
gen1 <- GenMatch(Tr=Tr, X=X, BalanceMat=BalanceMat, pop.size=50,
                 data.type.int=FALSE, print=1, replace=FALSE)

# Matching
mgen1 <- Match(Y=Y, Tr=Tr, X=X, Weight.matrix=gen1, replace=FALSE)
summary(mgen1)

# Balance
mbgen1 <- MatchBalance(Tr ~ age + I(age^2)+ educ + I(educ^2) + black +
                         hisp + married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +
                         u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*75),
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

mbgen2 <- MatchBalance(Tr ~ age + I(age^2)+ educ + I(educ^2) + black +
                         hisp + married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +
                         u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*75),
                       data=lalonde, match.out = mgen2, print=0)

stopifnot(gen1$value[1] == mbgen1$AMsmallest.p.value)

m = sum(gen2$matches[, 2] != mgen2$index.control)
n = nrow(gen2$matches)
print(sprintf("%d / %d are different matches", m, n))

print("yes :)")

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

ymat <- make_ymat(re78, matches)

library(sensitivitymv)
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

MatchBalance(Tr ~ age + I(age^2)+ educ + I(educ^2) + black +
                   hisp + married + nodegr + re74 + I(re74^2) + re75 + I(re75^2) +
                   u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*75),
                   data=lalonde, match.out = mout3, print.level=1)
