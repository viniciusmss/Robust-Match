# Retrieve source code
source("utils.R")
##### CPS-1 Lalonde ######

attach(lalonde.cps1)
X <-   as.matrix(cbind( age,
                        education,
                        black,
                        hispanic,
                        married,
                        nodegree,
                        I(re74/1000),
                        I(re75/1000)))

### ORTHOGONALIZE TO PROPENSITY SCORE

# From Review of Econ and Statistics
model.A = treat~I(age) + I(age^2) + I(age^3) + I(education) + I(education^2) + I(married) + I(nodegree) + I(black) + I(hispanic) + I(re74) + I(re75) + I(re74 == 0) + I(re75 == 0) + I(I(re74)*I(education))

pscores.A <- glm(model.A, family = binomial)

X2 <- X
for (i in 1:ncol(X2)) {
  lm1 = lm(X2[,i]~pscores.A$linear.pred)
  X2[,i] = lm1$residual
}


# VARIABLES WE MATCH ON BEGIN WITH PROPENSITY SCORE
orthoX2.plus.pscore = cbind(pscores.A$linear.pred, X2)

orthoX2.plus.pscore[,1] <- orthoX2.plus.pscore[,1] -
  mean(orthoX2.plus.pscore[,1])


# NORMALIZE ALL COVARS BY STANDARD DEVIATION
for (i in 1:(dim(orthoX2.plus.pscore)[2])) {
  orthoX2.plus.pscore[,i] <-
    orthoX2.plus.pscore[,i]/sqrt(var(orthoX2.plus.pscore[,i]))
}


# This part does not work because the model dependence calculation does not 
# accept the interaction terms

# # CREATE BALANCE MATRIX
# BalanceMat <- as.matrix(cbind(
#   lalonde.cps1$age,
#   lalonde.cps1$educ,
#   lalonde.cps1$black,
#   lalonde.cps1$hispan,
#   lalonde.cps1$married,
#   lalonde.cps1$nodegree,
#   I(lalonde.cps1$re74/1000),
#   I(lalonde.cps1$re75/1000),
#   I((lalonde.cps1$re74/1000)^2), # #
#   I((lalonde.cps1$re75/1000)^2), # #
#   
#   I((lalonde.cps1$age)^2),
#   I((lalonde.cps1$educ)^2),
#   I((lalonde.cps1$black)^2), # #
#   I((lalonde.cps1$hispan)^2), # #
#   I((lalonde.cps1$married)^2), # #
#   I((lalonde.cps1$nodegree)^2), # #
#   
#   I(lalonde.cps1$age*lalonde.cps1$educ),
#   I(lalonde.cps1$age*lalonde.cps1$black),
#   I(lalonde.cps1$age*lalonde.cps1$hispan),
#   I(lalonde.cps1$age*lalonde.cps1$married),    
#   I(lalonde.cps1$age*lalonde.cps1$nodegree),
#   I(lalonde.cps1$age*(lalonde.cps1$re74/1000)),
#   I(lalonde.cps1$age*(lalonde.cps1$re75/1000)),
#   
#   I(lalonde.cps1$educ*lalonde.cps1$black),
#   I(lalonde.cps1$educ*lalonde.cps1$hispan),
#   I(lalonde.cps1$educ*lalonde.cps1$married),    
#   I(lalonde.cps1$educ*lalonde.cps1$nodegree), # #
#   I(lalonde.cps1$educ*(lalonde.cps1$re74/1000)),
#   I(lalonde.cps1$educ*(lalonde.cps1$re75/1000)),
#   
#   I(lalonde.cps1$black*lalonde.cps1$hispan), # #
#   I(lalonde.cps1$black*lalonde.cps1$married),
#   I(lalonde.cps1$black*lalonde.cps1$nodegree),
#   I(lalonde.cps1$black*(lalonde.cps1$re74/1000)),
#   I(lalonde.cps1$black*(lalonde.cps1$re75/1000)),
#   
#   I(lalonde.cps1$hispan*lalonde.cps1$married),
#   I(lalonde.cps1$hispan*lalonde.cps1$nodegree),
#   I(lalonde.cps1$hispan*(lalonde.cps1$re74/1000)),
#   I(lalonde.cps1$hispan*(lalonde.cps1$re75/1000)),
#   
#   I(lalonde.cps1$married*lalonde.cps1$nodegree),
#   I(lalonde.cps1$married*(lalonde.cps1$re74/1000)),
#   I(lalonde.cps1$married*(lalonde.cps1$re75/1000)),
#   
#   I(lalonde.cps1$nodegree*(lalonde.cps1$re74/1000)), 
#   I(lalonde.cps1$nodegree*(lalonde.cps1$re75/1000)),                        
#   
#   I((lalonde.cps1$re74/1000)*(lalonde.cps1$re75/1000)) # #                                                 
# ) )
# 
# # ADD TREATMENT STATUS AND OUTCOME
# BalanceMat <- cbind(BalanceMat, lalonde.cps1$treat, lalonde.cps1$re78)
# 
# # ADD NAMES FOR FORMULA
# colnames(BalanceMat) <- c(1:(ncol(BalanceMat)-2), "treat", "re78")
# colnames(BalanceMat)

detach(lalonde.cps1)

### GENETIC MATCHING

genout.cps1 <- GenMatch(Tr=lalonde.cps1$treat, X=orthoX2.plus.pscore, 
                        BalanceMatrix = lalonde.cps1, pop.size=2000,
                        starting.values = diag(genout.cps1$Weight.matrix),
                        print=1, ties=TRUE, wait.generations = 25, 
                        fit.func = robust.fitfunc.wbal)

mout.cps1 <- Match(Y=lalonde.cps1$re78, Tr=lalonde.cps1$treat, 
                   X=orthoX2.plus.pscore, 
                   ties=TRUE, Weight.matrix=genout.cps1)
summary(mout.cps1)
# Estimate...  2579.7 
# AI SE......  869.98 
# T-stat.....  2.9652 
# p.val......  0.0030247

MatchBalance(treat ~ age + education  + black + hispanic + married + nodegree + re74 +re75,
             data=lalonde.cps1, match.out = mout.cps1, print.level=1)
# After Matching Minimum p.value: 0.11913 
# Variable Name(s): age  Number(s): 1

MatchBalance(treat ~ age + education  + black + hispanic + married + nodegree + re74 +re75,
             data=rbind(lalonde.cps1[mout.cps1$index.treated,],
                        lalonde.cps1[mout.cps1$index.control,]), 
             paired=TRUE, print.level=1)
# Before Matching Minimum p.value: 0.82194 
# Variable Name(s): age  Number(s): 1

robust.fitfunc.wbal(genout.cps1$matches, lalonde.cps1)
# 0.0000000  0.1780606 88.8178924
# Rosenbaum's Gamma = 2 p.val < 0.05 (fun coerces to 0)
# 1 - balance p.val on raw covar = 0.178 ===> balance p.val = 0.82
# Athey & Imbens model dependence std = 88


# What's the sensitivity?
y_mat <- make_ymat(genout.cps1$matches, lalonde.cps1$re78)
for (i in seq(1, 2.2, 0.1)) {
  p <- robust.fitfunc(gamma=i, genout.cps1$matches, lalonde.cps1)
  cat("gamma: ", i, "\t p-value:", p, "\n")
}
# gamma:  1.8 	 p-value: 0.01685837 
# gamma:  1.9 	 p-value: 0.02811927 
# gamma:  2 	   p-value: 0.04389178 
# gamma:  2.1  	 p-value: 0.06474296 
# gamma:  2.2 	 p-value: 0.09097257
