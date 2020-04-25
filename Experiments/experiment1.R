# based on oct.genmatch.pscore
# originally from the replication files for Diamond & Sekhon
# lapo:/home/adiamond/midwest/DW/forJas/oct.genmatch.pscore

set.seed(2423045)
rm(list=ls())
source("C:\\Users\\Vinic\\Research\\gensens\\utils.R")


X <-   as.matrix(cbind( lalonde.cps1$age,
                        lalonde.cps1$educ,
                        lalonde.cps1$black,
                        lalonde.cps1$hispan,
                        lalonde.cps1$married,
                        lalonde.cps1$nodegree,
                        I(lalonde.cps1$re74/1000),
                        I(lalonde.cps1$re75/1000)))

BalanceMat <- as.matrix(cbind(
  lalonde.cps1$age,
  lalonde.cps1$educ,
  lalonde.cps1$black,
  lalonde.cps1$hispan,
  lalonde.cps1$married,
  lalonde.cps1$nodegree,
  I(lalonde.cps1$re74/1000),
  I(lalonde.cps1$re75/1000),
  I((lalonde.cps1$re74/1000)^2), # #
  I((lalonde.cps1$re75/1000)^2), # #
  
  I((lalonde.cps1$age)^2),
  I((lalonde.cps1$educ)^2),
  I((lalonde.cps1$black)^2), # #
  I((lalonde.cps1$hispan)^2), # #
  I((lalonde.cps1$married)^2), # #
  I((lalonde.cps1$nodegree)^2), # #
  
  I(lalonde.cps1$age*lalonde.cps1$educ),
  I(lalonde.cps1$age*lalonde.cps1$black),
  I(lalonde.cps1$age*lalonde.cps1$hispan),
  I(lalonde.cps1$age*lalonde.cps1$married),    
  I(lalonde.cps1$age*lalonde.cps1$nodegree),
  I(lalonde.cps1$age*(lalonde.cps1$re74/1000)),
  I(lalonde.cps1$age*(lalonde.cps1$re75/1000)),
  
  I(lalonde.cps1$educ*lalonde.cps1$black),
  I(lalonde.cps1$educ*lalonde.cps1$hispan),
  I(lalonde.cps1$educ*lalonde.cps1$married),    
  I(lalonde.cps1$educ*lalonde.cps1$nodegree), # #
  I(lalonde.cps1$educ*(lalonde.cps1$re74/1000)),
  I(lalonde.cps1$educ*(lalonde.cps1$re75/1000)),
  
  I(lalonde.cps1$black*lalonde.cps1$hispan), # #
  I(lalonde.cps1$black*lalonde.cps1$married),
  I(lalonde.cps1$black*lalonde.cps1$nodegree),
  I(lalonde.cps1$black*(lalonde.cps1$re74/1000)),
  I(lalonde.cps1$black*(lalonde.cps1$re75/1000)),
  
  I(lalonde.cps1$hispan*lalonde.cps1$married),
  I(lalonde.cps1$hispan*lalonde.cps1$nodegree),
  I(lalonde.cps1$hispan*(lalonde.cps1$re74/1000)),
  I(lalonde.cps1$hispan*(lalonde.cps1$re75/1000)),
  
  I(lalonde.cps1$married*lalonde.cps1$nodegree),
  I(lalonde.cps1$married*(lalonde.cps1$re74/1000)),
  I(lalonde.cps1$married*(lalonde.cps1$re75/1000)),
  
  I(lalonde.cps1$nodegree*(lalonde.cps1$re74/1000)), 
  I(lalonde.cps1$nodegree*(lalonde.cps1$re75/1000)),                        
  
  I((lalonde.cps1$re74/1000)*(lalonde.cps1$re75/1000)) # #                                                 
) )


### ORTHOGONALIZE TO PROPENSITY SCORE

# From Review of Econ and Statistics
model.A = lalonde.cps1$treat~I(lalonde.cps1$age) + I(lalonde.cps1$age^2) + I(lalonde.cps1$age^3) + I(lalonde.cps1$educ) + I(lalonde.cps1$educ^2) + I(lalonde.cps1$married) + I(lalonde.cps1$nodegree) + I(lalonde.cps1$black) + I(lalonde.cps1$hispan) + I(lalonde.cps1$re74) + I(lalonde.cps1$re75) + I(lalonde.cps1$re74 == 0) + I(lalonde.cps1$re75 == 0) + I(I(lalonde.cps1$re74)*I(lalonde.cps1$educ))

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


# # #
# # # GENETIC MATCHING
# # #

# pvals <- c()
# effects <- c() 

cl <- makeCluster(detectCores()-2)

best_weights<- c(402.20358, 445.24018,  88.50716, 
                 635.81127, 249.42384, 702.71309,  
                 27.29682, 271.55110, 489.35055)

GM.out <- GenMatch(Tr = lalonde.cps1$treat,
                   X = orthoX2.plus.pscore,
                   BalanceMatrix = BalanceMat,
                   starting.values = best_weights,
                   pop.size = 1,           
                   max.generations = 1,
                   wait.generations = 5,
                   hard.generation.limit = TRUE,
                   print.level = 1)

# Solution Lexical Fitness Value:
#   3.143063e-01  3.173158e-01  3.173158e-01  3.173158e-01  3.173158e-01  3.173158e-01  3.173158e-01  3.173158e-01  3.173158e-01  3.215911e-01  3.270203e-01  3.339515e-01  3.508143e-01  3.570708e-01  3.610729e-01  3.658353e-01  4.373303e-01  4.435542e-01  4.793119e-01  4.864039e-01  5.023163e-01  5.070539e-01  5.324422e-01  5.492617e-01  5.640577e-01  5.773921e-01  6.401642e-01  6.570821e-01  6.605268e-01  6.886729e-01  7.327633e-01  7.391868e-01  7.714773e-01  8.390976e-01  8.390976e-01  8.390976e-01  8.471388e-01  8.801639e-01  9.004493e-01  9.004493e-01  9.004493e-01  9.004493e-01  9.333282e-01  9.471080e-01  9.471278e-01  9.471278e-01  9.631023e-01  9.774258e-01  9.774258e-01  9.774258e-01  9.986865e-01  9.998894e-01  9.998894e-01  9.998894e-01  9.998894e-01  9.998894e-01  9.999977e-01  9.999977e-01  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  
# 
# Parameters at the Solution:
#   
#   X[ 1] :	4.022036e+02
# X[ 2] :	4.452402e+02
# X[ 3] :	8.850716e+01
# X[ 4] :	6.358113e+02
# X[ 5] :	2.494238e+02
# X[ 6] :	7.027131e+02
# X[ 7] :	2.729682e+01
# X[ 8] :	2.715511e+02
# X[ 9] :	4.893505e+02
# 
# Solution Found Generation 1
# Number of Generations Run 1
# 
# Sat Apr 25 09:29:05 2020
# Total run time : 0 hours 0 minutes and 1 seconds

mout <- Match(Y = lalonde.cps1$re78,
              Tr = lalonde.cps1$treat,
              X = orthoX2.plus.pscore,
              Weight.matrix = GM.out)
summary(mout)
# Estimate...  1735.7 
# AI SE......  972.87 
# T-stat.....  1.7841 
# p.val......  0.074404 

mbout <- MatchBalance(treat ~ age+  education+  black+  hispanic+  married+  nodegree+  I(re74/1000)+  I(re75/1000)+  I((re74/1000)^2)+ 
                        I((re75/1000)^2)+  I((age)^2)+  I((education)^2)+  I((black)^2)+  I((hispanic)^2)+ I((married)^2)+
                        I((nodegree)^2)+ I(age*education)+ I(age*black)+ I(age*hispanic)+I(age*married)+ I(age*nodegree)+
                        I(age*(re74/1000))+ I(age*(re75/1000))+  I(education*black)+ I(education*hispanic)+ I(education*married)+ I(education*nodegree)+
                        I(education*(re74/1000))+ I(education*(re75/1000))+ I(black*hispanic)+ I(black*married)+I(black*nodegree)+I(black*(re74/1000))+
                        I(black*(re75/1000))+I(hispanic*married)+I(hispanic*nodegree)+ I(hispanic*(re74/1000))+ I(hispanic*(re75/1000))+ I(married*nodegree)+
                        I(married*(re74/1000))+ I(married*(re75/1000))+  I(nodegree*(re74/1000))+ I(nodegree*(re75/1000))+  I((re74/1000)*(re75/1000)), 
                      data = lalonde.cps1, match.out = mout) 

# After Matching Minimum p.value: 0.21 
# Variable Name(s): I(hispanic * (re74/1000))  Number(s): 37

robust.fitfunc.wbal(GM.out$matches, lalonde.cps1)
# 0.1796675  0.4120000 50.3158348
# Rosenbaum's Gamma = 2 p.val = 0.179
# 1 - balance p.val on raw covar = 0.412 ===> balance p.val = 0.59
# Athey & Imbens model dependence std = 50

y_mat <- make_ymat(GM.out$matches, lalonde.cps1$re78)
for (i in seq(1, 2, 0.1)) {
  p <- robust.fitfunc(gamma=i, GM.out$matches, lalonde.cps1)
  cat("gamma: ", i, "\t p-value:", p, "\n")
}
# gamma:  1.5 	 p-value: 0.0177421 
# gamma:  1.6 	 p-value: 0.03334341 
# gamma:  1.7 	 p-value: 0.05673043 
# gamma:  1.8 	 p-value: 0.08890003 

