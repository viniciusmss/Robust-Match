# Retrieve source code
source("utils.R")

##### DW Lalonde ######

genout.dw <- GenMatch(Tr=lalonde$treat, X=lalonde[,-which(names(lalonde) == "re78")], 
                      BalanceMatrix = lalonde, pop.size=100,
                      print=1, ties=TRUE, wait.generations = 5, 
                      fit.func = robust.fitfunc.plus)

mout.dw <- Match(Y=lalonde$re78, Tr=lalonde$treat, 
                 X=lalonde[,-which(names(lalonde) == "re78")], 
                 ties=TRUE, Weight.matrix=genout.dw)
summary(mout.dw)

MatchBalance(treat ~ age + education  + black + hispanic + married + nodegree + re74 +re75
             data=lalonde, match.out = mout.dw, print.level=1)


##### CPS-3 Lalonde ######

genout.csp3 <- GenMatch(Tr=lalonde.cps3$treat, X=lalonde.cps3[,-which(names(lalonde.cps3) == "re78")], 
                      BalanceMatrix = lalonde.cps3, pop.size=100,
                      print=1, ties=TRUE, wait.generations = 5, 
                      fit.func = robust.fitfunc.plus)

mout.csp3 <- Match(Y=lalonde.cps3$re78, Tr=lalonde.cps3$treat, 
                 X=lalonde.cps3[,-which(names(lalonde.cps3) == "re78")], 
                 ties=TRUE, Weight.matrix=genout.csp3)
summary(mout.csp3)

MatchBalance(treat ~ age + education  + black + hispanic + married + nodegree + re74 +re75
             data=lalonde.cps3, match.out = mouo.csp3, print.level=1)


##### CPS-2 Lalonde ######

genout.cps2 <- GenMatch(Tr=lalonde.cps2$treat, X=lalonde.cps2[,-which(names(lalonde.cps2) == "re78")], 
                        BalanceMatrix = lalonde.cps2, pop.size=100,
                        print=1, ties=TRUE, wait.generations = 5, 
                        fit.func = robust.fitfunc.plus)

mout.cps2 <- Match(Y=lalonde.cps2$re78, Tr=lalonde.cps2$treat, 
                   X=lalonde.cps2[,-which(names(lalonde.cps2) == "re78")], 
                   ties=TRUE, Weight.matrix=genout.cps2)
summary(mout.cps2)

MatchBalance(treat ~ age + education  + black + hispanic + married + nodegree + re74 +re75
             data=lalonde.cps2, match.out = mout.cps2, print.level=1)


##### CPS-1 Lalonde ######

genout.cps1 <- GenMatch(Tr=lalonde.cps1$treat, X=lalonde.cps1[,-which(names(lalonde.cps1) == "re78")], 
                        BalanceMatrix = lalonde.cps1, pop.size=100,
                        print=1, ties=TRUE, wait.generations = 5, 
                        fit.func = robust.fitfunc.plus)

mout.cps1 <- Match(Y=lalonde.cps1$re78, Tr=lalonde.cps1$treat, 
                   X=lalonde.cps1[,-which(names(lalonde.cps1) == "re78")], 
                   ties=TRUE, Weight.matrix=genout.cps1)
summary(mout.cps1)

MatchBalance(treat ~ age + education  + black + hispanic + married + nodegree + re74 +re75
             data=lalonde.cps1, match.out = mout.cps1, print.level=1)
