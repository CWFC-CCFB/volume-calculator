# 1 initial form
taper9.nlme <- nlme(dob2~dbh^2*(a1*dbh^a2-h)/(a1*dbh^a2-1.3)*(h/1.3)^(2-2*b),
                    data = taper9, fixed = a1+a2+b~1, random = b~1|idtree, start = c(a1 = 1, a2 = 0.2, b = 1))


# 2a random effects
taper9.nlme24<-update(taper9.nlme,
                      random = list (
                        prov = b~1,
                        idPlot =a2+b~1,
                        idtree = b~1
                      )
)
# 2b random variable varied by province 
taper9.nlme25<-update(taper9.nlme,
                      random = list (
                        prov = b~1,
                        idPlot =a2+b~1,
                        idtree = b~1|prov
                      )
)



# 3 variance function varExp
tapaer9.nlme310<-update(taper9.nlme24, weights = varExp(form = ~dbh)) 
taper9.nlme311<-update(taper9.nlme310,
                       random = list (
                         prov = b~1,
                         idPlot =a2~1,
                         idtree = b~1
                       ))
taper9.nlme321 = update(taper9.nlme311, weights = varExp(form = ~dbh|prov)) 

# 3b variance function varPower
taper9.nlme331<-update(taper9.nlme321, weights = varPower(form = ~fitted(.)|prov)) 


# 4 correlation
taper9.nlme4<-update(taper9.nlme321, correlation = corCAR1(form = ~h))