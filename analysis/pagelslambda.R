# (2) Testing for phylogenetic signal using Pagel's lambda with the fitDiscrete function in geiger.
library(geiger)
# Pagel's lambda is a tree transformation that assesses the degree of phylogenetic signal within the trait by multiplying the internal branches of the tree by values between 0 and 1 - which is the lambda parameter. A lambda of 1 indicates complete phylogenetic patterning as it returns the branch lengths untransformed while 0 indicates no patterning as it leads the tree collapsing to a single large polytomy. To visualize this we can transform the tree with lambda=0 using the lambdaTree function in geiger and compare it to the original tree

lambda0<-rescale(finaltree, "lambda",0)#transforms the tree topology to one that has all internal branch lengths multiplied by 0 (i.e. lambda=0) creating one giant basal polytomy
par(mfrow=c(1,2))#remember from day1 session2 that this sets the graphical parameters so that the plotting device has 1 row and 2 colums, so we can now plot two trees next to each other.
plot(finaltree)
plot(lambda0)

#Now find the maximum likelihood estimate of lambda for diet 
trance_lambda<-fitDiscrete(finaltree, bantudata$trance_binary, treeTransform="lambda")

# The output returns a list with 4 components. The two of interest are the maximum likelihood estimate of Lambda ($treeParam)  and the negative log likelihood ($lnl). 

# To see if this indicates significant phylogenetic signal we can compare the negative log likelihood when there is no signal i.e. using the tree transformed lambda=0, to that estimated from the original topology.

diet_lambda0<-fitDiscrete(lambda0, mydata$diet)

# You can then compare the negative log likelihood from this analysis to that when lambda was estimated using the original tree topology using a likelihood ratio test (or AIC etc.).

# Likelihood ratio test approximated by a chi-squared distribution
1-pchisq(2*(diet_lambda0$Trait1$lnl-diet_lambda$Trait1$lnl),1)

# (3) Testing whether rates have increased or slowed over evolutionary time with the fitDiscrete function in geiger. If rates decrease over time it is a signature of adaptive radiation if the traits are ecologically relevant.

