#secretly doing bantu phylogeny

glottocodes <- convertnames$glottocode[which(convertnames$taxon %in% bantu_phylo$tip.label)]



bantudata <- ac_nona[which(ac_nona$Glottolog.language.dialect.id %in% glottocodes),]
bantudata$treelabel[49] <- "S11_Shona"
bantudata$treelabel[17] <- "H16c_Yombe"

#add labels as in tree
dataglotto <- bantudata$Glottolog.language.dialect.id
corresponding <- convertnames$taxon[match(dataglotto, convertnames$glottocode)]
bantudata$treelabel <- corresponding

missing <- setdiff(glottocodes, ac_nona$Glottolog.language.dialect.id)
missing2 <- convertnames$taxon[which(convertnames$glottocode %in% missing)]

bantu_phylo_pruned <- drop.tip(bantu_phylo, missing2)
finaltree <- nnls.tree(cophenetic(bantu_phylo_pruned),bantu_phylo_pruned,rooted=TRUE)

prior.bantu<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

model0<-MCMCglmm(trance_binary~Variable..Precipitation.Constancy...,
                 random=~treelabel, 
                 ginverse=list(treelabel=inverseA(finaltree)$Ainv), 
                 prior = prior.bantu, 
                 verbose=TRUE, #this will let us see how fast the model is running
                 family="categorical",
                 data = bantudata,
                 nitt=55000, #number of iterations (this is a very small number! a real analysis would run for much longer)
                 thin=50, #sampling of iterations (will record every 50th iteration, for a total posterior sample size of 1000)
                 burnin=5000) #initial iterations to discard

summary(model0)

model1<-MCMCglmm(trance_binary~Variable..Precipitation.Constancy...,
                 random=~treelabel, 
                 ginverse=list(treelabel=inverseA(finaltree)$Ainv), 
                 prior = prior.bantu, 
                 verbose=F, #this will let us see how fast the model is running
                 family="categorical",
                 data = bantudata,
                 nitt=10000000, #number of iterations (this is a very small number! a real analysis would run for much longer)
                 thin=1000, #sampling of iterations (will record every 50th iteration, for a total posterior sample size of 1000)
                 burnin=10000) #initial iterations to discard

summary(model1)
plot(model1)

prior.bantu2<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,fix=1))

model2<-MCMCglmm(trance_binary~Variable..Precipitation.Constancy...,
                 random=~treelabel, 
                 ginverse=list(treelabel=inverseA(finaltree)$Ainv), 
                 prior = prior.bantu2, 
                 verbose=F, #this will let us see how fast the model is running
                 family="categorical",
                 data = bantudata,
                 nitt=1000000, #number of iterations (this is a very small number! a real analysis would run for much longer)
                 thin=500, #sampling of iterations (will record every 50th iteration, for a total posterior sample size of 1000)
                 burnin=100000) #initial iterations to discard

summary(model2)
plot(model2)

lambda<-model2$VCV[,"treelabel"]/(model2$VCV[,"treelabel"]+model2$VCV[,"units"])

mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)

prior.bantu3<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,fix=1))

model3<-MCMCglmm(trance_binary~Variable..Precipitation.Constancy...,
                 random=~treelabel, 
                 ginverse=list(treelabel=inverseA(finaltree)$Ainv), 
                 prior = prior.bantu3, 
                 verbose=F, #this will let us see how fast the model is running
                 family="categorical",
                 data = bantudata,
                 nitt=1000000, #number of iterations (this is a very small number! a real analysis would run for much longer)
                 thin=500, #sampling of iterations (will record every 50th iteration, for a total posterior sample size of 1000)
                 burnin=100000) #initial iterations to discard

summary(model3)
plot(model3)

lambda<-model3$VCV[,"treelabel"]/(model3$VCV[,"treelabel"]+model3$VCV[,"units"])

mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)



install.packages("ggtree")
library(ggtree)



prior.bantuChi<-list(G=list(G1=list(V=1,nu=1000, alpha.mu=0, alpha.V=1)),R=list(V=1,fix=1))

modelChi<-MCMCglmm(trance_binary~Variable..Precipitation.Constancy...,
                 random=~treelabel, 
                 ginverse=list(treelabel=inverseA(finaltree)$Ainv), 
                 prior = prior.bantuChi, 
                 verbose=F, #this will let us see how fast the model is running
                 family="categorical",
                 data = bantudata,
                 nitt=10000000, #number of iterations (this is a very small number! a real analysis would run for much longer)
                 thin=1000, #sampling of iterations (will record every 50th iteration, for a total posterior sample size of 1000)
                 burnin=1000000) #initial iterations to discard

summary(modelChi)
plot(modelChi)

lambda<-modelChi$VCV[,"treelabel"]/(modelChi$VCV[,"treelabel"]+modelChi$VCV[,"units"])

mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)


#plot traits on tree
plot(finaltree,show.tip.label=FALSE)
tiplabels(pie=bantudata$trance_binary, cex=0.4)
#add.simmap.legend(colors=cols) #click where you want to draw your legend!
