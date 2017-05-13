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

