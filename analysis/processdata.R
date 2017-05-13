
setwd("./QMSS2017/groupproject/QMSS2017trance-project")
#read in datafile from dplace for all societies
dplace <- read.csv("./rawData/dplace-societies_masterdata.csv", sep=",", skip=1, stringsAsFactors = F)

#select the atlantic congo groups from the data file
rows <- which(dplace$Language.family == "Atlantic-Congo")
ac <- dplace[rows,]

#remove groups that have missing values for our variables
ac_nona <- ac[!is.na(ac$Code..EA112.Trance.states),]

#remove first column with data source
ac_nona <- ac_nona[,-1]

ac_nona$Original.society.name <- as.character(ac_nona$Original.society.name)
ac_nona$Society.id <- as.character(ac_nona$Society.id)
ac_nona$Preferred.society.name <- as.character(ac_nona$Preferred.society.name)

#write to file
write.csv(ac_nona, file="atlanticCongoprocessed.csv", quote=FALSE)


bantu_phylo <- read.nexus(file="./rawData/bantu.trees-d-place.NEXUS")
ac_glotto <- read.nexus(file="./rawData/Atlantic-Congo.glotto.trees-d-place.NEXUS")


convertnames <- read.csv(file = "./rawData/grollemund_et_al2015/taxa.csv", stringsAsFactors = F)
#convertnames$glottocode[which(tipnames %in% convertnames$taxon)]
#ac_phylo$tip.label <- convertnames$glottocode

#remove society names to keep only society.id
tipnames.short <- sub(".*_", replacement="", x=tipnames)

setdiff(ac_nona$Society.id, tipnames.short)

#select tips that are present in datafile
#setdiff(ac_nona$Glottolog.language.dialect.id, convertnames$glottocode)

#tipstodrop <- setdiff(convertnames$glottocode, ac_nona$Glottolog.language.dialect.id)
#tipstodrop2 <- convertnames$taxon[which(convertnames$glottocode %in% tipstodrop)]

#ac_phylo2 <- drop.tip(ac_phylo, tipstodrop2)



#create binary code for presence/absence of trance, 8 is the current code for absence
ac_nona$trance_binary <- 1
ac_nona$trance_binary[which(ac_nona$Code..EA112.Trance.states == 8)] <- 0
ac_nona$trance_binary[which(ac_nona$Code..EA112.Trance.states == 2)] <- 0


library(MCMCglmm)
library(phangorn)
#make tree ultrametric
ac_glotto.u<-nnls.tree(cophenetic(ac_glotto),ac_glotto,rooted=TRUE)
ac_glotto.u.r <- root(ac_glotto.u, node=115, resolve.root = T)
ac_glotto.u.di <- multi2di(ac_glotto.u)
ac_glotto.u.di.r <- root(ac_glotto.u.di, node=115)


#set priors
prior.ac<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

model0<-MCMCglmm(trance_binary~Variable..Precipitation.Constancy...,
                 random=~treename, 
                 ginverse=list(treename=inverseA(ac_glotto.u.di.r)$Ainv), 
                 prior = prior.PN, 
                 verbose=TRUE, #this will let us see how fast the model is running
                 family="categorical",
                 data = ac_nona,
                 nitt=55000, #number of iterations (this is a very small number! a real analysis would run for much longer)
                 thin=50, #sampling of iterations (will record every 50th iteration, for a total posterior sample size of 1000)
                 burnin=5000) #initial iterations to discard
