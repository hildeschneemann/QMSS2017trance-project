
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

#write to file
write.csv(ac_nona, file="atlanticCongoprocessed.csv", quote=FALSE)


ac_phylo <- read.nexus(file="./rawData/grollemund_et_al2015/original/grollemund.mcct.trees")


convertnames <- read.csv(file = "./rawData/grollemund_et_al2015/taxa.csv", stringsAsFactors = F)
#convertnames$glottocode[which(tipnames %in% convertnames$taxon)]
#ac_phylo$tip.label <- convertnames$glottocode

#remove society names to keep only society.id
tipnames.short <- sub(".*_", replacement="", x=tipnames)

setdiff(ac_nona$Society.id, tipnames.short)

#select tips that are present in datafile
#setdiff(ac_nona$Glottolog.language.dialect.id, convertnames$glottocode)

tipstodrop <- setdiff(convertnames$glottocode, ac_nona$Glottolog.language.dialect.id)
tipstodrop2 <- convertnames$taxon[which(convertnames$glottocode %in% tipstodrop)]

ac_phylo2 <- drop.tip(ac_phylo, tipstodrop2)

ac_nona$Original.society.name <- as.character(ac_nona$Original.society.name)
ac_nona$Society.id <- as.character(ac_nona$Society.id)
ac_nona$Preferred.society.name <- as.character(ac_nona$Preferred.society.name)




library(MCMCglmm)

model0<-MCMCglmm(Code..EA112.Trance.states~Variable..Precipitation.Constancy...,
                 random=~treename, 
                 ginverse=list(treename=inverseA(PN.trimmed.u)$Ainv), 
                 prior = prior.PN, 
                 verbose=TRUE, #this will let us see how fast the model is running
                 family="gaussian",
                 data = ac_nona,
                 nitt=55000, #number of iterations (this is a very small number! a real analysis would run for much longer)
                 thin=50, #sampling of iterations (will record every 50th iteration, for a total posterior sample size of 1000)
                 burnin=5000) #initial iterations to discard
