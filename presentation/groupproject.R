societies <- read.csv("../QMSS2017/groupproject/dplace-societies-2017-05-09.csv", header=T)
soc_ids <- societies$Society.id

atlas <- read.csv("./dplace-data-master/datasets/EA/data.csv")
atlas_sub <- atlas[which(atlas$soc_id %in% soc_ids),]
summary(atlas_sub$code)
