# Plotting phylogenetic tree for group challenge

# You need two datafiles to produce the graph:
# 1. the phylogenetic tree, which here is bantu_phylo_pruned
# 2. the datafile that includes rainfall etc, which here is bantudata

# We first set a color palette (this is from the course notes)
cols<-setNames(palette()[1:length(unique(bantudata$trance_binary))],sort(unique(bantudata$trance_binary)))

#'ace' is known to stop working over long branch lengths, so we shink the tree.
bantu_phylo_pruned.scaled<-bantu_phylo_pruned
bantu_phylo_pruned.scaled$edge.length <- bantu_phylo_pruned$edge.length/100

# Add ancestral states probabilities for each node, using 'ace':
fitER<-ace(bantudata$trance_binary,bantu_phylo_pruned.scaled,model="ER",type="discrete") 

#ER here stands for 'equal-rates'; ace can estimate ancestral states under a few other models, and the packages 'diversitree' has even more flexibility. 
#Note diversitree requires  ultrametric trees, or trees with all of the tips being the same age, i.e., no extinct languages.


# We first set margins for the graph and fontsize for the tiplabels
par(mai = c(0,0,0.3,0), cex=0.9)

# Then we plot the tree and tilabels
# The tiplabels are offset from the trees
# x.lim and y.lim are used to position the tree and tiplabels against the background
plot(bantu_phylo_pruned,type="phylogram",show.tip.label=T, 
	label.offset = 0.02, x.lim =0.1, y.lim =50)

# Let's add the title and adjust its position
title(main = "Bantu", adj=0.3)

# Next we add the tiplabels
tiplabels(pie=to.matrix(bantudata$trance_binary,sort(unique(bantudata$trance_binary))),piecol=cols,cex=0.4)

# Then we add the labelnodes (the reconstructed ancestral states)
nodelabels(node=1:bantu_phylo_pruned$Nnode+Ntip(bantu_phylo_pruned),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)

# Next we add a horizontal thermo-type bar for rainfall and adjust its position and width.
# The width needs to be controlled and we need to leave enough space for the 
# horizontal bar between the tiplabels and the tipsymbols. For this, we need to
# carefully adjust label.offset and x.lim above.
# I tried to change the colors but was not able to do that yet.
tiplabels(thermo = bantudata[,21], col=cols[[1]], adj=0.51, horiz=T, width = 0.015)

# Add title for "Rain constancy"
# The first two numbers are x- and y- coordinates, cex is font-size
# pos left-aligns the text to the coordinates
# pch 16 is a symbol for the filled circle
text(0.067,50.5,label = "Rain constancy", cex=0.8)
text(0.056,51.6,label = "Absent", cex=0.8, pos = 2)
text(0.056,50.2,label = "Present", cex=0.8, pos = 2)
points(0.0565, 51.6, pch = 16, cex = 1.35, col = "black")
points(0.0565, 50.2, pch = 16, cex = 1.35, col = "red")
