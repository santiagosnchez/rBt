root.tree <- function(phy,rt){
	require(phytools)
	node <- getMRCA(phy, grep(paste(rt), phy$tip.label))
	edge <- min(which.edge(phy, grep(paste(rt), phy$tip.label)))-1
	midedgeL <- phy$edge.length[edge]/2
	rphy <- reroot(phy, node.number=node, position=midedgeL)
	rphy
}