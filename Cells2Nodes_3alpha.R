#This set of functions returns the node corresponding to the cell ids extracted by which
Cells2Nodes_3alpha <- function(cds,cell_ids){
	closest_vertex <- cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
	closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
	root_pr_nodes <- V(cds@minSpanningTree)$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
	return(root_pr_nodes)
}
