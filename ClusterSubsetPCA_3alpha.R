#This function takes a subset of the data based on a Boolean vector
#And then it does PCA and plots latent
ClusterSubsetPCA_3alpha <- function(cds,idents,genes=NULL){
	cell_type1_cells <- row.names(subset(pData(cds), orig.ident %in% idents))
	cds2 <- cds[,cell_type1_cells]
	cds2 <- setOrderingFilter(cds2,genes)
	cds2 <- estimateSizeFactors(cds2)
	cds2 <- preprocessCDS(cds2, num_dim=100)
	return(cds2)
}
