#This function plots a UMAP, partitions cells and learns the graphs
#Preprocessing has to be done before this step
UMAPLearnGraphFixPartition_3alpha <- function(cds){
	cds <- reduceDimension(cds,reduction_method = 'UMAP')
	cds <- partitionCells(cds,resolution = 0.0000001)
	cds <- learnGraph(cds,RGE_method = 'SimplePPT')
	return(cds)
}
