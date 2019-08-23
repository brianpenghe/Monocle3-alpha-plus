#This function converts a Seurat3 object to a Monocle3alpha object
#It carries over the metadata, idents and UMAP coordinates 
Seurat3ToMonocle3alpha <- function(Seurat){
	umap = Seurat[['umap']][[,]]
	obs = Seurat@meta.data
	obs["Seurat_ident"]<-as.vector(Seurat@active.ident)
	fData <- data.frame(gene_short_name = row.names(Seurat), 
                    row.names = row.names(Seurat))
	fd <- new('AnnotatedDataFrame', data = fData)
	X <- as(as.matrix(Seurat@assays$RNA@counts), 'sparseMatrix')
	pd = new('AnnotatedDataFrame',data=obs)
	mono = newCellDataSet(X, phenoData = pd,featureData = fd)
	umap = as.matrix(t(umap))
	colnames(umap) = rownames(obs)
	mono@reducedDimS = umap
	mono <- estimateSizeFactors(mono)
	mono <- estimateDispersions(mono)
	mono = partitionCells(mono)
	mono = learnGraph(mono,  RGE_method = 'SimplePPT')
	return(mono)
}
