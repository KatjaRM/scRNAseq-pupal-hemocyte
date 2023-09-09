library(monocle3)
library(Seurat)

my.so <- ProjectDim(hemocyte, reduction = "pca")

expression_matrix <- my.so@assays$RNA@counts

cell_metadata <- my.so@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

gene_annotation <- data.frame(gene_short_name = rownames(my.so@assays$RNA), row.names = rownames(my.so@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)

reducedDim(my.cds, type = "PCA") <- my.so@reductions$pca@cell.embeddings 
my.cds@preprocess_aux$prop_var_expl <- my.so@reductions$pca@stdev
plot_pc_variance_explained(my.cds)

my.cds@int_colData@listData$reducedDims$UMAP <- my.so@reductions$umap@cell.embeddings

my.cds@clusters$UMAP_so$clusters <- my.so@meta.data$gt_tp_cell_type_integrated_.0.9

my.cds <- cluster_cells(my.cds, reduction_method = "UMAP", resolution = 1e-3)

rownames(my.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(my.cds@int_colData@listData$reducedDims$UMAP) <- NULL

DimPlot(my.so, reduction = "umap")
DimPlot(my.so, reduction = "umap")
plot_cells(my.cds, color_cells_by="dataset", label_cell_groups=FALSE, cell_size = 1)

my.cds <- cluster_cells(my.cds, resolution = 0.002)
my.cds@clusters$UMAP$clusters <- hemocyte_named@active.ident
plot_cells(my.cds, label_cell_groups=FALSE, cell_size = 1)

my.cds <- learn_graph(my.cds, learn_graph_control = list(ncenter=480)) 
plot_cells(my.cds, label_cell_groups=FALSE, cell_size = 1)
my.cds<-order_cells(my.cds)
plot_cells(my.cds,color_cells_by="pseudotime",cell_size=2, label_roots = F, label_leaves = F, label_branch_points = F)
