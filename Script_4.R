#Data import:
#dataset = pupa: data published in this paper was prepared as described
#dataset = larva: Data from Cattenoz et al., 2020 (free larval hemocytes) was imported using Seurat::Read10X
#dataset = adult: Data from Li et al., 2022 (adult hemocytes) was imported from loom file (for details see: https://github.com/aaron-allen/VNC_scRNAseq/blob/master/R_preprocessing_brain_data.R)
#dataset = LG: Data from Cho et al., 2020 (lymph gland) was imported via read.table()
#Seurat objects were created using Seurat::CreateSeuratObject()

library(Seurat)
library(sctransform)
library(ggplot2)

hemocyte <- merge(x = pupa, y = c(larva, LG, adult),  merge.data = T, project = "hemocyte")

DefaultAssay(hemocyte) <- 'RNA'

hemocyte.list <- SplitObject(hemocyte, split.by = "dataset")

hemocyte.list <- lapply(X = hemocyte.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

anchors <- FindIntegrationAnchors(object.list = hemocyte.list, dims = 1:20, anchor.features = 2000)

hemocyte <- IntegrateData(anchorset = anchors)

DefaultAssay(hemocyte) <- "RNA"

hemocyte <- NormalizeData(hemocyte)

DefaultAssay(hemocyte) <- "integrated"
hemocyte <- FindVariableFeatures(hemocyte, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(hemocyte)
hemocyte <- ScaleData(hemocyte, features = all.genes)
hemocyte <- RunPCA(hemocyte, features = VariableFeatures(object = hemocyte))
hemocyte <- JackStraw(hemocyte, num.replicate = 100)
hemocyte <- ScoreJackStraw(hemocyte, dims = 1:20)
ElbowPlot(hemocyte, ndims = 50)

DimHeatmap(hemocyte, dims = 1:5, cells = 500, balanced = TRUE)
DimHeatmap(hemocyte, dims = 6:10, cells = 500, balanced = TRUE)
DimHeatmap(hemocyte, dims =11:15, cells = 500, balanced = TRUE)
DimHeatmap(hemocyte, dims =16:20, cells = 500, balanced = TRUE)
DimHeatmap(hemocyte, dims =21:25, cells = 500, balanced = TRUE)
DimHeatmap(hemocyte, dims =26:30, cells = 500, balanced = TRUE)
DimHeatmap(hemocyte, dims =31:35, cells = 500, balanced = TRUE)
DimHeatmap(hemocyte, dims =36:40, cells = 500, balanced = TRUE)

Y <- 26

hemocyte <- FindNeighbors(hemocyte, dims = 1:Y)

hemocyte <- RunUMAP(hemocyte, dims = 1:Y) 

UMAPPlot(hemocyte, group.by = 'dataset')

DimPlot(hemocyte, cells.highlight = WhichCells(pupa), pt.size = .5) + 
  scale_color_manual(labels = c("all other cells", "pupal cells"), values = c("#009E73", "magenta"))

DimPlot(hemocyte, cells.highlight = WhichCells(adult, idents = "crystal cell"), pt.size = .5) + 
  scale_color_manual(labels = c("all other cells", "adult crystal cell"), values = c("#009E73", "magenta"))

DimPlot(hemocyte, cells.highlight = WhichCells(larva, idents = "CC"), pt.size = .5) + 
  scale_color_manual(labels = c("all other cells", "larval CC"), values = c("#009E73", "magenta"))

DimPlot(hemocyte, cells.highlight = WhichCells(pupa, idents = "Secretory-PL"), pt.size = .5) + 
  scale_color_manual(labels = c("all other cells", "Secretory-PL"), values = c("#009E73", "magenta"))

DimPlot(hemocyte, cells.highlight = WhichCells(pupa, idents = "PSC"), pt.size = .5) + 
  scale_color_manual(labels = c("all other cells", "pupal PSC"), values = c("#009E73", "magenta"))

DimPlot(hemocyte, cells.highlight = WhichCells(adult, idents = "Hml+ with PSC markers"), pt.size = .5) + 
  scale_color_manual(labels = c("all cells", "adult Hml+ with PSC markers"), values = c("#009E73", "magenta"))

DimPlot(hemocyte, cells.highlight = WhichCells(larva, idents = "PL-Impl2"), pt.size = .5) + 
  scale_color_manual(labels = c("all cells", "larval PL-Impl2"), values = c("#009E73", "magenta"))


