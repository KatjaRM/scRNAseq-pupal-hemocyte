library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(sctransform)

dataset1_raw <- read.csv(file.choose())
dataset2_raw <- read.csv(file.choose())
dataset3_raw <- read.csv(file.choose())

dataset1_raw <- dataset1_raw %>% column_to_rownames(var="X")
dataset2_raw <- dataset2_raw %>% column_to_rownames(var="X")
dataset3_raw <- dataset3_raw %>% column_to_rownames(var="X")

dataset1 <- CreateSeuratObject(counts = dataset1_raw, project = "dataset1")
dataset2 <- CreateSeuratObject(counts = dataset2_raw, project = "dataset2")
dataset3 <- CreateSeuratObject(counts = dataset3_raw, project = "dataset3")

dataset1[["percent.mt"]] <- PercentageFeatureSet(dataset1, pattern = "mt:")
dataset2[["percent.mt"]] <- PercentageFeatureSet(dataset2, pattern = "mt:")
dataset3[["percent.mt"]] <- PercentageFeatureSet(dataset3, pattern = "mt:")

dataset1 <- subset(dataset1, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
dataset2 <- subset(dataset2, subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & percent.mt < 5)
dataset3 <- subset(dataset3, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 5)

dataset1 <- SCTransform(dataset1, vars.to.regress = "percent.mt")
dataset2 <- SCTransform(dataset2, vars.to.regress = "percent.mt")
dataset3 <- SCTransform(dataset3, vars.to.regress = "percent.mt")

dataset1 <- RenameCells(dataset1, add.cell.id = "fix2")
dataset2 <- RenameCells(dataset2, add.cell.id = "fix3")
dataset3 <- RenameCells(dataset3, add.cell.id = "live")

dataset1@meta.data[, "dataset"] <- "dataset1"
dataset2@meta.data[, "dataset"] <- "dataset2"
dataset3@meta.data[, "dataset"] <- "dataset3"

hemocyte <- merge(x = dataset1, y = c(dataset2, dataset3),  merge.data = T, project = "pupal hemocyte")

DefaultAssay(hemocyte) <- "RNA"
list <- SplitObject(hemocyte, split.by = "dataset")
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
anchors <- FindIntegrationAnchors(object.list = list)

hemocyte <- IntegrateData(anchorset = anchors)

all.genes <- rownames(hemocyte)
hemocyte <- ScaleData(hemocyte, features = all.genes)
hemocyte <- RunPCA(hemocyte, features = VariableFeatures(object = hemocyte))

hemocyte <- JackStraw(hemocyte, num.replicate = 100)
hemocyte <- ScoreJackStraw(hemocyte, dims = 1:20)
ElbowPlot(hemocyte)

hemocyte <- FindNeighbors(hemocyte, dims = 1:16)
hemocyte <- RunUMAP(hemocyte, dims = 1:16) 

hemocyte <- RunTSNE(hemocyte, dims = 1:16) 

hemocyte <- FindClusters(hemocyte, resolution = Z)

DefaultAssay(hemocyte) <- "RNA"
markers <- FindAllMarkers(hemocyte, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

SetIdent(hemocyte, cells = WhichCells(hemocyte, idents = c('0', '3'), '3'))

hemocyte_list <- SplitObject(hemocyte, split.by = "dataset")

hemocyte_list$dataset1 <- SCTransform(hemocyte_list$dataset1, vars.to.regress = "percent.mt")
hemocyte_list$dataset2 <- SCTransform(hemocyte_list$dataset2, vars.to.regress = "percent.mt")
hemocyte_list$dataset3 <- SCTransform(hemocyte_list$dataset3, vars.to.regress = "percent.mt")

hemocyte_SCT <- merge(x = hemocyte_list$dataset2, y = c(hemocyte_list$dataset2, hemocyte_list$dataset3),  merge.data = T, project = "pupal hemocyte")
hemocyte_SCT@reductions <- hemocyte@reductions 

hemocyte_SCT <- PrepSCTFindMarkers(hemocyte_SCT)

hemocyte_named <- SetIdent(hemocyte_SCT, cells = WhichCells(hemocyte, idents = "11"), "Muscle")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "6"), "Neuron")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "14"), "PSC")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "8"), "Secretory-PL")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "12"), "AMP-PL")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "10"), "Adhesive-PL")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "13"), "Chitinase-PL")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "2"), "OxPhos-PL")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "1"), "Lsp-Bomanin-PL")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "4"), "Spermatid-Marker-PL")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "7"), "Osiris-PL")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "5"), "transitory PL-2")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "9"), "transitory PL-1")
hemocyte_named <- SetIdent(hemocyte_named, cells = WhichCells(hemocyte, idents = "3"), "undifferentiated PL")
