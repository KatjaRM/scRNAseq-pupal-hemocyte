library(Seurat)
library(SCopeLoomR)
library(SCENIC)
library(pheatmap)
library(AUCell)

exprMatrix<- GetAssayData(hemocyte, slot="counts",assay = "RNA")
cellInfo <- data.frame(seuratCluster=Idents(hemocyte))

dir.create("SCENIC")
loom <- build_loom("hemocyte.loom", dgem=exprMatrix)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

options(width=200)
knitr::opts_knit$set(root.dir="SCENIC") 

loomPath <- "/hemocyte.loom"
library(SCopeLoomR)
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

scenicOptions <- initializeScenic(org="dmel", dbDir="/Users/katja/Documents/collaborations/hemocyte data/SCENIC hemocyte/", nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
org <- "dmel"
dbs <- defaultDbNames[[org]]
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

scenicOptions@settings$defaultTsne$perpl <- 50

cellInfo<- data.frame(seuratCluster=Idents(hemocyte))
nGene <- hemocyte@meta.data$nFeature_RNA
nGene <- as.data.frame(nGene)
nUMI <- hemocyte@meta.data$nCount_RNA
nUMI <- as.data.frame(nUMI)
cellInfo["nGene"] <- nGene
cellInfo["nUMI"] <- nUMI
names(cellInfo)[1] <- "CellType"

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)
saveRDS(savedSelections, "savedSelections.rds")
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_filtered)

regulonAUC <- readRDS(file.choose())

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "seuratCluster", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

#available as supplementary data
write.csv(topRegulators, "topRegulators.csv")
write.csv(regulonActivity_byCellType_Scaled, "regulonActivity_byCellType_Scaled.csv")

minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pheatmap::pheatmap(regulonActivity_byCellType_Binarized, fontsize_row=5, 
                   color = colorRampPalette(c("lightblue","white","magenta"))(100), breaks=seq(0, 1, length.out = 100),
                   cluster_rows = F, cluster_cols = F, border_color=NA)


tsneAUC(scenicOptions, aucType="AUC")

export2loom(scenicOptions, exprMat)

dr_coords <- Embeddings(hemocyte, reduction="umap")

par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, 'srp'), plots = "AUC")
close_loom(loom)