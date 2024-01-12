rm(list = ls())
setwd("../scRNA/")
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")

matrix_data <- read.table("./GSE149614_HCC.scRNAseq.S71915.count.txt.gz", sep="\t", header=T, row.names=1)
matrix_data[1:3,1:3]
table(substr(colnames(matrix_data),6,6) %in% "T" )
matrix_data = matrix_data[,substr(colnames(matrix_data),6,6) %in% "T"]
dim(matrix_data)

library(stringr)
Tissue = str_split(colnames(matrix_data),"_",simplify = T)
Tissue = as.data.frame(Tissue)
Tissue = as.character(Tissue$V1)
names(Tissue) = colnames(matrix_data)
seurat_obj <- CreateSeuratObject(counts = matrix_data,min.cells = 50,
                                 min.features = 800)

seurat_obj <- AddMetaData(object = seurat_obj, 
                          metadata = Tissue, 
                          col.name = 'Sample')

table(seurat_obj@meta.data$Sample)
pbmc = seurat_obj

minGene=800
maxGene=4000
pctMT=20

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


scRNA <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)

# Visualize QC metrics as a violin plot
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)




# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(scRNA, split.by = "Sample")
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", 
                            nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)


# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"
#save(immune.combined, file = "../data/HCC_immune.combined.Rdata")
pbmc = immune.combined


##如果内存不够，可以只对高变基因进行标准化
scale.genes <-  VariableFeatures(pbmc)
scRNA <- ScaleData(pbmc, features = scale.genes)

CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA))
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes) 

scRNAa <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
scRNAb <- ScaleData(scRNAa, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA))
DimPlot(scRNAb, reduction = "pca", group.by = "Phase")



scRNA <- RunPCA(scRNAb, features = VariableFeatures(scRNA)) 
ElbowPlot(scRNA, ndims=50, reduction="pca") 

pc.num=1:45

scRNA <- FindNeighbors(scRNA, dims = pc.num) 


res.used <- seq(0.1,1,by=0.1)
res.used
# Loop over and perform clustering of different resolutions 
sce = scRNA
for(i in res.used){
  sce <- FindClusters(object = sce, verbose = T, resolution = i)
}
# Make plot 
if(!require(clustree))install.packages("clustree")
library(clustree)
clus.tree.out <- clustree(sce) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")

clus.tree.out
scRNA <- FindClusters(scRNA, resolution = 0.5)
table(scRNA@meta.data$seurat_clusters)
#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(scRNA, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


scRNA <- RunTSNE(scRNA, dims = pc.num)
p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Sample")
p2 <- DimPlot(scRNA, reduction = "tsne", label = TRUE, repel = TRUE)
p1 + p2
save(immune.combined, file = "./data/HCC_immune.combined.Rdata")


rm(list = ls())
library(Seurat)
setwd("~/project/HCC/scRNA/")
load(file = "../data/HCC_immune.combined.Rdata")

# Visualize QC metrics as a violin plot
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(scRNA, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


library(Seurat)
library(ggplot2)
library(tidyverse)


# epi 10  # EPCAM SOX4 KRT19 KRT18
epi = c("EPCAM",  "KRT19", "KRT18",
        "DEFB1","CTSK","SOX4")  # 18

DotPlot(scRNA,  assay = "RNA",
        features = epi,
        group.by = 'seurat_clusters') + coord_flip()




Endo = c("CDH5","PECAM1","VWF",
         "CLDN5","FLT1","RAMP2",
         "SPARCL1","STC1","TM4SF1",
         "INSR")
DotPlot(scRNA,  assay = "RNA",
        features = Endo,
        group.by = 'seurat_clusters') + coord_flip()



myeloid = c("LYZ","CD68","FCGR3A","AIF1","RNASE1","C1QB","HLA-DRA")
DotPlot(scRNA,  assay = "RNA",
        features = myeloid,
        group.by = 'seurat_clusters') + coord_flip()




Bcell = c("CD79A","MS4A1","IGHM","IGHG3","BANK1","TNFRSF13C")
DotPlot(scRNA,  assay = "RNA",
        features = Bcell,
        group.by = 'seurat_clusters') + coord_flip()





Tcells = c("IFNG","PRF1","GZMK","GNLY","NKG7","GZMA",
           "LEF1","CCR7","SELL","TCF7","CCR6")
DotPlot(scRNA,  assay = "RNA",
        features = Tcells,
        group.by = 'seurat_clusters') + coord_flip()

NK = c("NKG7","GNLY","GZMB","CD7","KLRD1","NECAM1")
DotPlot(scRNA,  assay = "RNA",
        features = NK,
        group.by = 'seurat_clusters') + coord_flip()

FeaturePlot(scRNA, features = NK,
            reduction = "umap",  ncol=2,cols = c("lightgrey", "red"))




FeaturePlot(scRNA, features = Endo,
            reduction = "umap",  ncol=2,cols = c("lightgrey", "red"))

Fibro = c("DCN","LUM","COL1A1","COL1A2","THY1","RGS5","ACAT2","PDGFRB")
DotPlot(scRNA,  assay = "RNA",
        features = Fibro,
        group.by = 'seurat_clusters') + coord_flip()

FeaturePlot(scRNA, features = Fibro,
            reduction = "umap",  ncol=2,cols = c("lightgrey", "red"))





mast = c("ALB","APOA2","APOA1","AMBP","APOH","TTR")
DotPlot(scRNA,  assay = "RNA",
        features = mast,
        group.by = 'seurat_clusters') + coord_flip()



mast = c("IGLL1","MS4A2","GATA2")
DotPlot(scRNA,  assay = "RNA",
        features = mast,
        group.by = 'seurat_clusters') + coord_flip()

FeaturePlot(scRNA, features = mast,
            reduction = "umap",  ncol=2,cols = c("lightgrey", "red"))

DC = c("CLEC10A","CD1C","CLEC4C","PTCRA","CCR7","LAMP3","JCHAIN","TCF4","TCL1A")
DotPlot(scRNA,  assay = "RNA",
        features = DC,
        group.by = 'seurat_clusters') + coord_flip()


## 细胞注释信息添加到Seurat中
celltype = read.csv(file = "./main_celltype.csv")
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

scRNA@meta.data$main_celltype = scRNA@meta.data$celltype
## 汇总marker
marker=c(
  "FOXP3","IL2RA","CTLA4",
  "CD3D","CD2","CD3E","IL7R","CD3G","ITM2A",
  "LEF1","CCR7","SELL","TCF7","CCR6",
  "LYZ","CD68","FCGR3A","AIF1","C1QB","HLA-DRA",
  "KIT","MS4A2","GATA2",
  "DCN","LUM","COL1A1","COL1A2","THY1","RGS5","PDGFRB",
  "ALB","APOA2","APOA1","AMBP","APOH","TTR",
  "CDH5","PECAM1","FLT1","RAMP2","TM4SF1","INSR",
  
  "CLEC10A","CD1C","CD1E",
  "CD79A","MS4A1","BANK1"
)


DotPlot(scRNA,  assay = "RNA",
        features = marker,
        group.by = 'main_celltype') + coord_flip()


DimPlot(scRNA, reduction = "umap", group.by = "main_celltype",label = T)
## 细胞大类注释


#save(scRNA, file = "../data/HCC_immune.combined.Rdata")



## 绘制热图




marker=c(
  "FOXP3","IL2RA","CTLA4",
  
  "LEF1","CCR7","CCR6","IL7R",
  "LYZ","CD68","FCGR3A","AIF1","C1QB","HLA-DRA",
  "DCN","LUM","COL1A1","COL1A2","THY1","RGS5","PDGFRB",
  "ALB","APOA2","APOA1","AMBP","APOH","TTR",
  "CDH5","PECAM1","FLT1","RAMP2","TM4SF1","INSR",
  
  "CLEC10A","CD1C","PTCRA",
  "IFNG","PRF1","GZMK","GNLY","NKG7","GZMA",
  "CD79A","MS4A1","BANK1"
)

options(stringsAsFactors = FALSE)
heatmap = subset(scRNA, features = marker)
heatmap

DotPlot(scRNA,  assay = "RNA",
        features = marker,
        group.by = 'main_celltype') + coord_flip()






all_cell_types <- as.vector(heatmap@meta.data$main_celltype)
table(all_cell_types)

expr = heatmap@assays$RNA@counts
expr = as.matrix(expr)
dim(expr)
expr[1:3,1:3]
expr = t(scale(t(expr)))
mean_exp_eachCellType <- apply(expr, 1,
                               function(x)by(x, all_cell_types, mean))
dim(mean_exp_eachCellType)
ml = as.matrix(t(scale(mean_exp_eachCellType)))

ml[ml>2] = 2
ml[ml<-2] = -2
pheatmap(ml,
         cluster_cols = F,
         cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "black", "yellow"))(200))