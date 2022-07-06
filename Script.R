library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)

###Figure 3 ----
#Import Data From 10x Cell Ranger
matrix_dir = "/Volumes/Samsung_T5/SingleCell-TNBC-MetaAnalysis/Wu_EMBO/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat.set <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat.set) = barcode.names$V1
rownames(mat.set) = feature.names$V2

### Import $ Extract Meta Data 
setwd("/Volumes/Samsung_T5/SingleCell-TNBC-MetaAnalysis/Wu_EMBO/")
meta.data <- fread("Wu_EMBO_metadata.csv",data.table = FALSE)
meta.data <- as.data.frame(meta.data)
row.names(meta.data) <- meta.data[,1]
meta.data.1 <- as.data.frame(meta.data[,c("NAME","celltype_final")])
meta.data.1 <- meta.data.1[-1,]
meta.data.1 <- meta.data.1[-1]
colnames(meta.data.1) <- "celltype"

#Set up Seurat object 
TNBC.seurat <- CreateSeuratObject(counts = mat.set, meta.data = meta.data.1, min.cells = 1)
TNBC.seurat <- NormalizeData(TNBC.seurat, verbose = FALSE)
TNBC.seurat <- FindVariableFeatures(TNBC.seurat, selection.method = "vst", nfeatures = 2000)
TNBC.seurat <- ScaleData(TNBC.seurat, verbose = FALSE)
TNBC.seurat <- RunPCA(TNBC.seurat, npcs = 30, verbose = FALSE)
TNBC.seurat <- FindNeighbors(TNBC.seurat, reduction = "pca", dims = 1:20)
TNBC.seurat <- FindClusters(TNBC.seurat, resolution = 0.25)
TNBC.seurat <- RunUMAP(TNBC.seurat, dims = 1:20)
TNBC.seurat <- RunTSNE(TNBC.seurat, dims = 1:20, check_duplicates = FALSE)

# QC the Data Seurat
TNBC.seurat[["percent.mt"]] <- PercentageFeatureSet(TNBC.seurat, pattern = "^MT-")
VlnPlot(TNBC.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)

# Label clusters by metadata
Idents(TNBC.seurat) <- TNBC.seurat$celltype

#Visualization
TSNEPlot(TNBC.seurat, label = TRUE, pt.size=3, repel=TRUE)
ggsave("Fig3A.png",dpi=300)

RidgePlot(TNBC.seurat, features = c("SPARC","POSTN"), ncol = 2, log = TRUE, sort = FALSE)
ggsave("Fig3B.png",dpi=300)

#Stats
CAF.1 <- c("myCAFs","iCAFs","dPVL","imPVL",
         "Myoepithelial","Endothelial")

CAF.res.1 <- NULL
for (i in CAF.1) {
  for (j in CAF.1) {
    if (i!=j) CAF.res.1[[paste(i, j, sep="-")]] <- FindMarkers(TNBC.seurat, ident.1=i, ident.2=j, only.pos=T)  
  }
}
gene <- "SPARC"
CAF.res1.SPARC <- do.call(rbind, lapply(CAF.res.1, function(x) x[gene,]))

###Figure S6 ----
setwd("/Volumes/Samsung_T5/SingleCell-TNBC-MetaAnalysis/GSE118389/")

#Import Data
TNBC.matrix <- fread("GSE118389_counts_rsem.txt",data.table = FALSE)
row.names(TNBC.matrix) <- TNBC.matrix[,1]
TNBC.matrix2 <- TNBC.matrix[,-1]

#Set up Seurat object 
TNBC.seurat2 <- CreateSeuratObject(counts = TNBC.matrix2, min.cells = 1)
TNBC.seurat2 <- NormalizeData(TNBC.seurat2, verbose = FALSE)
TNBC.seurat2 <- FindVariableFeatures(TNBC.seurat2, selection.method = "vst", nfeatures = 2000)
TNBC.seurat2 <- ScaleData(TNBC.seurat2, verbose = FALSE)
TNBC.seurat2 <- RunPCA(TNBC.seurat2, npcs = 30, verbose = FALSE)
TNBC.seurat2 <- FindNeighbors(TNBC.seurat2, reduction = "pca", dims = 1:15)
TNBC.seurat2 <- FindClusters(TNBC.seurat2, resolution = 0.5)
TNBC.seurat2 <- RunUMAP(TNBC.seurat2, dims = 1:15)
TNBC.seurat2 <- RunTSNE(TNBC.seurat2, dims = 1:15, check_duplicates = FALSE)

# QC the Data Seurat
TNBC.seurat2[["percent.mt"]] <- PercentageFeatureSet(TNBC.seurat2, pattern = "^MT-")
VlnPlot(TNBC.seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)

# Label Clusters & Subset
TNBC.seurat2 <- RenameIdents(TNBC.seurat2,`0`="BadQC1",
                            `1`="Cancer_P1",`2`="Cancer_P2",`3`="BadQC2",
                            `4`="BadQC3",`5`="Cancer_P3",`6`="T_cells",
                            `7`="Cancer_P4",`8`="Macro_M2",`9`="CAF_A",
                            `10`="Cancer_P6",`11`="B_cells",`12`="Cancer_P5",
                            `13`="CAF_B",`14`="Endo",`15`="CAF_C")

TNBC.seurat2 <- subset(TNBC.seurat2, ident=c("Cancer_P1","Cancer_P2",
                                           "Cancer_P3","T_cells","Cancer_P4",
                                           "Macro_M2","CAF_A","Cancer_P5","Cancer_P6",
                                           "B_cells","CAF_B","Endo","CAF_C"))

levels(TNBC.seurat2) <- c("Cancer_P6","Cancer_P5",
                         "Cancer_P4","Cancer_P3","Cancer_P2","Cancer_P1",
                         "T_cells","Macro_M2","B_cells","Endo","CAF_C",
                         "CAF_B","CAF_A")

#Visualization
TSNEPlot(TNBC.seurat2, label = TRUE, pt.size=3, repel=TRUE)
ggsave("FigS6A.png",dpi=300)

RidgePlot(TNBC.seurat2, features = c("SPARC","POSTN"), ncol = 2, log = TRUE, sort = FALSE)
ggsave("FigS6B.png",dpi=300)

###Figure S7 ----
setwd("/Volumes/Samsung_T5/SingleCell-Breast-CAF-Fatima/")

#Load data
CAF.breast <- readRDS(file="S1_7patients_integrated.rds")

#Set up Seurat object
CAF.breast <- FindVariableFeatures(CAF.breast, selection.method = "vst", nfeatures = 10000)
CAF.breast <- ScaleData(CAF.breast, verbose = FALSE)
CAF.breast <- RunPCA(CAF.breast, npcs = 30, verbose = FALSE)
CAF.breast <- FindNeighbors(CAF.breast, reduction = "pca", dims = 1:20)
CAF.breast <- FindClusters(CAF.breast, resolution = 0.25)
CAF.breast <- RunUMAP(CAF.breast, dims = 1:20)
CAF.breast <- RunTSNE(CAF.breast, dims = 1:20, check_duplicates = FALSE)

#Fix identities
Idents(CAF.breast) <- CAF.breast$Final.class.Paper
CAF.breast <- RenameIdents(CAF.breast,`0`="ECM-myCAF",`1`="Detox-iCAF",
                           `2`="IL-iCAF",`3`="TGFβ-myCAF",`4`="Wound-myCAF",
                           `5`="IFNγ-iCAF",`6`="IFNαβ-myCAF",`7`="Acto-myCAF")

#Visualization
RidgePlot(CAF.breast, features = c("SPARC","POSTN"), ncol = 2, log = FALSE, sort = FALSE)
ggsave("FigS7.png",dpi=300)

#Stats
CAF.2 <- c("ECM-myCAF","Detox-iCAF","IL-iCAF","TGFβ-myCAF",
"Wound-myCAF","IFNγ-iCAF","IFNαβ-myCAF","Acto-myCAF")

CAF.res.2 <- NULL
for (i in CAF.2) {
  for (j in CAF.2) {
    if (i!=j) CAF.res.2[[paste(i, j, sep="-")]] <- FindMarkers(CAF.breast, ident.1=i, ident.2=j, only.pos=T)  
  }
}
gene <- "SPARC"
CAF.res2.SPARC <- do.call(rbind, lapply(CAF.res.2, function(x) x[gene,]))
