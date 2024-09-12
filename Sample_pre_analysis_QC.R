
# Pipeline for Seurat Object generation and Quality Control
## Selection of good quality cells - filter out low quality cells
### Pre-analysis before integration

# Upload packages

library(Seurat)
library(scDblFinder)
library(dplyr)

# Working directory
setwd(
  "")
getwd()

# Create Seurat object

A1 <- Read10X(data.dir="")
A1 <- CreateSeuratObject(counts = A1, project = "A1", min.cells = 3, min.features = 200)
A1

VlnPlot(A1, features = c("XIST", "UTY"), log= T)

### QC

A1$"percent.mt" <- PercentageFeatureSet(A1, pattern = "^MT-")
VlnPlot(A1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2, pt.size = 0.1)
FeatureScatter(A1, "nCount_RNA", "nFeature_RNA", pt.size = 0.5)

# Filter out doublets

A1 <- subset(A1, subset = nFeature_RNA > 200)
A1.SCE <- Seurat::as.SingleCellExperiment(A1)
Doublets.classification.A1 <- scDblFinder(A1.SCE)
table(Doublets.classification.A1$scDblFinder.class)
A1 <- A1[, Doublets.classification.A1$scDblFinder.class == "singlet"]
dim(A1)

### Selection of good quality cells

A1 <- subset(A1, subset = nFeature_RNA > 200 & percent.mt < 10)
A1

# SaveRDS
SeuratObject::saveRDS(A1, file = "A1.rds")

##### Same pipeline for the rest of the samples

