
##### Integration of scRNAseq data

# Upload packages

library(Seurat)
library(dplyr)
library(RColorBrewer)
library(reticulate)

set.seed(111)


# Working directory
setwd(
  "")
getwd()

# Upload data with good quality cells

A1 <- readRDS("A1.rds")
A2 <- readRDS("A2.rds")
A3 <- readRDS("A3.rds")
A4 <- readRDS("A4.rds")
A5 <- readRDS("A5.rds")
A6 <- readRDS("A6.rds")
A7 <- readRDS("A7.rds")
A8 <- readRDS("A8.rds")
A9 <- readRDS("A9.rds")
A10 <- readRDS("A10.rds")

# Colors
#for subclustering #RColorBrewer
cluster.col=brewer.pal(name = "Paired",n=12)
cluster.col=cluster.col[c(2,4,6,8,10,12,1)]
#for gradients
gradient.col = rev(brewer.pal(n = 11, name = "RdYlBu"))
half.gradient.col = brewer.pal(n = 9, name = "YlOrRd")
col.ramp<-colorRampPalette(gradient.col)

#### Integration Seurat v5 RPCA
# Normalization with SCT method
## This is quite important to make the clusters (if stromal cells are not present in all samples, their features - genes might not be identified for integration)
Patients.list <- list(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10)

Patients.list <- lapply(Patients.list, FUN = SCTransform, vars.to.regress = "percent.mt")
Integration.features <- SelectIntegrationFeatures(object.list= Patients.list, nfeatures= 3000)

Patients.list <- merge(Patients.list[[1]], y= Patients.list[2:length(Patients.list)], merge.data = TRUE)

VariableFeatures.Patients.merged <- Integration.features

Patients.list <- RunPCA(object = Patients.list, assay = "SCT", features = VariableFeatures.Patients.merged, npcs = 30)
Patients.list <- RunUMAP(object = Patients.list, dims = 1:30)
DimPlot(Patients.list, group.by = "Patient", reduction = "umap", pt.size = .1, alpha = 0.2)

##### Integration With RPCA method

Patients.list <- IntegrateLayers(
  object = Patients.list, method = RPCAIntegration,
  normalization.method = "SCT", new.reduction = "integrated.rpca",
  verbose = FALSE
)

# re-join layers after integration
Patients.integrated[["RNA"]] <- JoinLayers(Patients.integrated[["RNA"]])
Patients.integrated <- FindNeighbors(Patients.list, reduction = "integrated.rpca", dims = 1:30)
Patients.integrated <- FindClusters(Patients.integrated, resolution = 2, cluster.name = "rpca_clusters")
Patients.integrated <- RunUMAP(Patients.integrated, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

DimPlot(Patients.integrated, group.by = "rpca_clusters", reduction = "umap.rpca", pt.size = .1, alpha = 0.3, label = T)

saveRDS(Patients.integrated, file = "Patients_integrated.rds")



