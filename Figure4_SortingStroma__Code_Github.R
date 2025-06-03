

### Code to generate figures in Figure 4
# Sorting stromal cells and characterization
# Article. Leukemia Journal
# Upload packages

library(Seurat) 
library(dplyr)
library(RColorBrewer)
library(reticulate)
library(scater)
library(pals)
library(ggplot2)
library(Cairo)
library(ggpubr)
library(reshape2)
library(grid)
library(ComplexHeatmap)
library(tibble)
library(tidyr)
library(viridis)
set.seed(111)

# Working directory
setwd(
  '')
getwd()

# Colors
custom_colors <- c(
  "Adipogenic progenitors" =  "#FDCDAC",
  "Early mesenchymal progenitors" = "#B3E2CD" 
)

###

theme_leukemia <- function(base_size = 10, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      strip.text = element_text(size = 9, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(size = 0.5, color = "black"),
      axis.ticks = element_line(size = 0.5, color = "black"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

theme_no_axes <- function() {
  theme(
    plot.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )
}


# Figure 3A
# Upload data
Stromal <- readRDS("")
levels(Stromal)
markers_to_plot <- c("CD81","THY1","CD63", "NT5E", "ENG", "NGFR","ITGAV","CD164","MCAM", 
                     "CDH2","CALR",  "CD151",
                     "CADM1", "CD9", "ITGB1", "IGF1R","FGFR1",
                      "ICAM1", "DPP4","NCAM1",
                     "ITGA1","CD200","CD44",
                     "PDGFRA","PDGFRB",
                     "NCAM2", "LEPR", "CDH11", "VCAM1"# Other classical MSC markers
)

DotPlot(Stromal, group.by = "CellTypes_2Stromal",features = markers_to_plot,dot.scale = 6,assay = "SCT", scale = F, scale.min = 1, scale.max = 100)+ coord_flip()+
  scale_color_viridis_c(option = "C", direction = -1)+
  theme_leukemia()+
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 0)
  )

# Save the plot
ggsave("Fig_3A_Plot.pdf", width = 3, height = 8.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 3B
## Projection of sorted cells
library(SCP)

# Upload individual seurat object for each sorted stromal population
MSC1 <- readRDS("")
MSC2 <- readRDS("")
# Add info to metadata
MSC1$"Stromal population" <- "CD106_CDH11_StromalCells"
MSC2$"Stromal population" <- "CD90_CD81_StromalCells"

#### Downsample the objects to the same amount of cells before merging
MSC1<- subset(MSC1, downsample= 2300, seed= 111)
MSC2<- subset(MSC2, downsample= 2300, seed= 111)

# Normalization with SCT method
Sorted_Stromal.merged.list <- list(MSC1, MSC2)
Sorted_Stromal.merged <- lapply(Sorted_Stromal.merged.list, FUN = SCTransform, vars.to.regress = "percent.mt")
VariableFeatures.Stromal.merged <- SelectIntegrationFeatures(object.list= Sorted_Stromal.merged, nfeatures= 3000)
Sorted_Stromal.merged <- merge(Sorted_Stromal.merged[[1]], y= Sorted_Stromal.merged[2:length(Sorted_Stromal.merged)], merge.data = TRUE)
Sorted_Stromal.merged <- RunPCA(object = Sorted_Stromal.merged, assay = "SCT", features = VariableFeatures.Stromal.merged, npcs = 30)
Sorted_Stromal.merged <- RunUMAP(object = Sorted_Stromal.merged, dims = 1:30)
DimPlot(Sorted_Stromal.merged, group.by = "Stromal.population", reduction = "umap", pt.size = .5, alpha = 0.5)

# Project the datasets
srt_query <- RunKNNMap(srt_query = Sorted_Stromal.merged, srt_ref = Stromal, 
                       query_assay = "SCT", ref_assay = "SCT", ref_umap = "umap.rpca2",
                       ref_group = "CellTypes_2Stromal", features = VariableFeatures.Stromal.merged)

ProjectionPlot(
  srt_query = srt_query, 
  srt_ref = Stromal,
  query_group = "Stromal.population", 
  ref_group = "CellTypes_2Stromal",
  query_param = list(palette = "Pastel1", cells.highlight = TRUE),
  ref_param = list(palette = "Pastel2")) + theme_leukemia() + theme_no_axes()

# Save the plot
ggsave("Fig_3B_Plot.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Annotation with KNN
Sorted_Stromal.merged <- RunKNNPredict(
  srt_query = Sorted_Stromal.merged, srt_ref = Stromal,
  ref_group = "CellTypes_2Stromal", filter_lowfreq = 20, features = VariableFeatures.Stromal.merged
)

CellDimPlot(srt = Sorted_Stromal.merged, group.by = "KNNPredict_classification", label = TRUE)

# Ensure we have a clean data frame
df <- as.data.frame(Sorted_Stromal.merged@meta.data)

# Filter only for the sorted populations of interest
df_filtered <- df %>%
  filter(Stromal.population %in% c("CD106_CDH11_StromalCells", "CD90_CD81_StromalCells"))

# Count combinations
tab_counts <- as.data.frame(table(
  Stromal.population = df_filtered$Stromal.population,
  PredictedType = df_filtered$KNNPredict_classification
))

# Compute proportions within each Stromal.population
tab_counts <- tab_counts %>%
  group_by(Stromal.population) %>%
  mutate(
    proportion = Freq / sum(Freq),
    label = paste0(round(proportion * 100, 1), "%")
  )

# Statistics
# Create the contingency table from tab_counts
tab_counts2 <- tibble::tibble(
  Stromal.population = c("CD106_CDH11_StromalCells", "CD90_CD81_StromalCells", 
                         "CD106_CDH11_StromalCells", "CD90_CD81_StromalCells"),
  PredictedType = c("Adipogenic progenitors", "Adipogenic progenitors", 
                    "Early mesenchymal progenitors", "Early mesenchymal progenitors"),
  Freq = c(840, 210, 1460, 2090)
)

# Create a contingency table of frequencies
contingency_table <- tab_counts2 %>%
  pivot_wider(names_from = PredictedType, values_from = Freq, values_fill = list(Freq = 0)) %>%
  select(-Stromal.population) %>%
  as.matrix()

# Perform Chi-square test
chi_test <- chisq.test(contingency_table)

# Print the result
chi_test

# Your ggplot code with the p-value annotation
ggplot(tab_counts, aes(x = Stromal.population, y = proportion, fill = PredictedType)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Enrichment of progenitor types in sorted stromal populations",
    x = "Sorted population",
    y = "Proportion of predicted cell types",
    fill = "Predicted cell type"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  annotate("text", x = 1.8, y = 1.06, label = "p-value < 2.2e-16", size = 4, fontface = "italic", color = "black")+
  theme_leukemia() +  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, size = 10),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

# Save the plot
ggsave("Fig_3C_Plot.pdf", width = 5, height = 6, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()


# Figure 3D
## Cytotrace of sorted populations

# Upload data
Stromal <- readRDS('')

RidgePlot(Stromal,features = "Trayectory_CytoTrace",cols = c("#CCBB4480", "#228833"),
          sort = "deacreasing") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0,hjust =1,size = 10),axis.text.y = element_text(size = 10))+NoLegend()+ ggtitle("CytoTRACE") + ylab("")


# Save the plot
ggsave("Fig_3D_Plot.pdf", width = 9, height = 5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()



