

### Code to generate figures in Figure 3
# Spatial transcriptomics analysis
# Article. Leukemia Journal
# Upload packages

library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(RColorBrewer)
library(spacexr)
library(dplyr)
library(pals)
set.seed(111)

options(future.globals.maxSize = 8e+09)

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

# Colors

colors.celltypes <- c("Adipogenic progenitors"="#E31A1C" , "Early mesenchymal progenitors"= "#B5D33D", 
                      "Osteogenic-lineage cells"="#FF7F00", "Vascular endothelial cells"="#CAB2D6",
                      "Leukemic cells"= "#1F78B4", "T/ NK"= "#FB9A99",
                      "Monocytes/ Macrophages" ="#A63603", "Dendritic cells"= "#FDDACB",
                      "HSPCs"="#D579BA","Granulocytes"="#FC9778")


colors.celltypes2 <- c("Adipogenic progenitors"="#E31A1C" , "Early mesenchymal progenitors"= "#B5D33D", 
                       "Osteogenic-lineage cells"="#FF7F00", "Vascular endothelial cells"="#CAB2D6",
                       "Leukemic cells"= "#1F78B4")

colors.celltypes3 <- c("Adipogenic progenitors"="#E31A1C" , "Early mesenchymal progenitors"= "#B5D33D", 
                       "Osteogenic-lineage cells"="#FF7F00", "Vascular endothelial cells"="#CAB2D6")

colors.celltypes4 <- c("Adipogenic progenitors"="#E31A1C" , "Early mesenchymal progenitors"= "#B5D33D", 
                       "Osteogenic-lineage cells"="#FF7F00")

# Working directory
setwd(
  '')
getwd()


# Read object
xenium.obj <- readRDS('')
levels(xenium.obj)

#Reorder identities
levels(xenium.obj) <- c(
  "Adipogenic progenitors","Early mesenchymal progenitors","Osteogenic-lineage cells",
  "Vascular endothelial cells",
  "HSPCs", "T/ NK", "Monocytes/ Macrophages", "Dendritic cells", "Granulocytes",
  "Leukemic cells")

xenium.obj$celltypes <- Idents(xenium.obj)
Idents(xenium.obj) <- "celltypes"

# Figure 3 - Spatial -- Supplementary

DimPlot(xenium.obj,group.by = ,label = T, cols = colors.celltypes, reduction = "umap", pt.size = 0.6, alpha=0.8,label.size = 5,
) + NoAxes()


DimPlot(
  xenium.obj,
  label = F,
  group.by = "celltypes",
  cols = colors.celltypes,
  reduction = "umap",
  pt.size = 0.6,
  label.size = 5,
  repel = T,
  shuffle =T 
) +
  NoAxes() +
  theme_leukemia() +
  theme_no_axes()

# Save the plot
ggsave("Fig_3_Supplementary_UMAP_Spatial_AllPopulations_DimPlot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 3B - Spatial -- Supplementary

# Select 3 - 4 markers per celltype based on differential expression analysis
levels(xenium.obj)

canonical_markers <- c(
  "LPL", "FABP4","CDH11", "CXCL12", # Adipogenic progenitors
  "COL1A1", "VCAN","FBN1","SPARC","TNC", # Early mesenchymal progenitors
  "RUNX2", "CHAD","IBSP", "BMP2","FGFR1",# Osteogenic progenitors
  "MYH11","CNN1",# Pericytes - Smooth muscle cells
  "ENG", "PECAM1", "VWF", "CLEC14A",  # General vascular endothelial cells
  "LYVE1","SELE", # Sinusoidal
  "ANGPT2","SOX17",# Arteorioal
  "SPINK2", # HSPCs
  "CD3E","CD28","GZMB","PRF1","NKG7", # T/ NK cells
  "CD14", "CD68", "FCGR3A","CD163", # Monocytes/Macrophages
  "IL3RA","MRC1","LILRA4","CD1C","CLEC10A", # Dendritic cells
  "MMP8", "CEACAM8", # Granulocytes
  "CD79A", "TCL1A", "CD19", "BANK1" # Leukemic cells
)

### Perform this dotplot as a heatmap with the average expression per cluster. Select better markers (2 per cluster)
Xenium_average <- AverageExpression(xenium.obj, assays = "Xenium", return.seurat = F, slot = "counts")
zscore <- scale(t(Xenium_average$Xenium))
zscore[zscore > 2.5] <- 2.5
zscore[zscore < -2.5] <- -2.5

colors_pseudobulk <- c("white", "#d8f3dc", "#52b788", "#1b4332")
                      
library(ComplexHeatmap)
levels(xenium.obj)
levels_xenium <- levels(xenium.obj)
Heatmap(t(zscore[levels_xenium, intersect(canonical_markers, rownames(Xenium_average$Xenium))]),  
        cluster_rows = T, cluster_columns = F, col= colors_pseudobulk, column_names_rot = 55,
        column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8))

# Save the plot
pdf("Fig_3B_Supplementary_Heatmap_Spatial_AllPopulations.pdf", width = 6, height = 7)

Heatmap(t(zscore[levels_xenium, intersect(canonical_markers, rownames(Xenium_average$Xenium))]),  
        cluster_rows = T, cluster_columns = F, col= colors_pseudobulk, column_names_rot = 55,
        column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8))

dev.off()

### Figure 3. Spatial Stroma
# Subset stromal clusters and plot cells (UMAP, top markers in heatmap, one celltype per plot (zoom in blood vessel and bone) )
levels(xenium.obj)

### First plot stromal populations
levels(xenium.obj)
Stromal <- subset(xenium.obj, idents = c("Osteogenic-lineage cells", "Early mesenchymal progenitors", "Adipogenic progenitors"))
Stromal <- SCTransform(Stromal, assay = "Xenium")
Stromal <- RunUMAP(Stromal, dims = 1:30, reduction.name = "umap.stroma")

DimPlot(
  Stromal,
  label = F,
  group.by = "celltypes",
  cols = colors.celltypes,
  reduction = "umap.stroma",
  pt.size = 0.6,
  label.size = 5,
  repel = T,
  shuffle =T 
) +
  NoAxes() +
  theme_leukemia() +
  theme_no_axes()

# Save the plot
ggsave("Fig_3_UMAP_Spatial_Stromal_Populations_DimPlot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

## Top markers stromal populations

Cluster.markers_Stroma <- FindAllMarkers(Stromal, logfc.threshold = 0.1, min.pct = 0.01, 
                                         only.pos = T, test.use = "MAST")

table(Cluster.markers_Stroma$cluster)

# top marker genes that are significant up
Cluster.markers_Stroma <- subset(Cluster.markers_Stroma, c(Cluster.markers_Stroma$avg_log2FC >0 & Cluster.markers_Stroma$p_val_adj < 0.05))

# Sort markergene results by logFC and clusternumber
Cluster.markers_Stroma <- arrange(Cluster.markers_Stroma,avg_log2FC)
Cluster.markers_Stroma <- arrange(Cluster.markers_Stroma,cluster)

writexl::write_xlsx(Cluster.markers_Stroma, "Xenium_Stroma_ALL_Bone_markers_DEGs_MAST.xlsx")

#### Feature plot 

FeaturePlot(
  Stromal,
  features = c("LPL", "PPARG",
               "IBSP", "FGFR1",
               "COL1A1", "TNC"),
  reduction = "umap.stroma",
  pt.size = 0.6,
  label.size = 5
)

#### Heatmap
canonical_markers <- c(
  "LPL", "FABP4","CDH11", "CXCL12",# Adipogenic progenitors
  "MYH11","CNN1",# Pericytes - Smooth muscle cells
  "COL1A1", "VCAN","FBN1","SPARC","TNC", # Early mesenchymal progenitors
  "RUNX2", "CHAD","IBSP", "BMP2","FGFR1"# Osteogenic progenitors
)

### Perform this dotplot as a heatmap with the average expression per cluster. Select better markers (2 per cluster)
Xenium_average <- AverageExpression(Stromal, assays = "Xenium", return.seurat = F, slot = "counts")
zscore <- scale(t(Xenium_average$Xenium))
zscore[zscore > 2.5] <- 2.5
zscore[zscore < -2.5] <- -2.5

colors_pseudobulk <- c("white", "#d8f3dc", "#52b788", "#1b4332")

library(ComplexHeatmap)
levels(Stromal)
levels_xenium <- levels(Stromal)
Heatmap(t(zscore[levels_xenium, intersect(canonical_markers, rownames(Xenium_average$Xenium))]),  
        cluster_rows = T, cluster_columns = F, col= colors_pseudobulk, column_names_rot = 55,
        column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8))

# Save the plot
pdf("Fig_3B_Supplementary_Heatmap_Spatial_StromalPopulations.pdf", width = 6, height = 7)

Heatmap(t(zscore[levels_xenium, intersect(canonical_markers, rownames(Xenium_average$Xenium))]),  
        cluster_rows = T, cluster_columns = F, col= colors_pseudobulk, column_names_rot = 55,
        column_names_gp = grid::gpar(fontsize = 8), row_names_gp = grid::gpar(fontsize = 8))

dev.off()

### Figure 3. Barplot of cell populations
table(xenium.obj$celltypes)

# Generate data frame with counts and percentages
cell_df <- as.data.frame(table(xenium.obj$celltypes))
colnames(cell_df) <- c("CellType", "Count")
cell_df <- cell_df %>%
  mutate(Percent = Count / sum(Count) * 100)

# Create named vector of legend labels with percentages
legend_labels <- setNames(
  paste0(cell_df$CellType, " (", sprintf("%.1f", cell_df$Percent), "%)"),
  cell_df$CellType
)

# Plot
ggplot(cell_df, aes(x = "All Cells", y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(
    values = colors.celltypes,
    labels = legend_labels
  ) +
  labs(
    title = "Cell Type Composition",
    x = "",
    y = "Percentage of Cells",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )+ theme_leukemia()+ theme_no_axes()

# Save the plot
ggsave("Fig_3_Supplementary_Barplot_Spatial_All_Populations.pdf", width = 8, height = 3.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

getwd()

### Barplot showing only the abundance of the three stromal populations

# Custom colors for the 3 selected cell types
colors.selected <- c(
  "Adipogenic progenitors" = "#E31A1C",
  "Early mesenchymal progenitors" = "#B5D33D",
  "Osteogenic-lineage cells" = "#FF7F00"
)

# Filter the Seurat object for selected cell types
selected_types <- c("Adipogenic progenitors", "Early mesenchymal progenitors", "Osteogenic-lineage cells")
cell_df <- as.data.frame(table(xenium.obj$celltypes))
colnames(cell_df) <- c("CellType", "Count")

# Filter and compute percentages
cell_df <- cell_df %>%
  filter(CellType %in% selected_types) %>%
  mutate(Percent = Count / sum(Count) * 100)

# Format legend labels with percentages
legend_labels <- setNames(
  paste0(cell_df$CellType, " (", sprintf("%.1f", cell_df$Percent), "%)"),
  cell_df$CellType
)

# Plot
ggplot(cell_df, aes(x = "Selected Cells", y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = colors.selected, labels = legend_labels) +
  labs(
    title = "Selected Cell Type Composition",
    x = "",
    y = "Percentage of Cells",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  ) + theme_leukemia() + theme_no_axes()


# Save the plot
ggsave("Fig_3_Barplot_Spatial_Stromal_Populations.pdf", width = 8, height = 3.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()


##### Figure 3. Barplot of cell populations in the different niches
table(xenium.obj$niches, xenium.obj$celltypes)

library(tidyr)

# Reconstruct the table from data
niche_cell_table <- matrix(c(
  2173, 797, 0, 848, 11034, 11775, 7605, 8983, 2690, 62412,
  5586, 1956, 27, 1713, 2528, 4117, 3632, 2031, 1076, 15006,
  1875, 451, 3, 368, 4166, 13608, 16173, 3262, 2394, 24139,
  700, 241, 998, 154, 908, 1454, 1314, 709, 382, 6051
), 
nrow = 4, byrow = TRUE)

colnames(niche_cell_table) <- c(
  "Adipogenic progenitors", "Early mesenchymal progenitors", "Osteogenic-lineage cells",
  "Vascular endothelial cells", "HSPCs", "T/ NK", "Monocytes/ Macrophages",
  "Dendritic cells", "Granulocytes", "Leukemic cells"
)
rownames(niche_cell_table) <- paste0("Niche ", 1:4)

# Convert to long format for ggplot2
df_long <- as.data.frame(niche_cell_table) %>%
  mutate(Niche = rownames(.)) %>%
  pivot_longer(-Niche, names_to = "CellType", values_to = "Count") %>%
  group_by(Niche) %>%
  mutate(Percent = Count / sum(Count) * 100)

# Desired legend order
legend_order <- c(
  "Adipogenic progenitors", "Early mesenchymal progenitors", "Osteogenic-lineage cells",
  "Vascular endothelial cells", "HSPCs", "T/ NK", "Monocytes/ Macrophages",
  "Dendritic cells", "Granulocytes", "Leukemic cells"
)

# Apply factor level order to CellType
df_long$CellType <- factor(df_long$CellType, levels = legend_order)

# Plot
ggplot(df_long, aes(x = Niche, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = colors.celltypes) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Cell Type Composition by Niche",
    x = "Niche",
    y = "Percentage of Cells",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold")
  )+ theme_leukemia()

# Save the plot
ggsave("Fig_3_Supplementary_Barplot_Niche_Distribution_All_Populations.pdf", width = 6, height = 7, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

##### Figure 3. Barplot of cell populations in the different niches considering only stromal and vascular endothelial cells
# Reconstruct the full table
niche_cell_table <- matrix(c(
  2173, 797, 0, 848, 11034, 11775, 7605, 8983, 2690, 62412,
  5586, 1956, 27, 1713, 2528, 4117, 3632, 2031, 1076, 15006,
  1875, 451, 3, 368, 4166, 13608, 16173, 3262, 2394, 24139,
  700, 241, 998, 154, 908, 1454, 1314, 709, 382, 6051
), 
nrow = 4, byrow = TRUE)

colnames(niche_cell_table) <- c(
  "Adipogenic progenitors", "Early mesenchymal progenitors", "Osteogenic-lineage cells",
  "Vascular endothelial cells", "HSPCs", "T/ NK", "Monocytes/ Macrophages",
  "Dendritic cells", "Granulocytes", "Leukemic cells"
)
rownames(niche_cell_table) <- paste0("Niche ", 1:4)

# Convert to long format
df_long <- as.data.frame(niche_cell_table) %>%
  mutate(Niche = rownames(.)) %>%
  pivot_longer(-Niche, names_to = "CellType", values_to = "Count") %>%
  group_by(Niche) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup() %>%
  filter(CellType %in% c(
    "Adipogenic progenitors", 
    "Early mesenchymal progenitors", 
    "Osteogenic-lineage cells", 
    "Vascular endothelial cells"
  ))

# Set legend order and colors
legend_order <- c(
  "Adipogenic progenitors", 
  "Early mesenchymal progenitors", 
  "Osteogenic-lineage cells", 
  "Vascular endothelial cells"
)

colors.celltypes <- c(
  "Adipogenic progenitors" = "#E31A1C",
  "Early mesenchymal progenitors" = "#B5D33D",
  "Osteogenic-lineage cells" = "#FF7F00",
  "Vascular endothelial cells" = "#CAB2D6"
)

df_long$CellType <- factor(df_long$CellType, levels = legend_order)

# Plot
ggplot(df_long, aes(x = Niche, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = colors.celltypes) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Stromal Cell Type Composition by Niche",
    x = "Niche",
    y = "Percentage of Cells",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold")
  ) + theme_leukemia()

# Save the plot
ggsave("Fig_3_Barplot_a_Niche_Distribution_Stroma_Populations.pdf", width = 6, height = 7, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()


### Statistical analysis of niche distribution
# Reconstruct the 4x4 matrix with only the four stromal cell types

ct <- matrix(c(
  2173, 797,   0,  848,
  5586, 1956, 27, 1713,
  1875, 451,   3,  368,
  700,  241, 998,  154
), 
nrow = 4, byrow = TRUE)

colnames(ct) <- c(
  "Adipogenic progenitors", 
  "Early mesenchymal progenitors", 
  "Osteogenic-lineage cells", 
  "Vascular endothelial cells"
)
rownames(ct) <- paste0("Niche ", 1:4)

# Perform Chi-squared test
chisq_res <- chisq.test(ct)

# Extract standardized residuals
residuals_mat <- chisq_res$stdres

res_df <- as.data.frame(as.table(residuals_mat))
colnames(res_df) <- c("Niche", "CellType", "StdResidual")

# Plot heatmap of residuals
ggplot(res_df, aes(x = CellType, y = Niche, fill = StdResidual)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-max(abs(res_df$StdResidual)), max(abs(res_df$StdResidual))),
                       name = "Std. Residual") +
  geom_text(aes(label = round(StdResidual, 2)), size = 4) +
  theme_minimal() +
  labs(
    title = "Chi-squared Residuals (Post-hoc Test)",
    x = "Cell Type",
    y = "Niche"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +theme_leukemia()

# Save the plot
ggsave("Fig_3_Barplot_Statistics_Niche_Distribution_Stroma_Populations.pdf", width = 8, height = 7, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Compute p-values from residuals
pvals <- 2 * (1 - pnorm(abs(residuals_mat)))

# Adjust for multiple testing
pvals_adj <- p.adjust(pvals, method = "fdr")

# View as matrix
matrix(pvals_adj, nrow = 4, dimnames = dimnames(ct))



#### Where do cells of this type are located? do they form different niches?

subset.obj <- subset(xenium.obj, subset = celltypes %in% c("Vascular endothelial cells",  "Early mesenchymal progenitors", "Adipogenic progenitors","Osteogenic-lineage cells"))

table(subset.obj$niches, subset.obj$celltypes)

# Reconstruct table from subset.obj (without Leukemic cells)
ct <- matrix(c(
  2173, 797,   0,  848,
  5586, 1956, 27, 1713,
  1875, 451,   3,  368,
  700,  241, 998,  154
), 
nrow = 4, byrow = TRUE)

colnames(ct) <- c(
  "Adipogenic progenitors", 
  "Early mesenchymal progenitors", 
  "Osteogenic-lineage cells", 
  "Vascular endothelial cells"
)
rownames(ct) <- paste0("Niche ", 1:4)

# Convert to long format for ggplot2
df_long <- as.data.frame(ct) %>%
  mutate(Niche = rownames(.)) %>%
  pivot_longer(-Niche, names_to = "CellType", values_to = "Count") %>%
  group_by(CellType) %>%
  mutate(Percent = Count / sum(Count) * 100)

# Reorder cell types for x-axis
cell_order <- c(
  "Vascular endothelial cells",
  "Adipogenic progenitors",
  "Early mesenchymal progenitors",
  "Osteogenic-lineage cells"
)
df_long$CellType <- factor(df_long$CellType, levels = cell_order)

# Define pastel niche colors
pastel_colors <- c(
  "Niche 1" = "#FDBBD2",  # pastel pink
  "Niche 2" = "#B3E2CD",  # pastel teal
  "Niche 3" = "#FEE0B6",  # pastel orange
  "Niche 4" = "#D0BBFF"   # pastel purple
)

# Create grouped barplot
ggplot(df_long, aes(x = CellType, y = Percent, fill = Niche)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = pastel_colors) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Distribution of Stromal Cell Types Across Niches",
    x = "Cell Type",
    y = "Percentage per Cell Type",
    fill = "Niche"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )+ theme_leukemia()

# Save the plot
ggsave("Fig_3_Barplot_b_Niche_Distribution_Stroma_Populations.pdf", width = 8, height = 7, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()


#### Pseudobulk heatmap with the identities
# Visually assess the correspondence between the clusters in the Xenium object and the niches

tab <- table(cluster=xenium.obj$celltypes, label=xenium.obj$niches) 
tab
pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 

# Save the plot
pdf("Fig_3_Supplementary_Spatial_Niches_Heatmap_AllPopulations.pdf", width = 6, height = 6)

pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 

dev.off()






