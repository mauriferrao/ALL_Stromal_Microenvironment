
### Code to generate figures in Figure 5
# Niche factors
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
library(CellChat)
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


# Figure 4A
# Upload data
Stromal <- readRDS("")
levels(Stromal)

niche_factors <- c("CDH2","SPP1","GDF15","CXCL12", "IL7", "VCAM1", "KITLG",
                   "ANGPT1","IL6ST",
                   "CCN2",
                   "GAS6", "ITGA5", "ITGB1"
)

DotPlot(Stromal, group.by = "CellTypes_2Stromal",features = niche_factors,dot.scale = 6,assay = "SCT", scale = F, scale.min = 1, scale.max = 100)+ coord_flip()+
  scale_color_viridis_c(option = "C", direction = -1)+
  theme_leukemia()+
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 0)
  )

# Save the plot
ggsave("Fig_4A_Plot.pdf", width = 3, height = 8.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 4B
# Circle plot with cell-cell communication network
# Upload data

Patients.integrated_CellChat <- readRDS("Patients_integrated_CellChat.rds")

# Save the plot
pdf("Fig_4B_Plot.pdf", width = 8, height = 8)
netVisual_circle(Patients.integrated_CellChat@net$count, alpha.edge = 0.6,edge.curved = 0.4,edge.width.max = 0.5, shape = "none",vertex.label.cex = 0.8, color.use = ,  weight.scale = T,title.name = "")
dev.off()

# Figure 4C
# Communication axes scoring
# Scoring
scoreLeukemia <- function(sc.object){
  # Nourishing signals for scoring - all cell types
  LeukemiaNourishment = c("APP", "BTLA", "CD55", "CD99", "COL1A1",
                          "COL1A2", "COL4A1", "COL4A3", "COL6A1", "COL6A2",
                          "CTSG", "CXCL12", "EFNA5", "FGF2", "FGF7", "FLT3LG",
                          "FN1", "GDF15", "GZMA", "IGFBP3", "ITGAV", "ITGB1",
                          "LAMA1", "LAMA2", "LAMA4", "LAMB1", "LAMC1", "LGALS9",
                          "MDK", "MIF", "MPZL1", "NEGR1", "NRXN3", "PECAM1", 
                          "PGF", "PPIA", "PTN", "PTPRC", "RETN", "SEMA4D",
                          "SPP1", "TGFB1", "THBS1", "THY1", "VEGFA", "VEGFB")
  groups = list(LeukemiaNourishment)
  names(groups) = c("LeukemiaNourishment")
  
  #  scoring for the seurat object
  DefaultAssay(sc.object)="RNA"
  ctrl_genes = 35 #important
  
  for (gset in names(groups)){
    features = groups[gset]
    sc.object = AddModuleScore(object = sc.object, nbin = 22, features = features, name = gset, ctrl = ctrl_genes)
  }
  
  return(sc.object)
}#end of score


Patients.integrated <- scoreLeukemia(Patients.integrated)
DefaultAssay(Patients.integrated) <- "SCT"

Patients.integrated_withoutLeukemic <- subset(Patients.integrated, idents= "Leukemic cells",invert= T)

VlnPlot(Patients.integrated_withoutLeukemic, features = "LeukemiaNourishment1",pt.size = 0, cols = , sort = "decreasing", add.noise = F)+
  geom_boxplot(width=0.2,position = position_dodge(0.9),outlier.shape = NA,coef=0,lwd=0.3) +
  labs(y = "Leukemic cell interaction score")+ 
  labs(x = "")+ labs(title = NULL)+ NoLegend() + coord_flip()

# Save the plot
ggsave("Fig_4C_Plot.pdf", width = 7, height = 8, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 4D
###### Communication patterns
library(NMF)
library(ggalluvial)
# We run selectK to infer the number of patterns.
selectK(Patients.integrated_CellChat, pattern = "outgoing")
# Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 5
nPatterns = 5
dev.off()

# Save the plot
pdf("Fig_4D_Plot.pdf", width = 8, height = 8)
Patients.integrated_CellChat <- identifyCommunicationPatterns(Patients.integrated_CellChat, pattern = "outgoing", k = nPatterns, height = 12, font.size = 5)
dev.off()

# Figure 4E (Supplementary)
# river plot
# Save the plot
pdf("Fig_4E_Plot.pdf", width = 8, height = 8)
netAnalysis_river(Patients.integrated_CellChat, pattern = "outgoing")
dev.off()








