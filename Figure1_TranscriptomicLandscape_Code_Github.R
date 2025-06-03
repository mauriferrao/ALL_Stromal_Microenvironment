
### Code to generate figures in Figure 1
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
set.seed(111)

# Working directory
setwd(
  '')
getwd()

# Colors
colors.celltypes <- c("Stromal cells"="#FF7F00", 
                      "HSC"="#6A3D9A","LMPP"="#97D1F4","GMP"="#D579BA","Prog Mk"="#FC9778","Erythroblasts"="#E31A1C",
                      "pDC"="#40916c","cDC"="#BBD0E4",
                      "CD14 monocytes"="#A63603","CD16 monocytes"="#2ED1B5",
                      "CD4 naive"="#CAB2D6","Treg"="#AA1016","CD8 naive"="#A6CEE3", "CD8 memory"= "#8F9ECA", "CD8 effector"="#74c69d", "MAIT"="#FDDACB",
                      "GammaDelta"="#FDBF6F", "NK"= "#FB9A99",
                      "Mature B"= "#33A02C","Leukemic cells"="#B2DF8A")



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

#### Figure 1

# Figure 1C
# Upload data
Patients.integrated <- readRDS("")
levels(Patients.integrated)

#Reorder identities
levels(Patients.integrated) <- c("Stromal cells",
                                 "HSC","LMPP","GMP","Prog Mk","Erythroblasts",
                                 "pDC","cDC",
                                 "CD14 monocytes","CD16 monocytes",
                                 "CD4 naive","Treg","CD8 naive", "CD8 memory", "CD8 effector", "MAIT",
                                 "GammaDelta", "NK",
                                 "Mature B","Leukemic cells")
Patients.integrated$Annotation_level1 <- Idents(Patients.integrated)
Idents(Patients.integrated) <- "Annotation_level1"

# Subset leukemic niche object
Idents(Patients.integrated) <- "Condition"
Patients.integrated.Leukemic <- subset(Patients.integrated, idents= "Health", invert= T)


DimPlot(
  Patients.integrated.Leukemic,
  label = FALSE,
  group.by = "Annotation_level1",
  cols = colors.celltypes,
  reduction = "umap.rpca",
  pt.size = 0.1,
  label.size = 5
) +
  NoAxes() +
  theme_leukemia() +
  theme_no_axes()

# Save the plot
ggsave("Fig_1C_UMAP_DimPlot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 1D
# --- barplot of the contribution of cluster per patient
# How each patient contributed to each cluster?
breakdown<-table(Patients.integrated.Leukemic@meta.data$Annotation_level1, Patients.integrated.Leukemic@meta.data$patient)
breakdown<- breakdown[,c("A1", "A2", "A3","A4","A5","A6","A7","A8","A9")]
breakdown=t(breakdown)
breakdown
breakdown <- round(apply(breakdown, 2, function(x){x*100/sum(x)}),2)
breakdown.df = as.data.frame(breakdown)
breakdown.df = melt(t(breakdown.df))
breakdown.df

barplot.col=brewer.pal(name = "Set3",n=9)


ggplot(data = breakdown.df, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = barplot.col) +
  labs(
    x = "Cell Type",
    y = "Percentage of contribution [%]",
    fill = "Patient"  
  ) +
  theme_leukemia() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  )


# Save the plot
ggsave("Fig_1D_BarPlot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 1E
# Function - Dot Plots
dot_plot_leg <- function(seurat_object, features, group= group,lim=NULL, feats=NULL, min=-1.5, max=1.5){
  dot_plot <- DotPlot(object = seurat_object, dot.scale = 4,group.by= group, features=features,col="RdYlBu", col.min = min, col.max=max, scale.max = 50) + 
    theme(legend.position="right",axis.ticks.x = element_blank(), axis.text.x = element_text(size=10, angle=30, hjust=1), 
          axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    labs(x=feats, y=NULL) +  scale_y_discrete(limits=lim) + coord_flip()
  return(dot_plot)
}

canonical_markers <- c(
  "CD14", "LYZ",  # CD14_Mono
  "MS4A7","FCGR3A", # CD16_Mono
  "CD1C", "CLEC10A",  "IRF7",# cDC
  "TCF4", "LILRA4", # pDC 
  "CD34", "PROM1","AVP","SPINK2", # HSPC
  "MPO", "KIT", # GMP
  "FLT3", # LMPP
  "PF4", "ITGA2B", "GP9", # Prog_Mk
  "TCF7", "LEF1", "CD28", # CD4 naive T cells (CCR7 & SELL as well)
  "IL2RA", "FOXP3", "CTLA4", "CD27", # Tregs
  "CCR7", # CD8 naive T cells
  "GZMA","GZMB", "PRF1", "NKG7", "IFNG",  # CD8 effector
  "GZMK", # CD8_memory
  "CD8A", "CD8B", # CD8_general
  "TRDC", "TRGC1", "TRGC2", # gamma_delta_T
  "MR1","TBX21", # MAIT
  "CD3E", # General T Cell Markers
  "KLRK1", "CX3CR1", # NK-Like T Cells
  "NCAM1", "NCR1",# NK Cells
  "CD19", "MS4A1","CD22", # B-Cells (leukemic or not)
  "PTPRC", # General hematopoietic
  "HBD", "GYPA","GATA1", # Erythroblasts
  "CXCL12", "THY1", "FN1" # Stromal cells
)


dot_plot_leg(Patients.integrated, features = canonical_markers, group= "Annotation_level1") + 
  theme_leukemia()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# Save the plot
ggsave("Fig_1E_DotPlot.pdf", width = 6, height = 9, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()









