
### Code to generate figures in Figure 6
# Validation Targets. Cell-cell communication analysis
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

#

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

# Figure 6A, B, C
# Graphs and rawdata can be found in 
# '/Users/m.n.ferraoblanco/surfdrive - Mauricio Ferrao Blanco@surfdrive.surf.nl/HeterogeneousStromalCells/Pathway_Inhibition'
# Figure 6 E, F, G
# '/Users/m.n.ferraoblanco/surfdrive - Mauricio Ferrao Blanco@surfdrive.surf.nl/HeterogeneousStromalCells/ExVivoCoCulture_17February_2025/Analysis_DrugSensitivity'

# Figure 6D
# Upload data

Patients.integrated_CellChat <- readRDS('Patients_integrated_CellChat.rds')

# Compute the network centrality scores
Patients.integrated_CellChat <- netAnalysis_computeCentrality(Patients.integrated_CellChat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Plot interactions with VCAM1, ITGB1
netVisual_bubble(Patients.integrated_CellChat, signaling = c("VCAM", "COLLAGEN","SPP1", "FN1", "LAMININ", "TENASCIN"),color.heatmap = "viridis",direction = 1,sources.use = c("Adipogenic progenitors", "Early mesenchymal progenitors", "Leukemic cells"), targets.use = c("Adipogenic progenitors", "Early mesenchymal progenitors", "Leukemic cells"), angle.x = 30, vjust.x = 1, hjust.x = 1, sort.by.source = F,sort.by.source.priority = T, remove.isolate = T)
interactions_VCAM1_ITGB1 <- netVisual_bubble(Patients.integrated_CellChat, signaling = c("VCAM", "COLLAGEN","SPP1", "FN1", "LAMININ", "TENASCIN"),color.heatmap = "viridis",direction = 1,sources.use = c("Adipogenic progenitors", "Early mesenchymal progenitors", "Leukemic cells"), targets.use = c("Adipogenic progenitors", "Early mesenchymal progenitors", "Leukemic cells"), angle.x = 30, vjust.x = 1, hjust.x = 1, sort.by.source = F,sort.by.source.priority = T, remove.isolate = T, return.data = T)

# Extract the communication data
df <- interactions_VCAM1_ITGB1$communication

library(writexl)
write_xlsx(df,"interactions_VCAM1_ITGB1.xlsx")

library(readxl)
df <- read_excel('/Users/m.n.ferraoblanco/surfdrive - Mauricio Ferrao Blanco@surfdrive.surf.nl/Articles to submit/1_Distinct Mesenchymal Stromal Cell Populations Dictate Leukemic Cell Survival and are implicated in Relapse/Rebuttal_Leukemia/Figures/Figure6_PathwayInhibition/interactions_VCAM1_ITGB1.xlsx')  # Change filename accordingly
df

ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob.original)) +
  geom_point(size = 4) +  # Fixed size for all points
  scale_color_viridis_c(name = "Probability",
                        guide = guide_colorbar(title = "Probability",
                                               label = TRUE,
                                               barwidth = 10,
                                               barheight = 1.2,
                                               title.position = "top",
                                               label.position = "bottom",
                                               ticks = FALSE)) +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  guides(color = guide_colorbar(label = TRUE,
                                label.theme = element_text(angle = 0),
                                title = "Probability\n(min to max)"))+
  theme_leukemia()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


# Save the plot
ggsave("Fig_6D_BubblePlot_Interactions.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()







