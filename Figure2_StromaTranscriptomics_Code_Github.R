
### Code to generate figures in Figure 2
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


# Figure 2A
# Upload data
Stromal <- readRDS("")
DimPlot(Stromal,group.by = "CellTypes_2Stromal", reduction = "umap.rpca2", label = F, pt.size = 0.8, label.size = 8,cols = custom_colors,
) + NoAxes() + theme_leukemia() + theme_no_axes()

# Save the plot
ggsave("Fig_2A_Plot.pdf", width = 6.5, height = 5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 2B
# Dotplot
skeletal.cells.markers <- c("CXCL12", "THY1", "CD164", "PTGDS", "PRRX1","CD44","CDH2",# Multipotent Stromal Cells
                             "COL1A1", "COL1A2", "COL3A1","COL12A1", "FN1", "FBN1","BGLAP", "IBSP","VCAN","DCN",# Extracellular matrix proteins
                            "CHAD","PDGFRA",   "ALPL", "SPARCL1", "CDH11","RUNX2", "EBF3", "GAS6", "EBF1","SOX4","SOX9",# Osteochondroprogenitors
                            "CEBPD", "PPARG", "MGP", "LPL","LEPR", # adipogenic progenitors
                            "ACTA2", "FAP","CADM1","CALR","S100A4", "S100A6", "S100A8", "S100A9", "S100A10") # Fibroblasts


DotPlot(Stromal, group.by = "CellTypes_2Stromal",features = skeletal.cells.markers,dot.scale = 6,assay = "SCT", scale = F, scale.min = 1, scale.max = 100)+ coord_flip()+
  scale_color_viridis_c(option = "C", direction = -1)+
  theme_leukemia()+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 0)
  )


# Save the plot
ggsave("Fig_2B_Plot.pdf", width = 3, height = 8.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()


# Figure 2C
# Pathway activity - Progeny
library(progeny)

gradient.col = rev(brewer.pal(n = 8, name = "RdYlBu"))
col.ramp<-colorRampPalette(gradient.col)

DefaultAssay(Stromal) <- "progeny"
Stromal <- Seurat::ScaleData(Stromal, assay = "progeny")

progeny_scores_df <-
  as.data.frame(t(GetAssayData(Stromal, slot = "scale.data",
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

CellsClusters <- data.frame(Cell = names(Idents(Stromal)), 
                            CellType = as.character(paste0(Stromal$CellTypes_2Stromal)),
                            stringsAsFactors = FALSE)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

selected_pathways <- summarized_progeny_scores %>%
  semi_join(interesting_pathways,  by = "Pathway") %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

heat_data <- as.data.frame(t(selected_pathways))
heat_data$Pathway <- rownames(heat_data)

heat_data_long <- melt(heat_data, id.vars = "Pathway", variable.name = "Sample", value.name = "Score")

gradient.col = rev(brewer.pal(n = 8, name = "RdYlBu"))
col.ramp<-colorRampPalette(gradient.col)

# Plot heatmap
ggplot(heat_data_long, aes(x = Sample, y = Pathway, fill = Score)) +
  geom_tile(color = NA) +
  scale_fill_gradientn(colors = col.ramp(100), name = "Score") +
  labs(title = "Pathway activity score", x = NULL, y = NULL) +
  theme_leukemia() +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9)
  )+ theme_no_axes()

# Save the plot
ggsave("Fig_2C_Plot.pdf", width = 2, height = 3.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 2D
# Cytotrace
RidgePlot(Stromal,group.by = "CellTypes_2Stromal",features = "Trayectory_CytoTrace",cols = custom_colors,
          sort = "deacreasing") + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle = 0,hjust =1,size = 10),axis.text.y = element_text(size = 10))+ 
  NoLegend()+ ylab("") + theme_leukemia()

# Save the plot
ggsave("Fig_2D_Plot.pdf", width = 6.5, height = 5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

# Figure 2E
# SCENIC
scenicOptions <- readRDS("scenicOptions.Rds")
regulons <- loadInt(scenicOptions, "regulons")
cellInfo <- data.frame(seuratCluster=Idents(Stromal))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
aucMat <- getAUC(regulonAUC)
rownames(aucMat)
# Group by Seurat clusters and calculate mean activity per regulon
regulonActivity_byCluster <- sapply(
  split(rownames(cellInfo), cellInfo$seuratCluster),
  function(cells) rowMeans(aucMat[, cells, drop=FALSE])
)
regulonActivity_byCluster <- as.matrix(regulonActivity_byCluster)
# Top regulators
topRegulators <- reshape2::melt(regulonActivity_byCluster)
colnames(topRegulators) <- c("Regulon", "CellType", "RegulonActivity")
topRegulators <- topRegulators[which(topRegulators$RegulonActivity>0.01),]
viewTable(topRegulators)
topRegulonNames <- unique(topRegulators$Regulon)
topRegulonMatrix <- regulonActivity_byCluster[rownames(regulonActivity_byCluster) %in% topRegulonNames, ]

# Order by variance (most dynamic regulons first)
rowVar <- apply(topRegulonMatrix, 1, var)
topRegulonMatrix <- topRegulonMatrix[order(-rowVar), ]
topN <- 40
topRegulonMatrix <- head(topRegulonMatrix, topN)

library(viridis)
library(circlize)

col_fun <- colorRamp2(
  seq(0.001, 0.4, length.out = 5),
  viridis(5, option = "A", direction = -1)
)

library(Cairo)

CairoPDF("Fig_2E_Regulon_Heatmap.pdf", width = 3, height = 8)

Heatmap(topRegulonMatrix,
        name = "Regulon Activity",
        col = col_fun,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_rot = 75,        # Rotate cluster labels
        row_names_gp = gpar(fontsize = 6),  # Smaller row font
        column_names_gp = gpar(fontsize = 3),
        cluster_rows = F,
        cluster_columns = F)

dev.off()

#### Figure 2F
# Upload data
Stromal <- readRDS("")
Colors.Stromal = c("#EE6677","#CCBB44","#66CCEE","#228833","#4477AA","#AA3377")
DimPlot(Stromal,group.by = "AnnotationStromalCellArticle", reduction = "umap.rpca2", label = F, pt.size = 0.8, label.size = 8,cols = Colors.Stromal,
) + NoAxes() + theme_leukemia() + theme_no_axes()

# Save the plot
ggsave("Fig_2F_Plot.pdf", width = 6.5, height = 5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

#### Figure 2G
# Gene ontology

#----function to make a dotplot of gene sets found enriched in the marker genes------
plot.markergenes.GOs <- function(go.results,topGOs.number = 10,marker.genes.nr="ALL_up",custom.h=0,custom.w=0,provide.label=NULL,uniqueGO.filter=F){
  #only for overrepresented markergenes with positiv logFC
  #go.results-> calculated by get.markergenes.GOs function
  #topGOs.number-> how many top GOs should be included the final plot (ploting is more optimized for 10)
  #custom.h&w-> in case the perfect.width and height functions fail, this can be used to customize the graph dimensions for saving
  #provide.label-> in case x-axis label should be based on annotation provide a list of clusters+annotation e.g: "0"="Cap"
  #uniqueGO.filter -> in case the provided GOterms were filtered for unique ones per cluster
  GO.plot.width.heigt()
  if(uniqueGO.filter){
    per.cluster.unique.GOs<-NULL
    for (cluster in unique(go.results$cluster)) {
      all.GOs.found <- unique(go.results$term_id[go.results$cluster!=cluster])
      per.cluster.unique.GOs<-rbind(per.cluster.unique.GOs,go.results[go.results$cluster==cluster & !(go.results$term_id %in% all.GOs.found),])
    }#end of for
    go.results <- per.cluster.unique.GOs
    u="unique"
  }else{u=""}
  
  #----------------generate a dataframe of only the topXX GOs per cluster---------------
  go.results <- go.results[go.results$source %in% c("GO:BP","GO:MF","GO:CC","REAC","KEGG"),] #only take most relevant database
  go.results <- arrange(go.results,p_value)
  go.results <- arrange(go.results,cluster)
  
  #------extract and plot the top10 by p_value per cluster
  go.results.top10 <- go.results %>% group_by(cluster) %>% top_n(n = -topGOs.number,wt=p_value)
  go.results.top10<-arrange(go.results.top10,cluster)
  go.results.top10$cluster=factor(go.results.top10$cluster)
  order_desc_up <- unique(go.results.top10$term_name)
  go.results.top10.all <- go.results[go.results$term_name %in% order_desc_up,]
  x.angle=0
  x.hjust=0.5
  #in case labels are provided
  if (!(is.null(provide.label))) {
    go.results.top10.all$cluster=as.character(go.results.top10.all$cluster)
    for (cluster.nr in names(provide.label)) {go.results.top10.all[["cluster"]]<-ifelse(go.results.top10.all$cluster==cluster.nr,provide.label[[cluster.nr]],go.results.top10.all$cluster)}
    provide.label <- provide.label[unlist(provide.label)%in%go.results.top10.all$cluster]
    go.results.top10.all$cluster = factor(go.results.top10.all$cluster,levels = unlist(provide.label))
    x.angle=45
    x.hjust=1
  }#end of if
  print(
    ggplot(go.results.top10.all, aes(
      x = factor(term_name, level = order_desc_up),
      y = factor(cluster),
      size = precision,
      fill = -log10(p_value)
    )) +
      geom_point(shape = 21, color = "black", stroke = 0.2) +  # use shape 21 to apply fill
      scale_fill_viridis(option = "B", direction = -1, name = expression(-log[10](p))) +
      scale_size(range = c(3, 9)) +
      coord_flip() +
      ylab("Cluster") + xlab("") +
      ggtitle(paste0("Top ", topGOs.number, " GOs per cluster based on ", marker.genes.nr, " upregulated marker genes per cluster")) +
      theme(
        axis.text.x = element_text(angle = x.angle, hjust = x.hjust),
        legend.title = element_text(size = 8),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.key.size = unit(0.3, "cm")
      )
  )
} #end of plot.markergenes.GOs

plot.markergenes.GOs(go.results = GO.Stromal, topGOs.number = 6, provide.label = c("0"= "Early mesenchymal progenitors", "1"= "Adipogenic progenitors"),custom.h = 4,custom.w = 5.6)+
  theme_leukemia()+ theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
    axis.title.x = element_blank()
  )

# Save the plot
ggsave("Fig_2G_Plot.pdf", width = 4, height = 5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()



