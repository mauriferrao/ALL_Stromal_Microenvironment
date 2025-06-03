

### Code to generate figures in Figure 8
# Cell Deconvolution Analysis
# Article. Leukemia Journal
# Upload packages
library(readxl)  # For reading Excel files
library(dplyr)   # For data manipulation
library(stringr) # For string manipulation
library(writexl) # For writing Excel files
library(ggplot2)
library(corrplot)
library("PerformanceAnalytics")
library(RColorBrewer)
library(scater)
library(pals)
library(knitr)
library(kableExtra)
library(emmeans)
library(gt)
library(scales)
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

# Figure 8
# Read the Excel file
metadata <- read_excel(path = '')

#### Figure 8A. 
# Table of metadata
# Create summary
summary_table <- metadata %>%
  summarise(
    `Average Age` = round(mean(Age, na.rm = TRUE), 1),
    `Male (n)` = sum(Sex == "Male", na.rm = TRUE),
    `Female (n)` = sum(Sex == "Female", na.rm = TRUE),
    `Not Available (Sex)` = sum(Sex == "Not Available", na.rm = TRUE),
    `Diagnosis (n)` = sum(sample_type == "Diagnosis", na.rm = TRUE),
    `Relapse (n)` = sum(sample_type == "Relapse", na.rm = TRUE)
  ) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Metric") %>%
  rename(Value = V1)

summary_table_gt <- summary_table %>%
  gt() %>%
  tab_header(
    title = "Sample Metadata Summary"
  )

library(webshot2)

gtsave(summary_table_gt, filename = "Sample_Metadata_Summary.pdf")

# Figure 8B.
# Pie chart of molecular subtypes

# Prepare data
subtype_counts <- metadata %>%
  dplyr::count(molecular_subtype) %>%
  arrange(desc(n)) %>%
  mutate(
    fraction = n / sum(n),
    percent_label = paste0(molecular_subtype, " (", round(fraction * 100, 1), "%)"),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymin + ymax) / 2,
    label = paste0(molecular_subtype, "\n(", n, ")")
  )

# Generate pastel color palette (Pastel1 allows up to 9, we interpolate to 25)
pastel_colors <- colorRampPalette(brewer.pal(9, "Set3"))(nrow(subtype_counts))

# Plot
ggplot(subtype_counts, aes(ymin = ymin, ymax = ymax, xmin = 0.3, xmax = 1, fill = percent_label)) +
  geom_rect() +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(
    aes(y = label_pos, x = 1.2, label = label),
    size = 3,
    hjust = 0
  ) +
  geom_segment(
    aes(y = label_pos, yend = label_pos, x = 1, xend = 1.18),
    color = "gray30",
    size = 0.3
  ) +
  scale_fill_manual(values = pastel_colors) +
  labs(
    title = "Molecular Subtype Distribution",
    fill = "Subtype (% of total)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) + theme_leukemia() + theme_no_axes()

# Save the plot
ggsave("Fig_8B_PieChart.pdf", width = 12, height = 8, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()


# Correlation: Adipogenic Progenitors vs Tumor Cell Percentage (if available)
# Clean and convert Tumor_cell_percentage to numeric (remove % symbols if present)
metadata$Tumor_cell_percentage <- as.numeric(gsub("%", "", metadata$Tumor_cell_percentage))
# Ensure Adipogenic progenitors is numeric
metadata$Adipogenic_Progenitors <- as.numeric(metadata$`Adipogenic progenitors`)
# Spearman correlation
cor_tumor <- cor.test(metadata$Tumor_cell_percentage, metadata$Adipogenic_Progenitors, 
                      method = "spearman", use = "complete.obs")
print(cor_tumor)

sum(!is.na(metadata$Tumor_cell_percentage))
ggscatter(metadata, x = "Tumor_cell_percentage", y = "Adipogenic progenitors", 
          add = "reg.line", add.params = list(color= "#E5D8BD"),
          conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          color = "#DECBE4",
          xlab = "% Leukemic cells", ylab = "Adipogenic progenitors",
          ggtheme = theme_classic()
) + theme_leukemia()

# Save the plot
ggsave("Fig_Correlation_AdipogenicProgenitors_LeukemicCellPercentage_Plot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()


#### Correlation: Adipogenic Progenitors vs Age. (Supplementary)
ggplot(metadata, aes(x = Age, y = `Adipogenic progenitors`)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Adipogenic Progenitors vs Age")

cor_age <- cor.test(metadata$Age, metadata$`Adipogenic progenitors`, method = "spearman")
print(cor_age)

ggscatter(metadata, x = "Age", y = "Adipogenic progenitors", 
          add = "reg.line", add.params = list(color= "#E5D8BD"),
          conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          color = "#DECBE4",
          xlab = "Age", ylab = "Adipogenic progenitors",
          ggtheme = theme_classic()
)+ theme_leukemia()

# Save the plot
ggsave("Fig_Supplementary_Correlation_AdipogenicProgenitors_Age_Plot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

### Comparison with sex
# Check the unique values in the 'Sex' column
table(metadata$Sex)
# Remove rows with invalid values in the Sex column
metadata_clean <- metadata %>%
  filter(Sex %in% c("Male", "Female"))

# Re-run the Wilcoxon test
stat_test <- wilcox.test(`Adipogenic progenitors` ~ Sex, data = metadata_clean)
print(stat_test)

# Create box plot and add statistical annotations
ggplot(metadata_clean, aes(x = Sex, y = `Adipogenic progenitors`, fill = Sex)) +
  geom_boxplot(outlier.size = 0.1) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +  # Add p-value from Wilcoxon test
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Adipogenic Progenitors") +
  coord_cartesian(ylim = )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_leukemia()+
  theme(legend.position = "none")

# Save the plot
ggsave("Fig_Supplementary_AdipogenicProgenitors_Sex_Plot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()


### Difference in Adipogenic progenitors between samples at Diagnosis and Relapse
# while accounting for the covariates Age, Sex, Molecular Subtype, and Tumor Cell Percentage

# Ensure proper factor levels
metadata$sample_type <- factor(metadata$sample_type, levels = c("Diagnosis", "Relapse"))
metadata$Sex <- factor(metadata$Sex, levels = c("Male", "Female", "Not Available"))
metadata$molecular_subtype <- factor(metadata$molecular_subtype)

# Clean Age (blank to NA)
metadata$Age <- ifelse(metadata$Age == "", NA, metadata$Age)
metadata$Age <- as.numeric(metadata$Age)

# Run the linear model
model <- lm(`Adipogenic progenitors` ~ sample_type + Age + Sex + molecular_subtype, 
            data = metadata, na.action = na.omit)

summary(model)

# Calculate estimated marginal means for sample_type (Diagnosis vs Relapse)
emm <- emmeans(model, specs = "sample_type")
print(emm)

# Get pairwise comparison (p-value)
contrast_result <- contrast(emm, method = "pairwise")
print(contrast_result)

# Convert emmeans to dataframe for ggplot
emm_df <- as.data.frame(emm)

# Extract p-value for plot annotation
pval <- summary(contrast_result)$p.value

wilcox_diag_relapse <- wilcox.test(`Adipogenic progenitors` ~ sample_type, data = metadata)
print(wilcox_diag_relapse)

library(ggsignif)

## Plot
ggplot(metadata, aes(x = sample_type, y = `Adipogenic progenitors`, fill = sample_type)) +
  geom_boxplot(outlier.size = 0.1) +
  labs(title = "", x = "") +
  coord_cartesian(ylim = c(0, 10)) +
  theme_leukemia() +
  theme(
    legend.position = "none"
  ) +
  geom_signif(
    comparisons = list(c("Diagnosis", "Relapse")),
    map_signif_level = FALSE,
    annotations = "p < 0.0001",
    y_position = 8,   # Adjust this based on your real max values
    tip_length = 0.01,
    textsize = 5
  ) +
  scale_fill_manual(values = c("Diagnosis" = "#DEC197", "Relapse" = "#B3C76B")) 


# Save the plot
ggsave("Fig_8_AdipogenicProgenitors_Diagnosis_Relapse_Plot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

#### Determine adipogenic progenitors between molecular subtypes (exclude relapse samples)
metadata_diagnosis <- metadata %>%
  filter(sample_type == "Diagnosis")

nrow(metadata_diagnosis)
table(metadata$sample_type)


#### ✅ 1. Correlation: Adipogenic Progenitors vs Molecular Subtype
# Count samples per molecular subtype
subtype_counts <- metadata_diagnosis %>%
  count(molecular_subtype) %>%
  mutate(label = paste0(molecular_subtype, " (", n, ")"))  # Add count in parentheses

# Merge the new labels into the metadata
metadata_diagnosis <- metadata_diagnosis %>%
  left_join(subtype_counts, by = "molecular_subtype")

# Create the pastel color palette
pastel_colors <- colorRampPalette(brewer.pal(9, "Set2"))(length(unique(metadata_diagnosis$molecular_subtype)))

# Plot with custom pastel colors
ggplot(metadata_diagnosis, aes(x = label, y = `Adipogenic progenitors`, fill = molecular_subtype)) +
  geom_boxplot(outlier.size = 0.1) +
  labs(title = "", x = "") +
  coord_cartesian(ylim = c(0, 10)) +
  theme_leukemia() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = pastel_colors)


kruskal_molecular <- kruskal.test(`Adipogenic progenitors` ~ molecular_subtype, data = metadata_diagnosis)
print(kruskal_molecular)

# Output: Kruskal-Wallis result shows a highly significant difference in Adipogenic Progenitors across molecular subtypes (p = 4.6e-12)
# Identify which molecular subtypes differ significantly

library(FSA)
dunn_molecular <- dunnTest(`Adipogenic progenitors` ~ molecular_subtype, data = metadata_diagnosis, method = "bonferroni")
print(dunn_molecular)

# Convert to dataframe for easy filtering
dunn_df <- as.data.frame(dunn_molecular$res)

# Filter significant comparisons (Bonferroni adjusted p < 0.05)
sig_dunn <- dunn_df %>% filter(P.adj < 0.05)
print(sig_dunn)

# Save the plot
ggsave("Fig_8_AdipogenicProgenitors_MolecularSubtype_Plot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()

### Comparison between risk stratification
##### 3) Genomic Subtype Risk Group assignments. As discussed with Josef and related to article Florent Malard, Mohamad Mohty. Lancet 2020; 395: 1146–62
# High_risk: BCR-ABL, BCR-ABL-like, KMT2Ar, Hypodiploid, Near-haploid, MEF2Dr, MYC, iAMP21, TCF3-HLF, ABL
# Not_High_risk: ETV6-RUNX1, ETV6-RUNX1-like, Hyperdiploid, TCF3-PBX1, DUX4r, CRLF2r, CEBPr, KMT2Ar-like, ZNF384r-like, PAX5, ZNF384r, IGF2BP1r, IKFZ1, Other


metadata <- metadata %>%
  mutate(Risk_Stratification = case_when(
    molecular_subtype %in% c("ETV6-RUNX1", "Hyperdiploid", "TCF3-PBX1", "DUX4r", 
                             "PAX5", "ZNF384r","ETV6-RUNX1-like", "CRLF2r", "CEBPr", "KMT2Ar-like",
                             "ZNF384r-like", "IGF2BP1r", "TCF3-FLI1","NUTM1","PAX5P80R","IKZF1","NOS") ~ "Not_High_risk",
    molecular_subtype %in% c("BCR-ABL", "BCR-ABL-like", "KMT2Ar", "Hypodiploid", 
                             "Near-haploid", "MEF2Dr", "MYC", "iAMP21", 
                             "TCF3-HLF", "IKFZ1", "ABL") ~ "High_risk",
    TRUE ~ "Unclassified"  # Default fallback if none of the above matches
  ))

table(metadata$Risk_Stratification)


# Adipogenic progenitors between risk groups
# while accounting for the covariates Age, Sex, Molecular Subtype

# Ensure proper factor levels
metadata$Risk_Stratification <- factor(metadata$Risk_Stratification, levels = c("Not_High_risk", "High_risk"))
metadata$Sex <- factor(metadata$Sex, levels = c("Male", "Female", "Not Available"))
metadata$molecular_subtype <- factor(metadata$molecular_subtype)

# Clean Age (blank to NA)
metadata$Age <- ifelse(metadata$Age == "", NA, metadata$Age)
metadata$Age <- as.numeric(metadata$Age)

# Run the linear model
model <- lm(`Adipogenic progenitors` ~ Risk_Stratification + Age + Sex, 
            data = metadata, na.action = na.omit)

summary(model)

# Calculate estimated marginal means for sample_type
emm <- emmeans(model, specs = "Risk_Stratification")
print(emm)

# Get pairwise comparison (p-value)
contrast_result <- contrast(emm, method = "pairwise")
print(contrast_result)

# Convert emmeans to dataframe for ggplot
emm_df <- as.data.frame(emm)

# Extract p-value for plot annotation
pval <- summary(contrast_result)$p.value

library(ggsignif)

## Plot
ggplot(metadata, aes(x = Risk_Stratification, y = `Adipogenic progenitors`, fill = Risk_Stratification)) +
  geom_boxplot(outlier.size = 0.1) +
  labs(title = "", x = "") +
  coord_cartesian(ylim = c(0, 10)) +
  theme_leukemia() +
  theme(
    legend.position = "none"
  ) +
  geom_signif(
    comparisons = list(c("Not_High_risk", "High_risk")),
    map_signif_level = FALSE,
    annotations = "p = 0.074",
    y_position = 8,   # Adjust this based on your real max values
    tip_length = 0.01,
    textsize = 5
  ) +
  scale_fill_manual(values = c("Not_High_risk" = "#DECBE4", "High_risk" = "#A89CB0"))


# Save the plot
ggsave("Fig_8_AdipogenicProgenitors_RiskStratification_Plot.pdf", width = 6, height = 5.5, units = "in", dpi = 300, bg = "transparent", device = cairo_pdf)
dev.off()




