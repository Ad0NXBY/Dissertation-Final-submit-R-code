#DEG============================================================================
##Load necessary libraries------------------------------------------------------
library(DESeq2)
library(limma)
library(dplyr)

##Reading in the normalized counts file----------------------------------------
experimental_data <- read.csv("C:/Users/Brandon/Documents/MRes Big Data Biology/R studio/List_for_analysis with combined code//normalized_counts.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(experimental_data) <- experimental_data$X
experimental_data <- experimental_data[,-1]

##Filter low-expressed genes using limma's filterByExpr()----------------------
keep_genes <- filterByExpr(experimental_data, 
                           group = condition, 
                           min.count = 15, 
                           min.total.count = 70, 
                           min.prop = 0.2)

##Subset the data (correct way)-------------------------------------------------
experimental_data_filtered <- experimental_data[keep_genes, ]

condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23", 3)))


##Creating dds object-----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(experimental_data),
  colData = data.frame(condition = condition),
  design = ~ condition
)

dds <- DESeq(dds)

##Creating .csv DEG list for FIHKO22--------------------------------------------
results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti"))
results_plenti_vs_KO22_df <- as.data.frame(results_plenti_vs_KO22)
results_plenti_vs_KO22_df$gene <- rownames(results_plenti_vs_KO22_df)

#Save full DEG results
write.csv(results_plenti_vs_KO22_df, 
          file.path(data_dir, "DEG_KO22_vs_Plenti_filterByExpr.csv"), 
          row.names = TRUE)

##Filter for significant DEGs for KO22 GSEA (padj < 0.05)-----------------------
GSEA_FIH22KO <- results_plenti_vs_KO22_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(log2FoldChange))

#Save the filtered DEG list
write.csv(GSEA_FIH22KO,
          file.path(data_dir, "GSEA_plenti_vs_KO22.csv"),
          row.names = TRUE)

##Prepare ranked list for GSEA--------------------------------------------------
GSEA_ranked_list_KO22 <- GSEA_FIH22KO %>%
  dplyr::select(gene, log2FoldChange)

#Write .rnk file (no headers or row names as required by GSEA)
write.table(
  GSEA_ranked_list_KO22,
  file.path(data_dir, "GSEA_ranked_list_KO22_filter_by_padj_lessthan0.05.rnk"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)



##Creating .csv DEG list for FIHKO23--------------------------------------------
results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti"))
results_plenti_vs_KO23_df <- as.data.frame(results_plenti_vs_KO23)
results_plenti_vs_KO23_df$gene <- rownames(results_plenti_vs_KO23_df)

#Save full DEG results
write.csv(results_plenti_vs_KO23_df, 
          file.path(data_dir, "DEG_KO23_vs_Plenti_filterByExpr.csv"), 
          row.names = TRUE)

##Filter for significant DEGs for KO23 GSEA(padj < 0.05)------------------------
GSEA_FIH23KO <- results_plenti_vs_KO23_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(log2FoldChange))

#Save the filtered DEG list
write.csv(GSEA_FIH23KO,
          file.path(data_dir, "GSEA_plenti_vs_KO23.csv"),
          row.names = TRUE)

##Prepare ranked list for GSEA--------------------------------------------------
GSEA_ranked_list_KO23 <- GSEA_FIH23KO %>%
  dplyr::select(gene, log2FoldChange)

#Write .rnk file (no headers or row names as required by GSEA)
write.table(
  GSEA_ranked_list_KO23,
  file.path(data_dir, "GSEA_ranked_list_KO23_filter_by_padj_lessthan0.05.rnk"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

#Scree & PCA Plot===========================================================================
#Load necessary libraries------------------------------------------------------
library(DESeq2)
library(factoextra)
library(ggpubr)
library(ggplot2)

##Transform the data for PCA----------------------------------------------------
vsd <- vst(dds, blind = FALSE) #variance stabilizing transformation

condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23",3)))
##Scree plot -------------------------------------------------------------------
DS1.svd <- assay(vsd) |> 
  t() |> 
  prcomp(scale = FALSE) # PCA using prcomp()
summary(DS1.svd)

pScree <- fviz_eig(DS1.svd, addlabels = TRUE) + 
  theme_pubr(base_size = 20)

##PCA Plot with colors----------------------------------------------------------
pPCA <- fviz_pca_ind(DS1.svd, 
                     label = "all",  # Ensure all labels are visible
                     habillage = condition,  # Color by condition
                     repel = TRUE, # Prevent label overlap
                     mean.point = FALSE,  #remove centroid marker
                     pointsize = 4, 
                     geom.ind = "point",
                     labelsize = 6
)+  
  labs(title = "PCA Plot",
       x = "PC1",
       y = "PC2") + 
  theme(
    plot.title = element_text(size = 20),  # Title size
    axis.title = element_text(size = 18),                 # Axis titles
    axis.text = element_text(size = 14),                  # Axis tick labels
    legend.title = element_text(size = 16),               # Legend title
    legend.text = element_text(size = 14)                 # Legend text
  )

#Arrange plots
pScreePCA <- ggarrange(pScree, pPCA,
                       labels = c("A", "B"),
                       ncol = 1, nrow = 2)

print(pScreePCA)
ggsave(file.path(plot_dir, "PCA_ScreePlot.png"), plot = pScreePCA, width = 12, height = 20)




#VENN PLOT=================================================================================
##Load necessary libraries------------------------------------------------------
library(VennDiagram)
library(grid)
library(ggplot2)

##Filtering out genes and creating a gene list----------------------------------
up_genes_KO22 <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange > 0, ])
down_genes_KO22 <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange < 0, ])

up_genes_KO23 <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange > 0, ])
down_genes_KO23 <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange < 0, ])

gene_lists <- list(
  "FIH22KO_up" = up_genes_KO22,
  "FIH22KO_down" = down_genes_KO22,
  "FIH23KO_up" = up_genes_KO23,
  "FIH23KO_down" = down_genes_KO23
)

##Creating Venn plot------------------------------------------------------------
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("KO22 Upregulated", "KO22 Downregulated", "KO23 Upregulated", "KO23 Downregulated"),
  filename = NULL,
  imagetype = "jpg",  
  height = 750,
  width = 750,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = "transparent",
  fill = c("#F8766D", "#00BA38", "#619CFF", "purple4"),
  alpha = 0.50,
  label.col = "darkblue",
  cex = 2,
  fontfamily = "Arial",  
  fontface = "bold",
  cat.col = c("darkred", "darkgreen", "darkblue", "purple4"),
  cat.cex = 1.3,
  cat.fontfamily = "Arial", 
  cat.default.pos = "outer",
  cat.pos = c(-40, 40, -30, 30),
  cat.dist = c(0.3, 0.3, 0.2, 0.2),
  cat.fontface = "bold",
  rotation.degree = 0,
  margin = 0.2
)

grid.draw(venn.plot)

ggsave(file.path(plot_dir, "vennplot.jpg"), plot = venn.plot, width = 12, height = 8, dpi = 800)


#VOLCANO PLOT===================================================================
##Load necessary libraries------------------------------------------------------
library(ggplot2)
library(ggrepel)

##Volcano plot of Plent vs KO22-------------------------------------------------
#Ensure gene names are set as rownames
volcano_data_22 <- as.data.frame(results_plenti_vs_KO22)
volcano_data_22$Gene <- rownames(volcano_data_22)

#Remove NA rows
volcano_data_22 <- na.omit(volcano_data_22)

#Define significance thresholds
volcano_data_22$Significance <- "Not Significant"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange < -log2(2)] <- "Downregulated"

#Separate top 25 upregulated and downregulated genes
top_upregulated_22 <- volcano_data_22 %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(25)

top_downregulated_22 <- volcano_data_22 %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(25)

#Combine top 30 genes
top_genes_22 <- rbind(top_upregulated_22, top_downregulated_22)


#Create volcano plot with geom_text_repel
Volcano_Plot_Plenti_v_KO22 <- ggplot(volcano_data_22, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  geom_vline(xintercept = c(-log2(2), 0, log2(2)), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    data = subset(volcano_data_22, padj < 0.05 & abs(log2FoldChange) > log2(2)),
    aes(label = Gene),
    size = 3,
    color = "black",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = 25
  ) +
  xlim(-10, 10) +
  labs(title = "Volcano Plot (Plenti vs KO22)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),  # bigger title
    axis.title = element_text(size = 16),  # bigger x and y axis labels
    axis.text = element_text(size = 14)    # bigger tick labels
  )


#Save plot
ggsave(file.path(plot_dir, "Volcano_Plot_Plenti_vs_KO22.jpg"), plot = Volcano_Plot_Plenti_v_KO22,
       width = 10, height = 6, dpi = 800)
##Volcano plot of Plent vs KO23-------------------------------------------------
#Ensure gene names are set as rownames
volcano_data_23 <- as.data.frame(results_plenti_vs_KO23)
volcano_data_23$Gene <- rownames(volcano_data_23)

#Remove NA rows
volcano_data_23 <- na.omit(volcano_data_23)

#Define significance thresholds
volcano_data_23$Significance <- "Not Significant"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange < -log2(2)] <- "Downregulated"

#Separate top 25 upregulated and downregulated genes
top_upregulated_23 <- volcano_data_23 %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(25)

top_downregulated_23 <- volcano_data_23 %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(25)

#Combine top 30 genes
top_genes_23 <- rbind(top_upregulated_23, top_downregulated_23)


#Create volcano plot with geom_text_repel
Volcano_Plot_Plenti_v_KO23 <- ggplot(volcano_data_23, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  geom_vline(xintercept = c(-log2(2), 0, log2(2)), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    data = subset(volcano_data_23, padj < 0.05 & abs(log2FoldChange) > log2(2)),
    aes(label = Gene),
    size = 3,
    color = "black",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = 25
  ) +
  xlim(-10, 10) +
  labs(title = "Volcano Plot (Plenti vs KO23)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),  # bigger title
    axis.title = element_text(size = 16),  # bigger x and y axis labels
    axis.text = element_text(size = 14)    # bigger tick labels
  )


#Save plot
ggsave(file.path(plot_dir, "Volcano_Plot_Plenti_vs_KO23.jpg"), plot = Volcano_Plot_Plenti_v_KO23,
       width = 10, height = 6, dpi = 800)

# Combine KO22 and KO23 volcano plots vertically
Volcano_Combined <- ggarrange(
  Volcano_Plot_Plenti_v_KO22, 
  Volcano_Plot_Plenti_v_KO23,
  labels = c("A", "B"),        # Panel labels
  ncol = 1, nrow = 2,          # Vertical layout (A above B)
  common.legend = TRUE,        # One shared legend
  legend = "right"             # Shared legend on the right
)

# Save combined figure
ggsave(file.path(plot_dir, "Volcano_Plots_Combined.jpg"), plot = Volcano_Combined,
       width = 10, height = 12, dpi = 800)