## GSVA ========================================================================
## Load necessary libraries ----------------------------------------------------
library(GSVA)
library(msigdbr)
library(GSEABase)
library(ComplexHeatmap)
library(ggplot2)
library(grid)   # Needed for gpar()
library(gridExtra)  # Optional, in case needed for future multi-panel plots

## Get human Hallmark gene sets -----------------------------------------------
msigdbr_species <- msigdbr(species = "Homo sapiens")
hallmark_genesets_df <- msigdbr_species[msigdbr_species$gs_cat == "H", ]

## Create GeneSetCollection ----------------------------------------------------
gset.idx.list <- split(hallmark_genesets_df$gene_symbol, hallmark_genesets_df$gs_name)

geneSets <- lapply(names(gset.idx.list), function(name) {
  GeneSet(geneIds = unique(gset.idx.list[[name]]), 
          setName = name, 
          geneIdType = SymbolIdentifier())
})
geneSetCollection <- GeneSetCollection(geneSets)

expr_matrix <- assay(vsd)

## Run GSVA --------------------------------------------------------------------
param <- gsvaParam(expr = expr_matrix, 
                   geneSets = geneSetCollection, 
                   kcdf = "Poisson")
gsva_results <- gsva(param)

write.csv(gsva_results, file.path(data_dir, "GSVA_results.csv"), row.names = TRUE)

## T-tests: Control vs KO22 ---------------------------------------------------
ttest_results_KO22 <- t(apply(gsva_results, 1, function(row) {
  test <- t.test(row[c("FIH22KO1_Count", "FIH22KO2_Count", "FIH22KO3_Count")],
                 row[c("Control1_Count", "Control2_Count", "Control3_Count")])
  c(t_value = test$statistic, p_value = test$p.value)
}))
ttest_results_KO22 <- as.data.frame(ttest_results_KO22)
colnames(ttest_results_KO22) <- c("t_value", "p_value")

## T-tests: Control vs KO23 ---------------------------------------------------
ttest_results_KO23 <- t(apply(gsva_results, 1, function(row) {
  test <- t.test(row[c("FIH23KO1_Count", "FIH23KO2_Count", "FIH23KO3_Count")],
                 row[c("Control1_Count", "Control2_Count", "Control3_Count")])
  c(t_value = test$statistic, p_value = test$p.value)
}))
ttest_results_KO23 <- as.data.frame(ttest_results_KO23)
colnames(ttest_results_KO23) <- c("t_value", "p_value")


## Bar plot: KO22 --------------------------------------------------------------
bar_data_KO22 <- data.frame(
  GeneSet = sub("HALLMARK_", "", rownames(ttest_results_KO22)),
  T_Value = ttest_results_KO22$t_value,
  P_Value = ttest_results_KO22$p_value,
  
  FillColor = ifelse(ttest_results_KO22$p_value < 0.05,
                     ifelse(ttest_results_KO22$t_value > 0, "red", "blue"),
                     "grey")
)
bar_data_KO22 <- bar_data_KO22[order(bar_data_KO22$T_Value, decreasing = TRUE), ]

GSVA_BarPlot_KO22 <- ggplot(bar_data_KO22, aes(x = reorder(GeneSet, T_Value), y = T_Value, fill = FillColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  coord_flip() +
  labs(x = "Gene Set", y = "t value of GSVA", title = "GSVA: FIH KO22 vs Control") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "GSVA_BarPlot_Plenti_vs_KO22.png"),
       plot = GSVA_BarPlot_KO22, width = 10, height = 6, dpi = 800)

## Bar plot: KO23 --------------------------------------------------------------
bar_data_KO23 <- data.frame(
  GeneSet = sub("HALLMARK_", "", rownames(ttest_results_KO23)),
  T_Value = ttest_results_KO23$t_value,
  P_Value = ttest_results_KO23$p_value,
  FillColor = ifelse(ttest_results_KO23$p_value < 0.05,
                     ifelse(ttest_results_KO23$t_value > 0, "red", "blue"),
                     "grey")
)
bar_data_KO23 <- bar_data_KO23[order(bar_data_KO23$T_Value, decreasing = TRUE), ]

GSVA_BarPlot_KO23 <- ggplot(bar_data_KO23, aes(x = reorder(GeneSet, T_Value), y = T_Value, fill = FillColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  coord_flip() +
  labs(x = "Gene Set", y = "t value of GSVA", title = "GSVA: FIH KO23 vs Control") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

ggsave(file.path(plot_dir, "GSVA_BarPlot_Plenti_vs_KO23.png"),
       plot = GSVA_BarPlot_KO23, width = 10, height = 6, dpi = 800)


# Combine KO22 and KO23 GSVA bar plots vertically
GSVA_Combined <- ggarrange(
  GSVA_BarPlot_KO22,
  GSVA_BarPlot_KO23,
  labels = c("A", "B"),
  ncol = 1, nrow = 2,
  common.legend = TRUE,
  legend = "right"
)

# Save combined GSVA figure
ggsave(
  file.path(plot_dir, "GSVA_BarPlots_KO22_KO23_Combined.jpg"),
  plot = GSVA_Combined,
  width = 14, height = 20, dpi = 800
)




##Extract relevant GSVA scores for all 3 conditions for heat map----------------
gsva_subset <- gsva_results[, c("Control1_Count", "Control2_Count", "Control3_Count",
                                "FIH22KO1_Count", "FIH22KO2_Count", "FIH22KO3_Count",
                                "FIH23KO1_Count", "FIH23KO2_Count", "FIH23KO3_Count")]

rownames(gsva_subset) <- sub("HALLMARK_", "", rownames(gsva_subset))
gsva_scaled <- t(scale(t(gsva_subset)))  # mean-center

## Annotation: column (conditions) ---------------------------------------------
annotation_col <- data.frame(
  Condition = c(rep("Control", 3), rep("KO22", 3), rep("KO23", 3))
)
rownames(annotation_col) <- colnames(gsva_scaled)

ha <- HeatmapAnnotation(
  Condition = annotation_col$Condition,
  col = list(Condition = c("Control" = "#4DBBD5", "KO22" = "#E64B35", "KO23" = "#00A087"))
)

## Annotation: row (significance status) ---------------------------------------
sig_table <- data.frame(
  KO22 = factor(
    ifelse(ttest_results_KO22$p_value < 0.05,
           ifelse(ttest_results_KO22$t_value > 0, "Up", "Down"), "NS"),
    levels = c("Up", "Down", "NS")
  ),
  KO23 = factor(
    ifelse(ttest_results_KO23$p_value < 0.05,
           ifelse(ttest_results_KO23$t_value > 0, "Up", "Down"), "NS"),
    levels = c("Up", "Down", "NS")
  ),
  row.names = sub("HALLMARK_", "", rownames(ttest_results_KO22))
)
sig_table <- sig_table[rownames(gsva_scaled), ]

sig_cols <- list(
  KO22 = c(Up = "firebrick3", Down = "royalblue", NS = "grey90"),
  KO23 = c(Up = "firebrick3", Down = "royalblue", NS = "grey90")
)

row_sig_ha <- rowAnnotation(
  df = sig_table,
  col = sig_cols,
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

## Heatmap: Create & Save ------------------------------------------------------
GSVA_heatmap <- Heatmap(
  gsva_scaled,
  name               = "GSVA\nscore",
  top_annotation     = ha,
  left_annotation    = row_sig_ha,
  cluster_rows       = TRUE,
  cluster_columns    = FALSE,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  row_names_gp       = gpar(fontsize = 6),
  column_names_gp    = gpar(fontsize = 10, fontface = "bold"),
  col                = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  heatmap_legend_param = list(title = "GSVA score")
)

jpeg(file.path(plot_dir,"GSVA_heatmap_with_sig.jpg"), width = 1200, height = 1000, res = 150)
draw(GSVA_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()