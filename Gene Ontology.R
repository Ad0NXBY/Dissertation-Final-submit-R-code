#GENE ONTOLOGY=================================================================
##Load necessary libraries------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(forcats)
library(ggplot2)

##FIH22KO --------------------------------------------------------------------
up_genes_FIH22KO <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange > 0, ])
down_genes_FIH22KO <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange < 0, ])
background_genes <- rownames(experimental_data)

###Processing of GO analysis results of upregulated gene expression-------------
ego_up_KO22 <- enrichGO(gene         = up_genes_FIH22KO,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = "SYMBOL",
                        ont          = "ALL",  
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

ego_up_df_KO22 <- as.data.frame(ego_up_KO22)

sorted_up_bp_KO22 <- ego_up_df_KO22 %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_up_cc_KO22 <- ego_up_df_KO22 %>%
  filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_up_mf_KO22 <- ego_up_df_KO22 %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_ego_up_df_KO22 <- bind_rows(sorted_up_bp_KO22, sorted_up_cc_KO22, sorted_up_mf_KO22)


#Sort within each ONTOLOGY group by Count
sorted_ego_up_df_KO22 <- sorted_ego_up_df_KO22 %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


###Creating upregulated Gene ontology graph-----------------------------------
Gene_ontology_KO22_UP <- ggplot(sorted_ego_up_df_KO22, aes(x = Count, y = Description)) +
  geom_point(aes(color = p.adjust), size = 8) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", nrow = 3, strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "GO Terms Up-regulated Genes KO22",
    x = "Gene Count",
    y = "GO Term",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Up-regulated Genes KO22.jpg"), plot = Gene_ontology_KO22_UP, width = 12, height = 8, dpi = 800)


###Processing of GO analysis results of downregulated gene expression-------------
ego_down_KO22 <- enrichGO(gene         = down_genes_FIH22KO,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

ego_down_df_KO22 <- as.data.frame(ego_down_KO22)

sorted_down_bp_KO22 <- ego_down_df_KO22 %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_down_cc_KO22 <- ego_down_df_KO22 %>%
  filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_down_mf_KO22 <- ego_down_df_KO22 %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_ego_down_df_KO22 <- bind_rows(sorted_down_bp_KO22, sorted_down_cc_KO22, sorted_down_mf_KO22)

#Sort within each ONTOLOGY group by Count
sorted_ego_down_df_KO22 <- sorted_ego_down_df_KO22 %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


###Creating downregulated Gene ontology graph-----------------------------------
Gene_ontology_KO22_DOWN <- ggplot(sorted_ego_down_df_KO22, aes(x = Count, y = Description)) +
  geom_point(aes(color = p.adjust), size = 8) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", nrow = 3, strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "GO Terms for Down-regulated Genes KO22",
    x = "Gene Count",
    y = "GO Term",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Down-regulated Genes KO22.jpg"), plot = Gene_ontology_KO22_DOWN, width = 12, height = 8, dpi = 800)

# Combine Upregulated (A) and Downregulated (B) GO plots vertically
GO_Combined_KO22 <- ggarrange(
  Gene_ontology_KO22_UP, 
  Gene_ontology_KO22_DOWN,
  labels = c("A", "B"),   # Figure A = Upregulated, Figure B = Downregulated
  ncol = 1, nrow = 2,     # Vertical layout
  common.legend = TRUE,   # Shared legend
  legend = "right"        # Legend on the right
)

# Save combined figure
ggsave(file.path(plot_dir, "GO_Terms_KO22_Combined.jpg"), plot = GO_Combined_KO22,
       width = 12, height = 14, dpi = 800)


##FIH23KO------------------------------------------------------------------------
up_genes_FIH23KO <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange > 0, ])
down_genes_FIH23KO <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange < 0, ])
background_genes <- rownames(experimental_data)

###Processing of GO analysis results of upregulated gene expression-------------
ego_up_KO23 <- enrichGO(gene         = up_genes_FIH23KO,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = "SYMBOL",
                        ont          = "ALL",  
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

ego_up_df_KO23 <- as.data.frame(ego_up_KO23)

sorted_up_bp_KO23 <- ego_up_df_KO23 %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_up_cc_KO23 <- ego_up_df_KO23 %>%
  filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_up_mf_KO23 <- ego_up_df_KO23 %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_ego_up_df_KO23 <- bind_rows(sorted_up_bp_KO23, sorted_up_cc_KO23, sorted_up_mf_KO23)

#Sort within each ONTOLOGY group by Count
sorted_ego_up_df_KO23 <- sorted_ego_up_df_KO23 %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

###Creating upregulated Gene ontology graph-----------------------------------
Gene_ontology_KO23_UP <- ggplot(sorted_ego_up_df_KO23, aes(x = Count, y = Description)) +
  geom_point(aes(color = p.adjust), size = 8) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", nrow = 3, strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "GO Terms Up-regulated Genes KO23",
    x = "Gene Count",
    y = "GO Term",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Up-regulated Genes KO23.jpg"), plot = Gene_ontology_KO23_UP, width = 12, height = 8, dpi = 800)

###Processing of GO analysis results of downregulated gene expression-------------
ego_down_KO23 <- enrichGO(gene         = down_genes_FIH23KO,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

ego_down_df_KO23 <- as.data.frame(ego_down_KO23)

sorted_down_bp_KO23 <- ego_down_df_KO23 %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_down_cc_KO23 <- ego_down_df_KO23 %>%
  filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_down_mf_KO23 <- ego_down_df_KO23 %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_ego_down_df_KO23 <- bind_rows(sorted_down_bp_KO23, sorted_down_cc_KO23, sorted_down_mf_KO23)

#Sort within each ONTOLOGY group by Count
sorted_ego_down_df_KO23 <- sorted_ego_down_df_KO23 %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


###Creating downregulated Gene ontology graph-----------------------------------
Gene_ontology_KO23_DOWN <- ggplot(sorted_ego_down_df_KO23, aes(x = Count, y = Description)) +
  geom_point(aes(color = p.adjust), size = 8) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", nrow = 3, strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "GO Terms for Down-regulated Genes KO23",
    x = "Gene Count",
    y = "GO Term",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Down-regulated Genes KO23.jpg"), plot = Gene_ontology_KO23_DOWN, width = 12, height = 8, dpi = 800)

# Combine Upregulated (A) and Downregulated (B) GO plots vertically
GO_Combined_KO23 <- ggarrange(
  Gene_ontology_KO23_UP, 
  Gene_ontology_KO23_DOWN,
  labels = c("A", "B"),   # Figure A = Upregulated, Figure B = Downregulated
  ncol = 1, nrow = 2,     # Vertical layout
  common.legend = TRUE,   # Shared legend
  legend = "right"        # Legend on the right
)

# Save combined figure
ggsave(file.path(plot_dir, "GO_Terms_KO23_Combined.jpg"), plot = GO_Combined_KO23,
       width = 12, height = 14, dpi = 800)