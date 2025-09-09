#KEGG Enrichment Analysis=======================================================
##Load necessary libraries------------------------------------------------------
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

##Function to Convert Gene Symbols to Entrez IDs--------------------------------
convert_to_entrez <- function(gene_list) {
  entrez_ids <- mapIds(org.Hs.eg.db, 
                       keys = gene_list, 
                       column = "ENTREZID", 
                       keytype = "SYMBOL", 
                       multiVals = "first")
  
  # Remove NA values
  return(entrez_ids[!is.na(entrez_ids)])
}

##Function to Perform KEGG Enrichment-------------------------------------------
plot_kegg_ggplot <- function(kegg_res, title, filename, top_n = 20) {
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res@result)) > 0) {
    kegg_df <- as.data.frame(kegg_res@result)
    kegg_df <- head(kegg_df[order(kegg_df$p.adjust), ], top_n)
    
    # Calculate GeneRatio as numeric
    kegg_df$GeneRatio <- sapply(kegg_df$GeneRatio, function(x) eval(parse(text = x)))
    
    p <- ggplot(kegg_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
      geom_point() +
      scale_color_gradient(low = "red", high = "blue") +
      labs(title = title, x = "Gene Ratio", y = "KEGG Pathway", 
           color = "Adjusted p-value", size = "Gene Count") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)
      )
    
    # Save plot
    ggsave(file.path(plot_dir, filename), plot = p, width = 10, height = 6, dpi = 800)
    
    return(p)  # return ggplot object so we can combine later
  } else {
    print(paste("No significant KEGG pathways found for", title))
    return(NULL)
  }
}

## Create KEGG plots for KO22
KEGG_Up_Plot_KO22   <- plot_kegg_ggplot(kegg_up_KO22, "KEGG Pathway for KO22 - Upregulated Genes", "KEGG_Upregulated_KO22.png")
KEGG_Down_Plot_KO22 <- plot_kegg_ggplot(kegg_down_KO22, "KEGG Pathway for KO22 - Downregulated Genes", "KEGG_Downregulated_KO22.png")

# Combine KO22 Up & Down vertically
KEGG_Combined_KO22 <- ggarrange(
  KEGG_Up_Plot_KO22, 
  KEGG_Down_Plot_KO22,
  labels = c("A", "B"),
  ncol = 1, nrow = 2,
  common.legend = TRUE,
  legend = "right"
)

ggsave(file.path(plot_dir, "KEGG_KO22_Combined.jpg"), plot = KEGG_Combined_KO22,
       width = 16, height = 20, dpi = 800)


## Create KEGG plots for KO23
KEGG_Up_Plot_KO23   <- plot_kegg_ggplot(kegg_up_KO23, "KEGG Pathway for KO23 - Upregulated Genes", "KEGG_Upregulated_KO23.png")
KEGG_Down_Plot_KO23 <- plot_kegg_ggplot(kegg_down_KO23, "KEGG Pathway for KO23 - Downregulated Genes", "KEGG_Downregulated_KO23.png")

# Combine KO23 Up & Down vertically
KEGG_Combined_KO23 <- ggarrange(
  KEGG_Up_Plot_KO23, 
  KEGG_Down_Plot_KO23,
  labels = c("A", "B"),
  ncol = 1, nrow = 2,
  common.legend = TRUE,
  legend = "right"
)

ggsave(file.path(plot_dir, "KEGG_KO23_Combined.jpg"), plot = KEGG_Combined_KO23,
       width = 16, height = 20, dpi = 800)

##Creating .csv-----------------------------------------------------------------
#Convert Gene Symbols to Entrez IDs
up_entrez_KO22 <- convert_to_entrez(up_genes_FIH22KO)
down_entrez_KO22 <- convert_to_entrez(down_genes_FIH22KO)

#Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO22 <- enrichKEGG(gene = up_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO22 <- enrichKEGG(gene = down_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)

#Save Results
write.csv(as.data.frame(kegg_up_KO22@result), file.path(data_dir, "KEGG_Upregulated_KO22.csv"))
write.csv(as.data.frame(kegg_down_KO22@result), file.path(data_dir, "KEGG_Downregulated_KO22.csv"))

#Convert Gene Symbols to Entrez IDs
up_entrez_KO23 <- convert_to_entrez(up_genes_FIH22KO)
down_entrez_KO23 <- convert_to_entrez(down_genes_FIH22KO)

#Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO23 <- enrichKEGG(gene = up_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO23 <- enrichKEGG(gene = down_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)

#Save Results
write.csv(as.data.frame(kegg_up_KO23@result), file.path(data_dir, "KEGG_Upregulated_KO23.csv"))
write.csv(as.data.frame(kegg_down_KO23@result), file.path(data_dir, "KEGG_Downregulated_KO23.csv"))