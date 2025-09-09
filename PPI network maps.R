# Load libraries
library(STRINGdb)
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Set basic parameters
species_tax_id <- 9606
score_cutoff <- 900
img_w <- 10; img_h <- 8; dpi <- 300



##Prepare KO22 upregulated genes------------------------------------------------
df <- as.data.frame(results_plenti_vs_KO22)
df$Gene <- rownames(df)
df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]

df$Significance <- "Not Significant"
df$Significance[df$padj < 0.05 & df$log2FoldChange > log2(2)] <- "Upregulated"
df$Significance[df$padj < 0.05 & df$log2FoldChange < -log2(2)] <- "Downregulated"

deg_up <- subset(df, Significance == "Upregulated")
gene_vec <- unique(deg_up$Gene)

if (length(gene_vec) < 2) stop("Not enough DEGs to build network.")

# Initialize STRINGdb
string_db <- STRINGdb$new(
  version = "11.5",
  species = species_tax_id,
  score_threshold = score_cutoff,
  input_directory = ""
)

# Map to STRING IDs
mapped <- string_db$map(
  data.frame(GENE = gene_vec),
  "GENE",
  removeUnmappedRows = TRUE
)

if (nrow(mapped) < 2) stop("Too few mapped genes.")

# Get interactions
edges <- string_db$get_interactions(mapped$STRING_id)
if (nrow(edges) == 0) stop("No interactions found.")

# Save edge data
write_csv(edges, file.path(data_dir, "KO22_Upregulated_STRING_edges.csv"))

# Create igraph object
g <- graph_from_data_frame(edges[, 1:2], directed = FALSE) %>%
  set_vertex_attr("symbol", value = mapped$GENE[match(V(.)$name, mapped$STRING_id)])

# Detect communities and convert to tbl_graph
g_tbl_up22 <- as_tbl_graph(g) %>%
  mutate(cluster = as.factor(group_louvain()))

# Plot network with clusters
p_up22 <- ggraph(g_tbl_up22, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = cluster), size = 3) +
  geom_node_text(aes(label = symbol), repel = TRUE, size = 3) +
  theme_void() + ggtitle("KO22 Upregulated PPI Network (STRING)") +
  guides(color = "none")

ggsave(file.path(plot_dir, "KO22_Upregulated_STRING_network.png"),
       plot = p_up22, width = img_w, height = img_h, dpi = dpi)

# GO enrichment analysis
ego_up22 <- enrichGO(
  gene = V(g_tbl_up22)$symbol,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write_csv(as_tibble(ego_up22@result),
          file.path(data_dir, "KO22_Upregulated_GO_enrichment.csv"))

##Prepare KO22 downregulated genes-----------------------------------------------
df <- as.data.frame(results_plenti_vs_KO22)
df$Gene <- rownames(df)
df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]

df$Significance <- "Not Significant"
df$Significance[df$padj < 0.05 & df$log2FoldChange > log2(2)] <- "Upregulated"
df$Significance[df$padj < 0.05 & df$log2FoldChange < -log2(2)] <- "Downregulated"

deg_down <- subset(df, Significance == "Downregulated")
gene_vec <- unique(deg_down$Gene)

if (length(gene_vec) < 2) stop("Not enough DEGs to build network.")

# Initialize STRINGdb
string_db <- STRINGdb$new(
  version = "11.5",
  species = species_tax_id,
  score_threshold = score_cutoff,
  input_directory = ""
)

# Map to STRING IDs
mapped <- string_db$map(
  data.frame(GENE = gene_vec),
  "GENE",
  removeUnmappedRows = TRUE
)

if (nrow(mapped) < 2) stop("Too few mapped genes.")

# Get interactions
edges <- string_db$get_interactions(mapped$STRING_id)
if (nrow(edges) == 0) stop("No interactions found.")

# Save edge data
write_csv(edges, file.path(data_dir, "KO22_Downregulated_STRING_edges.csv"))

# Create igraph object
g <- graph_from_data_frame(edges[, 1:2], directed = FALSE) %>%
  set_vertex_attr("symbol", value = mapped$GENE[match(V(.)$name, mapped$STRING_id)])

# Detect communities and convert to tbl_graph
g_tbl_down22 <- as_tbl_graph(g) %>%
  mutate(cluster = as.factor(group_louvain()))

# Plot network with clusters
p_down22 <- ggraph(g_tbl_down22, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = cluster), size = 3) +
  geom_node_text(aes(label = symbol), repel = TRUE, size = 3) +
  theme_void() + ggtitle("KO22 Downregulated PPI Network (STRING)") +
  guides(color = "none")

ggsave(file.path(plot_dir, "KO22_Downregulated_STRING_network.png"),
       plot = p_down22, width = img_w, height = img_h, dpi = dpi)

# GO enrichment analysis
ego_down22 <- enrichGO(
  gene = V(g_tbl_down22)$symbol,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write_csv(as_tibble(ego_down22@result),
          file.path(data_dir, "KO22_Downregulated_GO_enrichment.csv"))

# Save GO enrichment results
write_csv(as_tibble(ego@result),
          file.path(data_dir, "KO22_Downregulated_GO_enrichment.csv"))

p_combined22 <- ggarrange(
  p_up22, p_down22,
  ncol = 1, nrow = 2
)

ggsave(
  filename = file.path(plot_dir, "KO22_STRING_network_combined.png"),
  plot = p_combined22, width = 12, height = 14, dpi = dpi
)

#KO23

##Prepare KO23 upregulated genes------------------------------------------------
df <- as.data.frame(results_plenti_vs_KO23)
df$Gene <- rownames(df)
df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]

df$Significance <- "Not Significant"
df$Significance[df$padj < 0.05 & df$log2FoldChange > log2(2)] <- "Upregulated"
df$Significance[df$padj < 0.05 & df$log2FoldChange < -log2(2)] <- "Downregulated"

deg_up <- subset(df, Significance == "Upregulated")
gene_vec <- unique(deg_up$Gene)

if (length(gene_vec) < 2) stop("Not enough DEGs to build network.")

# Initialize STRINGdb
string_db <- STRINGdb$new(
  version = "11.5",
  species = species_tax_id,
  score_threshold = score_cutoff,
  input_directory = ""
)

# Map to STRING IDs
mapped <- string_db$map(
  data.frame(GENE = gene_vec),
  "GENE",
  removeUnmappedRows = TRUE
)

if (nrow(mapped) < 2) stop("Too few mapped genes.")

# Get interactions
edges <- string_db$get_interactions(mapped$STRING_id)
if (nrow(edges) == 0) stop("No interactions found.")

# Save edge data
write_csv(edges, file.path(data_dir, "KO23_Upregulated_STRING_edges.csv"))

# Create igraph object
g <- graph_from_data_frame(edges[, 1:2], directed = FALSE) %>%
  set_vertex_attr("symbol", value = mapped$GENE[match(V(.)$name, mapped$STRING_id)])

# Detect communities and convert to tbl_graph
g_tbl_up23 <- as_tbl_graph(g) %>%
  mutate(cluster = as.factor(group_louvain()))

# Plot network with clusters
p_up23 <- ggraph(g_tbl_up23, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = cluster), size = 3) +
  geom_node_text(aes(label = symbol), repel = TRUE, size = 3) +
  theme_void() + ggtitle("KO23 Upregulated PPI Network (STRING)") +
  guides(color = "none")

ggsave(file.path(plot_dir, "KO23_Upregulated_STRING_network.png"),
       plot = p_up23, width = img_w, height = img_h, dpi = dpi)

# GO enrichment analysis
ego_up23 <- enrichGO(
  gene = V(g_tbl_up23)$symbol,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write_csv(as_tibble(ego_up23@result),
          file.path(data_dir, "KO23_Upregulated_GO_enrichment.csv"))

##Prepare KO23 downregulated genes-----------------------------------------------
df <- as.data.frame(results_plenti_vs_KO23)
df$Gene <- rownames(df)
df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]

df$Significance <- "Not Significant"
df$Significance[df$padj < 0.05 & df$log2FoldChange > log2(2)] <- "Upregulated"
df$Significance[df$padj < 0.05 & df$log2FoldChange < -log2(2)] <- "Downregulated"

deg_down <- subset(df, Significance == "Downregulated")
gene_vec <- unique(deg_down$Gene)

if (length(gene_vec) < 2) stop("Not enough DEGs to build network.")

# Initialize STRINGdb
string_db <- STRINGdb$new(
  version = "11.5",
  species = species_tax_id,
  score_threshold = score_cutoff,
  input_directory = ""
)

# Map to STRING IDs
mapped <- string_db$map(
  data.frame(GENE = gene_vec),
  "GENE",
  removeUnmappedRows = TRUE
)

if (nrow(mapped) < 2) stop("Too few mapped genes.")

# Get interactions
edges <- string_db$get_interactions(mapped$STRING_id)
if (nrow(edges) == 0) stop("No interactions found.")

# Save edge data
write_csv(edges, file.path(data_dir, "KO23_Downregulated_STRING_edges.csv"))

# Create igraph object
g <- graph_from_data_frame(edges[, 1:2], directed = FALSE) %>%
  set_vertex_attr("symbol", value = mapped$GENE[match(V(.)$name, mapped$STRING_id)])

# Detect communities and convert to tbl_graph
g_tbl_down23 <- as_tbl_graph(g) %>%
  mutate(cluster = as.factor(group_louvain()))

# Plot network with clusters
p_down23 <- ggraph(g_tbl_down23, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = cluster), size = 3) +
  geom_node_text(aes(label = symbol), repel = TRUE, size = 3) +
  theme_void() + ggtitle("KO23 Downregulated PPI Network (STRING)") +
  guides(color = "none")

ggsave(file.path(plot_dir, "KO23_Downregulated_STRING_network.png"),
       plot = p_down23, width = img_w, height = img_h, dpi = dpi)

# GO enrichment analysis
ego_down23 <- enrichGO(
  gene = V(g_tbl_down23)$symbol,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write_csv(as_tibble(ego_down23@result),
          file.path(data_dir, "KO23_Downregulated_GO_enrichment.csv"))

# Save GO enrichment results
write_csv(as_tibble(ego@result),
          file.path(data_dir, "KO23_Downregulated_GO_enrichment.csv"))

p_combined23 <- ggarrange(
  p_up23, p_down23,
  ncol = 1, nrow = 2
)

ggsave(
  filename = file.path(plot_dir, "KO23_STRING_network_combined.png"),
  plot = p_combined23, width = 12, height = 14, dpi = dpi
)

