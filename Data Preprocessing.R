#PREPROCESSING==================================================================
##Installing the necessary packages---------------------------------------------
install.packages("ggplot2")       
install.packages("tidyverse")      
install.packages("factoextra")     
install.packages("ggpubr")         
install.packages("data.table")     
install.packages("limma")          
install.packages("ggrepel")        
install.packages("VennDiagram")    
install.packages("scales")         

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")           
BiocManager::install("biomaRt")           
BiocManager::install("edgeR")             
BiocManager::install("clusterProfiler")   
BiocManager::install("org.Hs.eg.db")      
BiocManager::install("AnnotationDbi")     
BiocManager::install("GSVA")             
BiocManager::install("msigdbr")          
BiocManager::install("GSEABase")          
BiocManager::install("ComplexHeatmap")   


##Create output directories-----------------------------------------------------
plot_dir <- "Analysis plots with combined code using Rowsum"
data_dir <- "List_for_analysis with combined code using Rowsum"

dir.create(plot_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)

##Load necessary libraries------------------------------------------------------
library(data.table)  # For fread()
library(dplyr)       # For data manipulation (select, left_join, distinct, etc.)
library(biomaRt)     # For mapping Ensembl IDs to gene symbols
library(edgeR)
library(DESeq2)
##Read files into R-------------------------------------------------------------
Control1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A1_featurecounts.txt")
Control2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A2_featurecounts.txt")
Control3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A3_featurecounts.txt")
FIH22KO1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B1_featurecounts.txt")
FIH22KO2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B2_featurecounts.txt")
FIH22KO3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B3_featurecounts.txt")
FIH23KO1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C1_featurecounts.txt")
FIH23KO2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C2_featurecounts.txt")
FIH23KO3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C3_featurecounts.txt")

##Set column names--------------------------------------------------------------
colnames(Control1)[7] <- "Control1_Count"
colnames(Control2)[7] <- "Control2_Count"
colnames(Control3)[7] <- "Control3_Count"
colnames(FIH22KO1)[7] <- "FIH22KO1_Count"
colnames(FIH22KO2)[7] <- "FIH22KO2_Count"
colnames(FIH22KO3)[7] <- "FIH22KO3_Count"
colnames(FIH23KO1)[7] <- "FIH23KO1_Count"
colnames(FIH23KO2)[7] <- "FIH23KO2_Count"
colnames(FIH23KO3)[7] <- "FIH23KO3_Count"
colnames(Control1)[1] <- "GeneID"
colnames(Control2)[1] <- "GeneID"
colnames(Control3)[1] <- "GeneID"
colnames(FIH22KO1)[1] <- "GeneID"
colnames(FIH22KO2)[1] <- "GeneID"
colnames(FIH22KO3)[1] <- "GeneID"
colnames(FIH23KO1)[1] <- "GeneID"
colnames(FIH23KO2)[1] <- "GeneID"
colnames(FIH23KO3)[1] <- "GeneID"

##Merge data--------------------------------------------------------------------
combined_data <- Control1[, .(GeneID, Control1_Count)]
combined_data <- merge(combined_data, Control2[, .(GeneID, Control2_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, Control3[, .(GeneID, Control3_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH22KO1[, .(GeneID, FIH22KO1_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH22KO2[, .(GeneID, FIH22KO2_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH22KO3[, .(GeneID, FIH22KO3_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH23KO1[, .(GeneID, FIH23KO1_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH23KO2[, .(GeneID, FIH23KO2_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH23KO3[, .(GeneID, FIH23KO3_Count)], by = "GeneID", all = TRUE)
combined_data <- as.data.frame(combined_data)
rownames(combined_data) <- combined_data$GeneID

##Filter genes with zero counts across all samples------------------------------ (USING ROW SUMS, REMOVES ROW IF THE ENTIRE ROW = ZERO)
zero_counts <- rowSums(combined_data == 0)
combined_data <- combined_data[zero_counts <= 8, ]

##Using Ensembl as the database-------------------------------------------------
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensemblID_to_genesymbol <- combined_data[, 1]

gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensemblID_to_genesymbol,
                      mart = mart)

colnames(gene_mapping) <- c("GeneID", "GeneSymbol")

combined_data <- as.data.frame(combined_data)
combined_data$GeneID <- rownames(combined_data)

#Merge with gene mapping
combined_data <- left_join(combined_data, gene_mapping, by = "GeneID")

#Replace GeneID with GeneSymbol if available
combined_data$FinalName <- ifelse(is.na(combined_data$GeneSymbol), combined_data$GeneID, combined_data$GeneSymbol)

# Keep only the first occurrence of each gene symbol
combined_data <- combined_data %>%
  distinct(FinalName, .keep_all = TRUE)

# Set row names and remove unnecessary columns
rownames(combined_data) <- combined_data$FinalName
combined_data <- combined_data %>% dplyr::select(-GeneID, -GeneSymbol, -FinalName)



##Convert the data to a DGEList object for edgeR to work------------------------
counts <- as.matrix(combined_data)
condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23",3)))
dge <- DGEList(counts = counts, group = condition)

#Filter out low-expressed genes
keep <- filterByExpr(dge, group = condition, min.count = 15, min.total.count = 70, min.prop = 0.2)
dge_filtered <- dge[keep,]
filtered_data <- as.data.frame(dge_filtered$counts)

##Conver the data to the format for DESeq2--------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(filtered_data),
  colData = data.frame(condition),
  design = ~condition
)

#Run DESeq2 normalization and differentail expression analysis
dds <- DESeq(dds)


##Testing for normal distribution-----------------------------------------------
normalized_counts <- assay(rlog(dds))
apply(normalized_counts, 2, function(column_data) {
  ks.test(column_data, "pnorm", mean = mean(column_data), sd = sd(column_data))
})

##Test for Poisson distribution-------------------------------------------------
#Fit a Poisson distribution for each sample and calculate the likelihood
poisson_fit <- function(column_data) {
  lambda <- mean(column_data)  # The λ for Poisson distribution is the mean of the data
  log_likelihood <- sum(dpois(column_data, lambda, log = TRUE))  # Calculate log-likelihood of Poisson distribution
  return(log_likelihood)
}

#Perform the fit for each sample
log_likelihood_results <- apply(normalized_counts, 2, poisson_fit)

#View results
log_likelihood_results

#Save results
write.csv(normalized_counts, file.path(data_dir, "normalized_counts.csv"))