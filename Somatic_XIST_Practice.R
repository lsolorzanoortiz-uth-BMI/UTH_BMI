# Set CRAN mirror to avoid prompts during installation #####
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install WGCNA and dependencies from CRAN #####
install.packages("WGCNA", dependencies = TRUE)

install.packages(c(
  "BiocManager",      # Manages the installation of Bioconductor 
  "ggplot2",          # For vizualization
  "gplots",           # For heatmap.2
  "reshape2",         # For data manipulation
  "igraph",           # For network analysis and visualization
  "pheatmap",         # For heatmap visualization
  "RColorBrewer",     # For color palettes
  "corrplot",         # For correlation plots
  "ggrepel" ,         # For non-overlapping labels in plots
  "doParallel",       # For parallel processing
  "dplyr",            # Data manipulation and transformation
  'GOplot',           # Enrichment Analysis
  "gggenes",          # Plotting 
  "circlize",
  "viridisLite",
  "ggridges",         # For visualizing the GSEA 
  'VennDiagram'
))


# Install Bioconductor packages #####
BiocManager::install(c(
  "limma",            # For normalization and QC
  "DESeq2",           # For data transformation
  "clusterProfiler",  # For functional enrichment analysis
  "org.Hs.eg.db",     # Human gene annotations
  "GO.db",            # Gene Ontology database
  "STRINGdb",         # For protein-protein interaction networks
  "impute",           # Missing value imputation
  "preprocessCore",   # Low-level Preprocessing and Normalization
  "GO.db",            # Gene Ontology (GO) Database
  "AnnotationDbi",    # Database Interface/Query Tool
  "GEOquery",         # Downloads Data from NCBI Geodataset
  "biomaRt",          # To retrieve ESEMBL IDs
  "GO.db",             # To retrieve Go ontology enrichment analysis
  "karyoploteR",
  "Gviz",
  "EnhancedVolcano",
  "ComplexHeatmap"
))

# Load required libraries #####
library(WGCNA)
library(ggplot2)
library(gplots)
library(ggrepel)
library(reshape2)
library(pheatmap)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(readr)
library(doParallel)
library(GEOquery)
library(biomaRt)
library(dplyr) 
library(limma)
library(GO.db)
library(GOplot)
library(karyoploteR)
library(Gviz)
library(gggenes)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(ggridges)
library(VennDiagram)

# Increasing timeout cutoff 
options(timeout = 1200)

# Set working directory #####
setwd("/Users/luzsmac/Desktop/Intro_to_Bioinformatics /Practice_Winter_Break/SLE_Practice _Somatic_Cells/GSE305810")

# Loading Data #####
getGEOSuppFiles("GSE305810")

# List all files in the downloaded folder #####
files <- list.files("GSE305810", full.names = TRUE)

# Read the count data file (In RPKM) #####
count_data <- read.delim(files[1])

# Extract sample information (metadata) #####
sample_info <- getGEO("GSE305810", GSEMatrix = TRUE)
sample_info <- pData(sample_info[[1]])

# Number of columns #####
last_col <- ncol(count_data)

# Filter rows where the sum of columns 5 through the end is > 0 #####
count_data <- count_data[rowSums(count_data[, 5:last_col], na.rm = TRUE ) >0, ] 
  
# Naming the rows #####                                                               
rownames(count_data) <- count_data$Gene

# Initializing useMart #####
mart <- useMart("ensembl", dataset ="hsapiens_gene_ensembl")

# Retrieving Ensembl Ids 
ensembl_ids <- getBM(
  filters = "external_gene_name",
  attributes = c("external_gene_name","ensembl_gene_id"),
  values =count_data$Gene,
  mart = mart
)

# Sub-setting for the Genes with ensembl id #####
count_data <- count_data[ensembl_ids$external_gene_name, ]

# Checking to see if they are identical #####
identical(count_data$Gene, ensembl_ids$external_gene_name)

# Renaming #####
Ensembl_ID <- ensembl_ids$ensembl_gene_id

# Adding the Ensmbl ID to the dataframe #####
count_data <- cbind(ensembl_ids$ensembl_gene_id, count_data)

# Cleaning sample info #####
sample_info_cleaned <- data.frame(
  sample_id = sample_info$title,
  treatment = gsub("^[^,]+,([^,]+),.+", "\\1" ,sample_info$title)
)
rownames(sample_info) = sample_info$sample_id

######### LIMMA########

# Subsetting to get the RNA-seq reads only #####
matrix <- count_data[ ,6:ncol(count_data)]  

# Metadata #####
gene_metadata <- count_data[, 1:5]

group_names = factor(c(
  rep("WT",3),
  rep("XIST_KD",4)
))

# Log Transform #####
log_matrix <- as.matrix(log2(matrix+1))

# Construct Matrix #####
design <- model.matrix(~0 + group_names)
colnames(design) <- levels(group_names)

# Fitting #####
fit <- lmFit(log_matrix, design, genes = gene_metadata)

# Contrast #####
contrast_matrix <- makeContrasts(
  WT_vs_XISTKD = WT-XIST_KD,
  levels = design
)

# Fitting with contrasts #####
fit2 <- contrasts.fit(fit, contrast_matrix)

# Fitting Ebayes #####
fit3 <- eBayes(fit2, trend = TRUE, robust = TRUE)
# Summary of results #####
results <- decideTests(fit3)
summary(results)
test <- as.data.frame(fit3)
# Getting chrom metadata #####
chrom_metadata <- count_data$Chromosome


# Results for Limma #####
all_results <- topTable(
  fit3,
  coef = "WT_vs_XISTKD",
  sort.by = "none",
  n = Inf,
  adjust.method = "BH"
)

# Plotting the differential analysis results 
all_results$status <- "Not significant"

all_results$status[all_results$adj.P.Val < 0.05 & all_results$logFC > 0] <- "Upregulated"
all_results$status[all_results$adj.P.Val < 0.05 & all_results$logFC < 0] <- "Downregulated"

counts <- table(all_results$status)

png("plots/DEG_counts_text.png",
    width = 1200, height = 800, res = 150)

# blank canvas with nicer limits
plot(1, type = "n",
     xlim = c(0, 1), ylim = c(0, 1),
     axes = FALSE, xlab = "", ylab = "",
     main = "Differential Expression Summary GSE305810 WT vs XIST KD hTERT RPE-1")

# optional subtitle-style total
total_genes <- sum(counts)

text(0.5, 0.85,
     paste0("Total genes analyzed: ", total_genes),
     cex = 1.4, col = "black")

# Upregulated
text(0.5, 0.65,
     paste0("Upregulated: ", counts["Upregulated"]),
     col = "#D62728", cex = 2.2, font = 2)

# Downregulated
text(0.5, 0.45,
     paste0("Downregulated: ", counts["Downregulated"]),
     col = "#1F77B4", cex = 2.2, font = 2)

# Not significant
text(0.5, 0.25,
     paste0("Not significant: ", counts["Not significant"]),
     col = "gray40", cex = 2.0, font = 2)

dev.off()


# Transforming into X, Y or Other #####
custom_chrom <- function(x){
  if(x == "X"){
    return("X")
  }
  else if(x =="Y"){
    return("Y")
  }
  else{
    return("Other")
  }
}

all_results['Chrom'] <- sapply(chrom_metadata, custom_chrom)

# Plotting DEGS####
png("plots/Differentially_Enhaced_Volcano_GSE305810.png", width = 1000, height =500, res =100 )
EnhancedVolcano(
  all_results,
  x = 'logFC',
  y = 'adj.P.Val',
  lab = rownames(all_results),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 2,
  labSize = 3,
  title = 'Differentially Expressed Genes in GSE305810',
  subtitle = 'WT vs XISTKD using limma'
)
dev.off()

# Color coding by chrom ###
keyvals.col <- ifelse(all_results$Chrom == "X", "forestgreen","lightcoral")

names(keyvals.col) <- all_results$Chrom


# Odering ddf by X so that the volcano plot will show the X color on top 
all_results <- all_results[order(all_results$Chrom == "X"),]
keyvals.col <- sort(keyvals.col, decreasing= TRUE)

# Plotting, chrom by color ####                                
 png("plots/Differentially_Enhaced_Volcano_Colored_by_Chrom_GSE305810.png", width = 1000, height =500, res =100)
 p <-EnhancedVolcano(
  all_results,
  x = 'logFC',
  y = 'adj.P.Val',
  lab = rownames(all_results),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 2,
  labSize = 3,
  colCustom = keyvals.col,
  shape = 6,
  colAlpha = 1,
  title = 'Differentially Expressed Genes in GSE305810 Colored by Chromosome',
  subtitle = 'WT vs XISTKD using limma'
)
p + theme(plot.title = element_text(hjust = 0))
dev.off()

# Plotting the differential analysis results 
all_results_x_chrom <- subset(all_results, Chrom == "X")

all_results_x_chrom$status <- "Not significant"

# Upregulated (significant + positive logFC)
all_results_x_chrom$status[
  all_results_x_chrom$adj.P.Val < 0.05 & all_results_x_chrom$logFC > 0
] <- "Upregulated (FDR + logFC)"

# Downregulated (significant + negative logFC)
all_results_x_chrom$status[
  all_results_x_chrom$adj.P.Val < 0.05 & all_results_x_chrom$logFC < 0
] <- "Downregulated (FDR + logFC)"

counts <- table(all_results_x_chrom$status)
total_genes <- sum(counts)

png("plots/DEG_counts_text_X_chom.png",
    width = 1200, height = 800, res = 150)
plot(1, type = "n",
     xlim = c(0, 1), ylim = c(0, 1),
     axes = FALSE, xlab = "", ylab = "",
     main = "Differential Expression GSE305810 WT vs XIST KD hTERT RPE-1 (X Chrom)")

text(0.5, 0.85,
     paste0("Total genes analysed: ", total_genes),
     cex = 1.4, col = "black")

# Upregulated
text(0.5, 0.65,
     paste0("Upregulated (FDR + logFC): ", counts["Upregulated (FDR + logFC)"]),
     col = "#D62728", cex = 2.0, font = 2)

# Downregulated
text(0.5, 0.45,
     paste0("Downregulated (FDR + logFC): ", counts["Downregulated (FDR + logFC)"]),
     col = "#1F77B4", cex = 2.0, font = 2)

# Not significant
text(0.5, 0.25,
     paste0("Not significant: ", counts["Not significant"]),
     col = "gray40", cex = 2.0, font = 2)

dev.off()

# Getting significant results #####
sig_results <- all_results[which(all_results$adj.P.Val<0.05), ]

# Subsetting for sig on count data #####
count_data_sig <- count_data[rownames(sig_results), ]

# Saving sig count_data #####
write.table(count_data_sig, "sig_count_data", sep ="\t")

# Matrix sig #####
matrix_sig <- matrix[rownames(sig_results), ]

matrix_sig_scaled <- t(scale(t(matrix_sig)))

# Create annotation dataframe for samples #####
annotation_col <- data.frame(
  Group = sample_info_cleaned$treatment,
  row.names = sample_info_cleaned$sample_id
)

# Hierarchical tree of the samples #####
hc_samples <- hclust(as.dist(1-cor(matrix_sig_scaled ,
                           method = "pearson")), method = "complete")
# Getting dendrogram
sampleTree = as.dendrogram(hc_samples, method = average)

# Plotting Samples #####
png("plots/Sample_clustering_GSE305810.png", width = 800, height = 600, res = 100)

plot(sampleTree,
     main = "Sample Clustering in GSE305810 dataset ")

dev.off()

# Hierarchical tree of genes #####
hc_genes <- hclust(as.dist(1-cor(t(matrix_sig_scaled) ,
                           method ="pearson")), method = "complete")
# Getting dendrogram #####
geneTree = as.dendrogram(hc_genes, method = average)

# Plotting Samples #####
png("plots/Gene_clustering_GSE305810.png", width = 800, height = 600, res = 100)

plot(geneTree,
     main = "Gene Clustering in GSE305810 dataset ",
     leaflab = "none",
     ylab = "Height")

dev.off()

# Heatmap #####
png("plots/Heatmap_GSE305810.png", width = 800, height = 600, res = 100)

heatmap.2(as.matrix(matrix_sig),
          Rowv = as.dendrogram(hc_genes),
          Colv = as.dendrogram(hc_samples),
          col = redgreen(100),
          scale ="row",
          margins = c(7,7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap of GSE305810 dataset",
          trace = "none"
)

dev.off()

# Clustering #####
n_clusters <- 2
gene_clusters <- cutree(hc_genes, k = n_clusters)

# Cluster dataframe #####
cluster_assignments <- data.frame(
  Gene = names(gene_clusters),
  Cluster = gene_clusters
)


# Reorder genes by cluster assignment #####
gene_order <- order(gene_clusters)
deg_scaled_ordered <- matrix_sig_scaled[gene_order, ]
cluster_ordered <- gene_clusters[gene_order]

# Create cluster annotation for genes #####
cluster_annotation <- data.frame(
  Cluster = factor(cluster_ordered),
  row.names = rownames(deg_scaled_ordered),
  SYMBOL = rownames(deg_scaled_ordered)
)

# Calculate gap positions (where clusters change) #####
gap_positions <- which(diff(cluster_ordered) != 0)

# Create annotated heatmap with gaps between clusters #####
png("plots/Heatmap_GSE305810_gapped.png", width = 800, height = 600, res = 100)
pheatmap(
  deg_scaled_ordered,
  cluster_rows = FALSE,  
  cluster_cols = hc_samples,
  annotation_row = NA,
  annotation_col = NA,
  gaps_row = gap_positions,  
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "DEG Hierarchical Clusters (Separated by Gaps)",
  fontsize = 10,
  filename = "hierarchical_clustering_DEGs_gapped.png",
  width = 10,
  height = 8
)
dev.off()

# Enrichment analysis
# Convert gene symbols to Entrez IDs 
gene_entrez <- bitr(cluster_annotation$SYMBOL, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Setting the row names to the symbols #####
row.names(gene_entrez) <- gene_entrez$SYMBOL

# 3.38% of genes not mapped. Getting the missing symbols #####
missing_genes <- setdiff(cluster_annotation$SYMBOL, gene_entrez$SYMBOL)

# Saving Symbols not mapped for documentation purposes #####
write.table(missing_genes, "Genes_not_mapped_to_ENTREZID.txt")

# subset for cluster 1 #####
cluster_1 <- subset(cluster_annotation, Cluster == 1)
cluster_2 <- subset(cluster_annotation, Cluster == 2)

gene_entrez_cluster_1 <- gene_entrez[cluster_1$SYMBOL, ]
gene_entrez_cluster_2 <- gene_entrez[cluster_2$SYMBOL, ]

# Perform GO enrichment for KD-associated module #####
go_enrichment_cluster_1 <- enrichGO(
  gene = gene_entrez_cluster_1$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",              
  pAdjustMethod = "BH",    
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

png("plots/kegg_cluster1.png", width = 800, height = 600, res = 100)
go_enrichment_kegg_cluster1 <- enrichKEGG(
  gene = gene_entrez_cluster_1$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",    
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
dev.off()

go_enrichment_cluster_2 <- enrichGO(
  gene = gene_entrez_cluster_2$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",              
  pAdjustMethod = "BH",    
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
count_data[,"SYMBOL"] <- count_data[,"Gene"]

#GO Enrichment Barplot #####
png("plots/GO_enrichment_barplot_cluster_1_GSE305810.png",
    width = 1000, height = 800, res = 100)
barplot(go_enrichment_cluster_1, 
        showCategory = 15,
        title = "GO Enrichment in Cluster 1")
dev.off()


png("plots/GO_enrichment_barplot_cluster_2_GSE305810.png",
    width = 1000, height = 800, res = 100)
barplot(go_enrichment_cluster_2, 
        showCategory = 15,
        title = "GO Enrichment in Cluster 2")
dev.off()


# Entrez ID for count_data 
count_data_gene_entrez <- bitr(count_data$SYMBOL, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Check for NA and duplicates 0.01% not mapped 
count_data_gene_entrez<- count_data_gene_entrez[!is.na(count_data_gene_entrez$ENTREZID),]
count_data_gene_entrez<- count_data_gene_entrez[!duplicated(count_data_gene_entrez$ENTREZID),]
count_data_gene_entrez<- count_data_gene_entrez[!duplicated(count_data_gene_entrez$SYMBOL),]

# Rename rows on gee entrez count data
rownames(count_data_gene_entrez) <- count_data_gene_entrez$SYMBOL

# Prepare for enrichment analysis 
count_data_gene_entrez <- cbind(count_data_gene_entrez[rownames(count_data_gene_entrez),], all_results[rownames(count_data_gene_entrez), ]) 
gene_entrex_t <- count_data_gene_entrez$t
names(gene_entrex_t) <- count_data_gene_entrez$ENTREZID

# Gse enrichment 
gsea_do <- gseDO(
  geneList = sort(gene_entrex_t, decreasing= TRUE),
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# Filtering by sig 
sig_gsea <- filter(gsea_do, p.adjust < 0.05)

# Plotting the gsea
png("plots/gsea_do_enrichment.png", width = 2000,
    height = 4000,res = 100)
ridgeplot(sig_gsea, showCategory = 35, orderBy= "NES", fill= "p.adjust")+ 
  scale_fill_continuous(low="#FF0000CC", high="#3182bdCC") +
  ggtitle("GSEA Results in GSE305810 WT vs XIST KDhTERT RPE-1") +
  xlab("NES")+
 theme(text = element_text(size = 26)) +
  theme(axis.text.y = element_text(size = 26))

dev.off()


# Df of differentially expressed genes and chroms #####
sig_genes_df <- data.frame(
  Symbol = rownames(count_data_sig),
  Chrom = count_data_sig$Chromosome
)
# Extracting coordinates #####
mart <- useMart("ensembl", dataset ="hsapiens_gene_ensembl")

# Retrieve coordinates #####
sig_genes_coordinates <- getBM(attributes = c("external_gene_name", "chromosome_name", "strand",
                     "start_position", "end_position"),
      filters = c("external_gene_name", "chromosome_name"),
      values = list(sig_genes_df$Symbol, c(1:22, "X", "Y", "MT")),
      mart = mart)
# Formatting to convert to granges #####
sig_genes_coordinates$chromosome_name <- paste0("chr", sig_genes_coordinates$chromosome_name)

# Formating strands #####
formatting_strands <- function(x){
  if(x==1){
    return("+")
  }
  else if (x==-1){
    return("-")
  }
  else{
    return("*")
  }
}

# Calling formatting_strands #####
sig_genes_coordinates$strand<-lapply(sig_genes_coordinates$strand, formatting_strands)

# Df of unique genes #####                        
sig_genes_coordinates_unique <- sig_genes_coordinates[!duplicated(sig_genes_coordinates$external_gene_name), ]

# Indices of duplicated genes #####
duplicated_indices <- which(duplicated(sig_genes_coordinates$external_gene_name))

# Df of duplicated genes #####
df_duplicated <- sig_genes_coordinates[duplicated_indices, ]

# Define chrom order #####
chromosome_order <- factor(paste0("chr", c(1:22, "X", "Y")))

# Order df #####
sig_genes_coordinates_unique$chromosome_name <- factor(
  sig_genes_coordinates_unique$chromosome_name, 
  levels = chromosome_order
)
# Plotting genes and coordinates #####
png("plots/genes_coordinates.png", width = 4500, height = max(2000, 23 * 250), res = 300)
ggplot(sig_genes_coordinates_unique, aes(xmin = start_position, xmax = end_position,
                                         y= chromosome_name)) +
  geom_gene_arrow() +
  facet_wrap(~ chromosome_name, ncol = 1) +
  theme(
    # Removes the text labels on the Y axis
    axis.text.y = element_blank(),
    # Removes the little tick marks on the Y axis
    axis.ticks.y = element_blank(),
    # Removes the labels on the side/top of the facets (strips)
    strip.text.y = element_blank()
  )+
  ggtitle("Differentially Expressed Genes In WT vs XISTKD In GSE305810 By Chromosome")

dev.off()

# subsetting and saving data with X chrom #####
sig_genes_coordinates_unique_X_chrom <- filter(data.frame(sig_genes_coordinates_unique), chromosome_name == "chrX")

# Formatting 
sig_genes_coordinates_unique_X_chrom[,'strand' ] <- as.character(sig_genes_coordinates_unique_X_chrom[,'strand' ])

# Saving 
write.table(sig_genes_coordinates_unique_X_chrom , file = "sig_genes_coordinates_unique_X_chrom",  col.names=TRUE, row.names= FALSE)


