
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
  "data.table"       # For data handling
))

# Install Bioconductor packages #####
BiocManager::install(c(
  "DiffBind",                             # For diff analysis
  "ChiPseeker",                           # For Genomic annotations
  "TxDb.Mmusculus.UCSC.mm10.knownGene",   # Mouse genome annotations 
  "org.Mm.eg.db",                         # Mouse gene annotations 
  "Limma"                                  # For diff analysis 
  
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
library(DiffBind)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(limma)

options(timeout = 1200)


setwd("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis")

if(!dir.exists("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis/GSE249505")){
getGEOSuppFiles("GSE249505")
}

# List all files in the downloaded folder #####
if(dir.exists("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis/GSE249505")){
files <- list.files(path="/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis/GSE249505", full.names = TRUE)

# Read the count data file (In RPKM) #####
count_data <- read.delim(files[1])

# Extract sample information (metadata) #####
sample_info <- getGEO("GSE249505", GSEMatrix = TRUE)
sample_info <- pData(sample_info[[1]])
}

# Cleaning sample info 
sample_info <- sample_info |>
  gsub("ATAC_iXist-ChrX-", "", x = _)

sample_info <- sapply(sample_info, function(x) gsub("ATAC_iXist-ChrX-", "", x))


# Setting rownames 
rownames(count_data) <- count_data$Peak

# Getting counts only 
attack_matrix <- count_data[,7:16]
rownames(attack_matrix) <- count_data$Peak

sample_info_metadata <- data.frame( 
  sample = factor(colnames(attack_matrix)),
  genotype = factor(rep(c("WT","SmD5"), each=5)),
  stage = factor(rep(c("ES_d0","NPC_d3","NPC_d5","NPC_d7","NPC_d15"), 2))
    )

log_attack <- log2( attack_matrix +1)
 
design <- model.matrix(~0 + genotype+stage, data= sample_info_metadata)

colnames(design) <- make.names(colnames(design))

fit <-lmFit(log_attack, design)

contrasts_matrix <- makeContrasts(
  WT_vs_SmD5 = genotypeWT - genotypeSmD5,
  levels = design
)


fit2<- contrasts.fit(fit, contrasts_matrix)
fit2<- eBayes(fit2)

all_results <- topTable(
  fit2,
  coef= "WT_vs_SmD5",
  sort.by="none",
  n=Inf,
  adjust.method="BH"
)
# Subset of sig_results 
sig_results<- subset(all_results, P.Value<0.05)

# Count data of sig results
count_data_sig <- count_data[rownames(sig_results), ]

write.table(count_data, file = "count_data_sig_Cast", row.names = TRUE, col.names=TRUE)



