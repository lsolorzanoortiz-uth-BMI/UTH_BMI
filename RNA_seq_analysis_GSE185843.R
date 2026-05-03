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
  "data.table",       # For data handling
  "tidyverse"         # For data manipulation

))

BiocManager::install("DESeq2")

# Install Bioconductor packages #####
BiocManager::install(c(
  "DiffBind",                             # For diff analysis
  "ChiPseeker",                           # For Genomic annotations
  "TxDb.Mmusculus.UCSC.mm10.knownGene",   # Mouse genome annotations 
  "org.Mm.eg.db",                         # Mouse gene annotations 
  "Limma",                                # For diff analysis 
  "DESeq2"
  
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
library(DESeq2)
library(tidyr)
library(tidyverse)

# Increasing time allowed for compiling ####
options(timeout = 1200)


setwd("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE185843_Analysis")

if(!dir.exists("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE185843_Analyais")){
  getGEOSuppFiles("GSE185843")
}

# List all files in the downloaded folder ####
if(dir.exists("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE185843_Analysis/GSE185843")){
  files <- list.files(path="/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE185843_Analysis/GSE185843/GSE185843_RAW/DOM_CAST", full.names = TRUE)


  # Extract sample information (metadata) ####
  sample_info <- getGEO("GSE185843", GSEMatrix = TRUE)
  sample_info <- pData(sample_info[[1]])
}

files_df <-as.data.frame(files)
sample_info_df <- as.data.frame(sample_info)

library(readr)

# Count data ####
count_data_list <- lapply(seq_along(files), function(i) {
    read_delim(
      files[i],
      delim = "\t",
      escape_double = FALSE,
      skip = 1
    )
  })
    
aggregation_of_count_data<- function(list){
  for( i in seq_along(list)){
    if( i==1){
      count_data <- as.data.frame(list[[i]])
      rownames(count_data) <- count_data$Geneid
    }
    else{
      current_count_data <- as.data.frame(list[[i]])
      equal_boolean <-all.equal(rownames(count_data), current_count_data$Geneid)
      
      if(!equal_boolean){
        rownames(current_count_data) <- current_count_data$Geneid
        current_count_data <- current_count_data[rownames(count_data), ]
        equal_row_names <- all.equal(rownames(count_data), rownames(current_count_data))
        
        if(!equal_row_names){
          print("Geneid column is not equal")
          break 
        }
      }
      if(equal_boolean){
      count_data <- cbind(count_data, current_count_data[, 7: ncol(current_count_data)])
      }
      
      }}
  return(count_data)
  }

count_data <-aggregation_of_count_data(count_data_list)

# Renaming stage cols ####
rename_d_columns <- function(df) {
  cols <- colnames(df)
  
  # Only modify columns 7+
  target_cols <- cols[7:length(cols)]
  
  # Extract d + number (e.g., d12)
  d_vals <- sub(".*(A11B2|C7H8).*(d[0-9]+).*", "\\1_\\2", target_cols)
  
  # Create replicate indices (1,2,3 repeating)
  reps <- rep(1:3, length.out = length(target_cols))
  
  # New names
  new_names <- paste0(d_vals, "_", reps)
  
  # Replace column names
  cols[7:length(cols)] <- new_names
  colnames(df) <- cols
  
  return(df)
}

count_data <-rename_d_columns(count_data)


# Cleaning the chr column ####
count_data$Chr <- sub("([chr0-9XY]*).*", "\\1", count_data$Chr)

# Cleaning the positional columns ####
start_cols_data <- count_data %>% separate_wider_delim(Start, delim = ";", names_sep = "",
                                                       too_few="align_start", names_repair = "universal" )
count_data[, "Start"]<-apply(start_cols_data[, 3:24], 1, min, na.rm = TRUE)

end_cols_data <- count_data %>% separate_wider_delim(End, delim = ";", names_sep = "",
                                                     too_few="align_start", names_repair = "universal" )

count_data[, "End"] <-apply(end_cols_data[, 4:25], 1, max, na.rm = TRUE)

# Converting positional cols to numeric for plotting ####
count_data$Start <- as.numeric(count_data$Start)
count_data$End <- as.numeric(count_data$End)

# Sorting the expression columns ####
names(count_data[, 7:ncol(count_data)]) <- sort(names(count_data[, 7:ncol(count_data)]))

# Extracting phenotype information ####
# Extracting # of each cell line
n_A11B2_cols <- length(grep("A11B2", names(count_data)))
n_C7H8_cols <- length(grep("C7H8", names(count_data)))

# Extracting the time data along with the cell line
matched_cols <- sub(".*(A11B2|C7H8).*(d\\d+).*", "\\2", names(count_data[, 7:ncol(count_data)]))

# Making it into a df fo easy view and use 
count <- 0
count_vector <- c()
match_vector <- c()

counting_matches <-function(matched_cols){
for( i in seq_along(matched_cols)){
  current_match <- matched_cols[i]
  
  if( i== 1){
  count = count+1 
  match_vector <- c(match_vector, current_match)
  }
  
  else if ((i != 1) && (current_match == matched_cols[i-1])){
    count = count+1
  }
  else if( (i != 1) && current_match != matched_cols[i-1]){
    count_vector <- c(count_vector, count)
    match_vector <- c( match_vector, current_match)
    count = 1
  }
}
count_vector <- c(count_vector, count)
match_info <- setNames(count_vector, match_vector)
return (match_info)
}
match_info <-counting_matches(matched_cols)

# Subsetting only the counts ####
count_matrix <- count_data[7:ncol(count_data)]

# Creating metadata info ####
sample_info_metadata <- data.frame( 
  sample = factor(colnames(count_matrix)),
  genotype = factor(rep(c("A11B2", "C7H8"), times = c(36, 42))),
  stage = factor(rep(names(match_info), times = unname(match_info))))

rownames(sample_info_metadata) <- sample_info_metadata$sample

# Deseq for DOM and CAST ####
dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                              colData = sample_info_metadata, 
                              design = ~genotype+stage)
dds <- DESeq(dds)

# Results for DOM and CAST ####
results_deseq2 <- as.data.frame(results(dds, contrast = c("genotype", "A11B2", "C7H8")))

results_deseq2_sig <-subset(results_deseq2, padj<0.05)
genomic_info_sig <- count_data[rownames(results_deseq2_sig), 1:6]

# Define chrom order #####
chromosome_order <- factor(paste0("chr", c(1:22, "X", "Y")))

# Order df #####
genomic_info_sig$Chr <- factor(
  genomic_info_sig$Chr, 
  levels = chromosome_order)




# Plotting DOM and Cast results ####
png("plots/DOM_CAST_CHR_RNA_seq.png", width = 5500, height = max(2500, 23 * 250), res = 300)
ggplot(genomic_info_sig, aes(xmin = Start, xmax = End, y= Chr )) +
  geom_gene_arrow() +
  facet_grid(Chr ~.)+
  theme(
    # Removes the text labels on the Y axis
    axis.text.y = element_blank(),
    # Removes the little tick marks on the Y axis
    axis.ticks.y = element_blank(),
    # Removes the labels on the side/top of the facets (strips)
    #strip.text.y = element_blank()
  )+
  ggtitle("Differential Positions in Chr-RNA-seq in DOM vs CAST NPC mice cells (GSE185843)")

dev.off()

# Subset for Xchrom (Dom Vs Cast ) ####
XChr_genomic_sig_DOM_CAST_subset <- subset(genomic_info_sig, Chr== "chrX")

# Save 
XChr_genomic_sig_DOM_CAST <- write.table(XChr_genomic_sig_DOM_CAST_subset, "XChr_genomic_sig_DOM_CAST", col.names = TRUE, , quote=FALSE, sep="\t", row.names = FALSE)

# Subset for Xchrom (Dom Vs Cast ) ####
XChr_deseq2_results<- results_deseq2_sig[rownames(XChr_genomic_sig_DOM_CAST_subset),]

# Save 
write.table(XChr_deseq2_results, "XChr_deseq2_results", col.names = TRUE, , quote=FALSE, sep="\t", row.names = TRUE)


# Differential by stage - DOM(A11B2) ####

dds_stage_A11B2 <- DESeqDataSetFromMatrix(countData = count_matrix[, 1:36], 
                                          colData = sample_info_metadata[1:36, ], 
                                          design = ~stage)

dds_stage_A11B2 <- DESeq(dds_stage_A11B2)

results_stage_A11B2 <- as.data.frame(results(dds_stage_A11B2, contrast = c("stage", "d35", "d2")))

subset_sig_A11B2_results <- subset(results_stage_A11B2, padj <0.05)


# Subset for Xchrom (Stages) ####

genomic_info_stages_results_A11B2 <- count_data[rownames(subset_sig_A11B2_results), 1:6]

write.table(genomic_info_stages_results_A11B2, "genomic_info_stages_results_DOM", row.names= FALSE, col.names= TRUE, quote = FALSE, sep= "\t")



table(count_data["Mir22hg", ])
dom_subset <- data.frame(count_data["Mir22hg", 7:42])
cast_subset <- data.frame(count_data["Mir22hg", 43:ncol(count_data)])
rowMeans(dom_subset)
rowMeans(cast_subset)

