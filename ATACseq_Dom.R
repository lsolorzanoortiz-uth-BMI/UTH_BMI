
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
  "tidyverse"         # For data handling
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
library(tidyverse)

options(timeout = 1200)


setwd("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis")

if(!dir.exists("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis/GSE249505")){
getGEOSuppFiles("GSE249505")
}

# List all files in the downloaded folder #####
if(dir.exists("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis/GSE249505")){
files <- list.files(path="/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis/GSE249505", full.names = TRUE)

# Read the count data file (In RPKM) #####
count_data <- read.delim(files[2])

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
count_data_sig_Dom <- count_data[rownames(sig_results), ]

count_data_sig_Cast <- read.table("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE249505_Allelic_Analysis/count_data_sig_Cast")

dom_pos <- data.frame("Dom_ATAC_seq", count_data_sig_Dom[, 3:5])
names(dom_pos)<- c("seq_type", "Start", "End", "Strand")
cast_pos <- data.frame("Cast_ATAC_seq", count_data_sig_Cast[, 3:5])
names(cast_pos)<- c("seq_type", "Start", "End", "Strand")
positions_by_cell <-rbind(dom_pos, cast_pos)

# Visualizing the peaks 
# Plotting genes and coordinates #####
png("plots/Differential_peaks_DOM_CAST.png", width = 1500, height = 500, res = 200)
ggplot(positions_by_cell, aes(xmin = Start, xmax = End, y= seq_type)) +
  geom_gene_arrow() +
  facet_grid(seq_type ~.)+
  theme(
    # Removes the text labels on the Y axis
    axis.text.y = element_blank(),
    # Removes the little tick marks on the Y axis
    axis.ticks.y = element_blank(),
    # Removes the labels on the side/top of the facets (strips)
    #strip.text.y = element_blank()
  )+
  scale_x_continuous(breaks = seq(0, max(positions_by_cell$End)+10000000, 9005000)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Differential Peaks In Dom And Cast WT and SMD5 Mice cells (GSE249505)")

dev.off()

# Reading in RNA-seq differential locus #####
chrom_X_diff_locus_RNAseq <- read_table("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE185843_Analysis/XChr_genomic_sig_DOM_CAST")

# Cleaning and converting to numerical #####
Strand <- chrom_X_diff_locus_RNAseq %>% separate_wider_delim(Strand, delim = ";", names_sep = "",
                                                         too_few="align_start", names_repair = "universal" )  %>%
                                                  select(Strand1)
# Reformatting the strand column #####
chrom_X_diff_locus_RNAseq[,"Strand"] <- Strand

# Renaming the position columns #####
position_cols <- c("Start", "End")

# Converting to numeric #####
chrom_X_diff_locus_RNAseq[, position_cols] <-
  lapply(chrom_X_diff_locus_RNAseq[, position_cols], function(x) {
    as.numeric(gsub('"', '', x))
  })

# Subsetting to join RNA-seq differntial locus and ATAC-seq differential peaks for comparison 
RNA_seq_differential_df <- chrom_X_diff_locus_RNAseq[, c(1,3:5)]

# Creating a column with the DomvsCast_Chrseq seq type #####
RNA_seq_differential_df <- data.frame(rep("DomVsCast_Chrseq", each = nrow(RNA_seq_differential_df)), RNA_seq_differential_df)

# Renaming cols DomvsCast_Chrseq and binding 
names(RNA_seq_differential_df) <- c("seq_type", "Geneid", "Start", "End", "Strand")

all_positions_by_type <- data.frame(rbind(positions_by_cell, RNA_seq_differential_df[,c(1,3,4,5)]))

# Plotting by seq type 
png("plots/Positions_by_seq_type.png", width = 2000, height = 500, res = 200)
ggplot(all_positions_by_type, aes(xmin = Start, xmax = End, y= seq_type)) +
  geom_gene_arrow() +
  facet_grid(seq_type ~.)+
  theme(
    # Removes the text labels on the Y axis
    axis.text.y = element_blank(),
    # Removes the little tick marks on the Y axis
    axis.ticks.y = element_blank(),
    # Removes the labels on the side/top of the facets (strips)
    #strip.text.y = element_blank()
  )+
  scale_x_continuous(breaks = seq(0, max(positions_by_cell$End)+10000000, 9005000)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.y= element_text(size=3))+
  ggtitle("Euchromatic Differential Positions in ATAC-seq and Differential chr-RNA-seq")

dev.off()

# Making a copy of df by seq type #####
all_positions_by_type_non_numeric <- all_positions_by_type

# Formatting the seq types to numeric #####
to_numeric_function <- function(x){
  if( x== "Cast_ATAC_seq"){
      return(1)
  }
  else if( x== "Dom_ATAC_seq"){
      return(2)
  }
   else if( x== "DomVsCast_Chrseq"){
      return(3)
   }
    else{
      return( "NA")
    }
  }

all_positions_by_type$seq_type <- sapply(all_positions_by_type$seq_type, to_numeric_function)

# Converting to numeric ####
all_positions_by_type[, 1:3]<-
  lapply(all_positions_by_type[, 1:3], function(x) {
    as.numeric(x)
  })

all_positions_by_type <- as.data.frame(all_positions_by_type)

# Renaming cols to make a granges object 
names(all_positions_by_type) <- c("seqnames", "start", "end", "strand")

granges_positions <-makeGRangesFromDataFrame(all_positions_by_type, 
                                             start.field = "start",
                                             end.field = "end",
                                            keep.extra.columns = TRUE)

# Subsetting the granges object ####
gr1 <- GRanges(seqnames="chrX", ranges=IRanges(start=all_positions_by_type$start[all_positions_by_type$seqnames==1],
                                               end=all_positions_by_type$end[all_positions_by_type$seqnames==1]), all_positions_by_type$strand[all_positions_by_type$seqnames==1])

gr2 <- GRanges(seqnames="chrX", ranges=IRanges(start=all_positions_by_type$start[all_positions_by_type$seqnames==2],
                                               end=all_positions_by_type$end[all_positions_by_type$seqnames==2]), all_positions_by_type$strand[all_positions_by_type$seqnames==2])

gr3 <- GRanges(seqnames="chrX", ranges=IRanges(start=all_positions_by_type$start[all_positions_by_type$seqnames==3],
                                               end=all_positions_by_type$end[all_positions_by_type$seqnames==3]), all_positions_by_type$strand[all_positions_by_type$seqnames==3])

# Finding overlaps #####
overlaps_gr1_gr2 <- findOverlaps(gr1, gr2, select = "all", type= "any", ignore.strand=TRUE, maxgap = 0L)

overlaps_gr3_gr2 <- findOverlaps(gr3, gr2, select = "all", type = "any", ignore.strand=TRUE, maxgap = 100000L)
overlaps_gr3_gr1 <- findOverlaps(gr3, gr1, select = "all", type = "any", ignore.strand=TRUE, maxgap = 100000L)

# Retrieving overlaps #####
overlap_indices_gr1_gr2 <- overlaps_gr1_gr2[which(!is.na(overlaps_gr1_gr2))]
overlap_indices_gr3_gr2 <- overlaps_gr3_gr2[which(!is.na(overlaps_gr3_gr2))]
overlap_indices_gr3_gr1 <- overlaps_gr3_gr1[which(!is.na(overlaps_gr3_gr1))]

# Overlaps between Cast(1) and Chr_RNA_seq(3) #####
hits_gr3_gr1_df <- as.data.frame(overlap_indices_gr3_gr1)
gr3_gr1_overlap_RNA_seq_differential_df <- unique(chrom_X_diff_locus_RNAseq[hits_gr3_gr1_df$queryHits, ])

# Plotting overlaps Cast and chr_RNA_seq #####
png("plots/Positions_gr3_gr1.png", width = 2000, height = 500, res = 200)
ggplot( gr3_gr1_overlap_RNA_seq_differential_df , aes(xmin = Start, xmax = End, y=Chr)) +
  geom_gene_arrow() +
  theme(
    # Removes the text labels on the Y axis
    axis.text.y = element_blank(),
    # Removes the little tick marks on the Y axis
    axis.ticks.y = element_blank(),
    # Removes the labels on the side/top of the facets (strips)
    #strip.text.y = element_blank()
  )+
  scale_x_continuous(breaks = seq(0, max(positions_by_cell$End)+10000000, 9005000)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.y= element_text(size=3))+
  ggtitle("Overlap Locus in chrRNA-seq(DOM vs CAST NPC cells) and ATAC-Seq (Cast NPC cells)")
dev.off()


gr3_gr1_overlap_df_with_seq<- data.frame(rep("CastATAC_RNAseq_Overlap", each = nrow(gr3_gr1_overlap_RNA_seq_differential_df)), gr3_gr1_overlap_RNA_seq_differential_df[,3:5])
names(gr3_gr1_overlap_df_with_seq) <- names(all_positions_by_type_non_numeric)

all_positions_by_type_non_numeric <- rbind(all_positions_by_type_non_numeric, gr3_gr1_overlap_df_with_seq)
seq_order <- factor(c("Dom_ATAC_seq", "Cast_ATAC_seq", "DomVsCast_Chrseq", "CastATAC_RNAseq_Overlap"))


# Order df #####
all_positions_by_type_non_numeric$seq_type <- factor(
  all_positions_by_type_non_numeric$seq_type, 
  levels = seq_order)


# Plotting all seq types and overlaps 
png("plots/All_seq_types_overlaps.png", width = 2500, height = 1500, res = 100)
ggplot( all_positions_by_type_non_numeric, aes(xmin = Start, xmax = End, y=seq_type)) +
  geom_gene_arrow() +
  facet_grid(seq_type ~.)+
  theme(
    # Removes the text labels on the Y axis
    axis.text.y = element_blank(),
    # Removes the little tick marks on the Y axis
    axis.ticks.y = element_blank(),
    # Removes the labels on the side/top of the facets (strips)
    #strip.text.y = element_blank()
  )+
  scale_x_continuous(breaks = seq(0, max(all_positions_by_type_non_numeric$End)+10000000, 9005000)) +
  theme(axis.text.x = element_text(size=16,angle = 90, hjust = 1)) +
  theme(strip.text.y= element_text(size=16))+
  ggtitle("Overlap Locus in chrRNA-seq(DOM vs CAST NPC cells) and ATAC-Seq (Cast NPC cells), max gap =100000bp")
dev.off()

# Subset to find percent overlap #####
DomVsCast_Chrseq_count <- dim(subset(all_positions_by_type_non_numeric,seq_type== "DomVsCast_Chrseq"))
DomVsCast_ATAC_seq_overlap_count<- dim(subset(all_positions_by_type_non_numeric,seq_type== "CastATAC_RNAseq_Overlap"))

percent_overlap <- DomVsCast_ATAC_seq_overlap_count[1]/DomVsCast_Chrseq_count[1]

# Rename geneid column to retrieve entrez id ######
names(gr3_gr1_overlap_RNA_seq_differential_df)[1] <- "SYMBOL"

# Retrieve Entrez id ######
gr3_gr1_overlap_entrez <- bitr(gr3_gr1_overlap_RNA_seq_differential_df$SYMBOL, 
                               fromType = "SYMBOL",
                               toType = "ENTREZID", 
                               OrgDb = org.Mm.eg.db)

# Read in statistics for chr-rna-seq results
XChr_chrRNA_seq_sig_results <-read.table("/Users/luzsmac/Desktop/5332_Stats_Analysis_of_Genomic_Data /Project/GSE185843_Analysis/XChr_deseq2_results")

# Subset for the gr3 gre 1 overlaps 
XChr_chrRNA_seq_sig_results_overlap_gr1_gr3 <- XChr_chrRNA_seq_sig_results[gr3_gr1_overlap_RNA_seq_differential_df$SYMBOL,]

XChr_chrRNA_seq_sig_results_overlap_gr1_gr3 <- XChr_chrRNA_seq_sig_results_overlap_gr1_gr3[gr3_gr1_overlap_entrez$SYMBOL,]
XChr_chrRNA_seq_sig_results_overlap_gr1_gr3 <- cbind(XChr_chrRNA_seq_sig_results_overlap_gr1_gr3,gr3_gr1_overlap_entrez$ENTREZID)
names(XChr_chrRNA_seq_sig_results_overlap_gr1_gr3)[7] <- "ENTREZID"
gene_entrex_log_fold <- XChr_chrRNA_seq_sig_results_overlap_gr1_gr3[,"log2FoldChange"]
names(gene_entrex_log_fold) <- XChr_chrRNA_seq_sig_results_overlap_gr1_gr3[,"ENTREZID"]

## read in the differential analysis log2 fold statistics from chrrnaseq of the overlaps gr3 and gr1 
gsea_do <- gseDO(
  geneList = sort(gene_entrex_log_fold, decreasing=TRUE),
  organism = "mmu",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  use_internal_data = TRUE
)

go_enrichment_cluster_2 <- enrichGO(
  gene = gr3_gr1_overlap_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",              
  pAdjustMethod = "BH",    
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
test<- as.data.frame(go_enrichment_cluster_2)