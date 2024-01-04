library(DESeq2)
library(ggpubr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(reshape)

run_deseq_analysis <- function(file_path, condition) {
  # Read count data from the file path
  counts <- read.csv(file_path, sep = "\t", row.names = 1)
  counts <- counts[, !(names(counts) %in% "Gene")]
  
  # Define sample information and experimental design
  samples <- colnames(counts)
  coldata <- data.frame(samples, condition)
  rownames(coldata) <- samples
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ condition)
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Use the second entry of resultsNames for lfcShrink
  coef_name <- resultsNames(dds)[2]
  res_apeglm <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  # Return the results
  return(res_apeglm)
}

process_counts <- function(filepath, sep=",", keep_str = NULL, remove_str = NULL, colData = NULL) {
  counts <- read.csv(filepath, sep = sep, row.names = 1)
  guide_names <- rownames(counts)
  
  counts <- counts[, !(names(counts) %in% "Gene")]

  counts <- data.frame(lapply(counts, as.numeric))
  

  # Preserve row names for later use
  
  column_sums <- colSums(counts, na.rm = TRUE)
  
  # Identify columns to drop
  columns_to_drop <- names(column_sums)[column_sums < 1000000]
  columns_to_drop_string <- paste(columns_to_drop, collapse = ", ")
  print(paste0("Dropping columns with counts < 10^6: ", columns_to_drop_string))
  
  # Keep only columns where the sum is >= 100,000
  counts <- counts[, column_sums >= 1000000]
  
  # Apply keep and remove filters if specified
  if (!is.null(keep_str)) {
    keep_indices <- grep(keep_str, names(counts))
    counts <- counts[, keep_indices]
  }
  if (!is.null(remove_str)) {
    remove_indices <- grep(remove_str, names(counts))
    counts <- counts[, -remove_indices]
  }
  
  # Reassigning the preserved gene names
  rownames(counts) <- guide_names
  
  # Create colData if not provided
  if (is.null(colData)) {
    colData <- data.frame(condition = colnames(counts))
    rownames(colData) <- colnames(counts)
  }
  
  # Create a DESeqDataSet object, ensuring row names are preserved
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ 1)
  
  count_matrix <- counts(dds)
  count_df <- as.data.frame(count_matrix)
  count_df$Guide <- rownames(guide_names)
  
  # Minimum number of samples and count threshold
  n <- 3  
  threshold <- 10
  
  # Filter out genes with low counts
  keep <- rowSums(counts(dds) >= threshold) >= n
  dds <- dds[keep,]
  
  rlog_counts <- rlog(dds, blind = TRUE)
  
  pca_plot <- plotPCA(rlog_counts, intgroup = c("condition")) +
    geom_text(aes(label = name), vjust = -0.5, nudge_x = -1) + 
    ggtitle("PCA of DESeq2 normalized log transform (rld function)") + 
    theme(legend.position = "bottom", legend.box = "horizontal") +
    labs(title = "PCA of DESeq2 normalized log transform (rld function)", color = "Samples")
  
  print(pca_plot)
  
  sampleDists <- dist(t(assay(rlog_counts)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rlog_counts$condition, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, main = "Sample Distance Matrix")
  
  return(list(dds = dds, count_df = count_df, rlog_counts = rlog_counts))
}

calculate_differences <- function(rlog_object, gene_name, prefix, high_suffix = "GFPhi", low_suffix = "GFPLow") {
  # Extract rows corresponding to the gene
  gene_rows <- assay(rlog_object[grep(gene_name, rownames(rlog_object)), ])
  
  # Convert to data frame if it's not already
  if (is.matrix(gene_rows)) {
    gene_rows <- as.data.frame(gene_rows)
  }
  
  # Check if any rows are found
  if (nrow(gene_rows) == 0) {
    cat("No rows found for gene:", gene_name, "\n")
    return(NULL)
  }
  
  # Calculate differences where both columns exist
  pairs <- c("A", "B", "C")  # Define pairs
  for (pair in pairs) {
    high_col <- paste0(prefix, pair, "_", high_suffix)
    low_col <- paste0(prefix, pair, "_", low_suffix)
    
    if (high_col %in% names(gene_rows) && low_col %in% names(gene_rows)) {
      diff_col <- paste0(prefix, pair, "_Diff")
      gene_rows[[diff_col]] <- gene_rows[[high_col]] - gene_rows[[low_col]]
    }
  }
  
  avg_row <- colMeans(gene_rows, na.rm = TRUE)
  gene_rows <- rbind(gene_rows, avg_row)
  rownames(gene_rows)[nrow(gene_rows)] <- "Average"
  
  # Return the modified data frame
  print(gene_rows)
}

process_genes <- function(dataframe, gene_names, prefix) {
  for (gene in gene_names) {
    cat("\nResults for", gene, ":\n")
    calculate_differences(dataframe, gene, prefix)
  }
}

# Example usage
# Assuming df1, df2, etc. are your dataframes with count data
result_df <- process_dataframes(list(df1, df2, df3), "substring", add_average = TRUE)

setwd("/home/eric/Projects/Screen/BMDM_CRISPRi_Screens")
davis_file_path = "./NFKB/iBMDM_Davis/ScreenCounts.count.txt"
qb3_file_path = "./NFKB/iBMDM_QB3/ScreenCounts.count.txt"

pos_ctrls = c("Myd88", "Tlr4", "Rela","Traf","mP3","gfp")
#process_counts(davis_file_path, keep_str = "NT|LPS", remove_str = "CS")
davis_result = process_counts(davis_file_path, keep_str = "LPS", remove_str = "CS"); davis_dds <- davis_result$ddsdavis_count ; davis_rlog <- davis_result$rlog_counts
qb3_result = process_counts(qb3_file_path, keep_str = "LPS", remove_str = "CS"); qb3_dds <- qb3_result$dds ; qb3_rlog <- qb3_result$rlog_counts

prefix = "LPS_"

process_genes(davis_rlog, pos_ctrls, prefix); process_genes(qb3_rlog, pos_ctrls, prefix)


davis_result = process_counts(davis_file_path, keep_str = "CS_LPS") ; davis_dds <- davis_result$ddsdavis_count ; davis_rlog <- davis_result$rlog_counts
qb3_result = process_counts(qb3_file_path, keep_str = "CS_LPS") ; qb3_dds <- qb3_result$dds ; qb3_rlog <- qb3_result$rlog_counts

prefix = "CS_LPS_"

process_genes(davis_rlog, pos_ctrls, prefix); process_genes(qb3_rlog, pos_ctrls, prefix)

davis_result = process_counts(davis_file_path, keep_str = "CS", remove_str = "LPS") ; davis_dds <- davis_result$ddsdavis_count ; davis_rlog <- davis_result$rlog_counts
qb3_result = process_counts(qb3_file_path, keep_str = "CS", remove_str = "LPS") ; qb3_dds <- qb3_result$dds ; qb3_rlog <- qb3_result$rlog_counts

prefix = "CS_"

process_genes(davis_rlog, pos_ctrls, prefix); process_genes(qb3_rlog, pos_ctrls, prefix)



##################################################
##################################################
### 
###  Preparing counts for MAUDE
###
##################################################
##################################################

combine_and_normalize <- function(file_paths, substring) {
  df_list <- lapply(file_paths, function(file) {
    df <- read.csv(file, sep = ",", row.names = 1)
    print(paste("Dimensions of", file, "before conversion:", dim(df)))
    
    # Apply numeric conversion column-wise
    for (col in colnames(df)) {
      df[[col]] <- as.numeric(as.character(df[[col]]))
    }
    
    print(paste("Dimensions of", file, "after conversion:", dim(df)))
    return(df)
  })
  
  # Combine dataframes
  if (length(df_list) > 1) {
    combined_df <- Reduce(`+`, df_list)
    print(paste("Dimensions of combined dataframe:", dim(combined_df)))
  } else {
    combined_df <- df_list[[1]]
    print(paste("Dimensions of single dataframe:", dim(combined_df)))
  }
  
  # Subset combined_df
  subset_combined_df <- combined_df[, grep(substring, colnames(combined_df))]
  head(dim(subset_combined_df))
  head(subset_combined_df)
  
  
  colData <- data.frame(condition = colnames(subset_combined_df))
  rownames(colData) <- colnames(subset_combined_df)
  dds <- DESeqDataSetFromMatrix(countData = subset_combined_df, colData = colData, design = ~ 1)
  
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  
  if (!is.data.frame(norm_counts)) {
    norm_counts <- as.data.frame(norm_counts)
  }
  norm_counts$guide <- rownames(norm_counts)
  colnames(norm_counts) <- tolower(colnames(norm_counts))
  print(head(norm_counts,3))
  
  return(norm_counts)
}

counts_to_maude <- function(count_input, guide_lib) {
  counttable <- merge(guide_lib, count_input, by = "guide")
  counttable$isNontargeting = grepl("non-targeting", counttable$symbol);head(counttable,2);tail(counttable,2)
  
  counttable = melt(counttable, id.vars = c("isNontargeting","guide","symbol","name"));head(counttable,2);tail(counttable,2)
  counttable$expt <- gsub(".*_(.)_.*", "\\1", counttable$variable)
  counttable$Bin = gsub("(.*)(hi|lo|low|high)","\\2", counttable$variable);head(counttable,2);tail(counttable,2)
  counttable$value= floor(as.numeric(counttable$value)); head(counttable,2);tail(counttable,2)
  binReadMat = data.frame(cast(counttable, 
                               guide+symbol+name+isNontargeting+expt ~ Bin, value="reads"))
  head(binReadMat,2);tail(binReadMat,2)
  binReadMat$background  <- apply(binReadMat[,6:7],1, function(x) floor(mean(x)));head(binReadMat,2);tail(binReadMat,2)
  
  return(binReadMat)
}
  

# Need to merge with Gene on sgRNA
# add isNonTargeting Column

guideLib <- "./guidelib/lib8_allGFP.tsv.txt"
guideDF <- read.csv(guideLib, sep = "\t", header=FALSE) ; names(guideDF) <- c("guide","sequence","symbol")
symbolToGeneName <- "./guidelib/lib8_trunc_cap_uniqseq_column3ToGeneName_v3.xref"
geneNameDF <- read.csv(symbolToGeneName, sep = "\t", header=FALSE) ; names(geneNameDF) <- c("symbol","name")

#####################
# Want a df like this:
#
# guide symbol  name
# 50732493.23-P1P2 ENSMUSG00000022639               Dubr
# 50732382.23-P1P2 ENSMUSG00000022639               Dubr
# 96082692.23-P1P2 ENSMUSG00000022915 ENSMUSG00000022915
# 96082689.23-P1P2 ENSMUSG00000022915 ENSMUSG00000022915
#####################

# Merging guideDF and geneNameDF on 'symbol'
merged_df <- merge(guideDF, geneNameDF, by = "symbol")
# Selecting only the desired columns: 'guide', 'symbol', and 'name'
final_df <- merged_df[c("guide", "symbol", "name")]
# For some Gfp entries are lost:
final_df[rowSums(apply(merged_df, 2, function(x) grepl("Gfp", x, ignore.case = TRUE))) > 0, ]
# Adding GFP rows to the dataframe
gfp_rows <- data.frame(guide=c("gfp_em01","gfp_sc01"), symbol=c("Gfp_em01","Gfp_sc01"),name=c("Gfp_em01","Gfp_sc01"))
final_df <- rbind(final_df, gfp_rows)
final_df[rowSums(apply(final_df, 2, function(x) grepl("Gfp", x, ignore.case = TRUE))) > 0, ]


##################################################
##################################################
### 
###  Preparing matrices for MAUDE
###
##################################################
##################################################

# Save dataframes for MAUDE
davis_path <- c(davis_file_path)
qb3_path <- c(qb3_file_path)
file_paths <- c(davis_file_path, qb3_file_path)


experiments <-c("LPS_A","LPS_B","LPS_C",
                "CS_LPS_A","CS_LPS_B",
                "CS_A","CS_B","CS_C")
for (exp in experiments) {
  outdir<-paste0(getwd(),"/MAUDE/",exp,"/Counts/")
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  print(exp)
  normalized_qb3 <- combine_and_normalize(file_paths = qb3_file_path,
                                          substring = paste0("^",exp,"_GFP"))
  normalized_davis <- combine_and_normalize(file_paths = davis_file_path,
                                              substring = paste0("^",exp,"_GFP"))
  normalized_combined <- combine_and_normalize(file_paths = file_paths,
                                    substring = paste0("^",exp,"_GFP"))
  
  binreadmat_qb3 <- counts_to_maude(count_input=normalized_qb3, guide_lib=final_df)
  binreadmat_davis <- counts_to_maude(count_input=normalized_davis, guide_lib=final_df)
  binreadmat_combined <- counts_to_maude(count_input=normalized_combined, guide_lib=final_df)
  
  write.table(binreadmat_qb3, file = paste0(outdir,"/qb3_reads.tsv"), quote = FALSE,
              sep = "\t", row.names = FALSE)
              
  write.table(binreadmat_davis, file = paste0(outdir,"/davis_reads.tsv"), quote = FALSE, 
              sep = "\t", row.names = FALSE)
              
  write.table(binreadmat_combined, file = paste0(outdir,"/combined_reads.tsv"), quote = FALSE, 
              sep = "\t", row.names = FALSE)
  
}



## binStats matrix, like this:
# Bin binStartQ binEndQ fraction  binStartZ    binEndZ expt
# 1  lo     0.001   0.200    0.199 -3.0902323 -0.8416212    A
# 2  hi     0.799   0.999    0.200  0.8380547  3.0902323    A
# 3  lo     0.001   0.200    0.199 -3.0902323 -0.8416212    B
# 4  hi     0.799   0.999    0.200  0.8380547  3.0902323    B

make_bin_stat_matrix <- function(binStartQ=c(0.001,0.799), binEndQ=c(0.20, 0.999), Bin = c('lo','hi')) {
  
  binStats <- data.frame(Bin,binStartQ,binEndQ)
  binStats$fraction = binStats$binEndQ - binStats$binStartQ
  binStats$binStartZ = qnorm(binStats$binStartQ)
  binStats$binEndZ = qnorm(binStats$binEndQ)
  
  return(binStats)
}

binStats <- make_bin_stat_matrix(binStartQ=c(0.001,0.799), binEndQ=c(0.20, 0.999), Bin = c('low','hi'))

write.table(binStats, file = "./MAUDE/binStats.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE) 

#binStats = rbind(binStats, binStats) #duplicate data
#binStats$expt = c(rep("A",2),rep("B",2)) 

