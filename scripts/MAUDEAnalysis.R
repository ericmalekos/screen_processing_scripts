#!/usr/bin/env Rscript

# Activate correct conda environment in terminal:
# conda activate maude_r3.6

library(MAUDE)

base_dir <- "/home/eric/Projects/Screen/BMDM_CRISPRi_Screens/MAUDE"

# Find all .tsv files in "Counts" subfolders
all_tsv_files <- list.files(path = base_dir, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
counts_tsv_files <- all_tsv_files[grep("/Counts/", all_tsv_files)]

#print(counts_tsv_files)


options(width = 1000)

# Print the list of .tsv files
for (file_path in counts_tsv_files) {
  # Print the file path
  cat("File:", file_path, "\n")
  binReadMat <- read.table(file = file_path, sep = '\t', header = TRUE)
  
  print(head(binReadMat,2))
  
  binStats <- read.table(file = '/home/eric/Projects/Screen/BMDM_CRISPRi_Screens/MAUDE/binStats.txt', sep = '\t', header = TRUE)
  
  expt_string <- as.character(binReadMat$expt[1])

  binStats$expt = c(rep(expt_string,2))
  
  print(head(binStats,2))

  path_components <- strsplit(file_path, "/")[[1]]
  reduced_path_components <- path_components[-(length(path_components)-(0:1))]
  reduced_path <- paste(reduced_path_components, collapse = "/")
  print(reduced_path)
  
  file_name <- basename(file_path)
  file_name_without_extension <- tools::file_path_sans_ext(file_name)
  print(file_name_without_extension)
  
  output <- paste0(reduced_path,"/",file_name_without_extension)
  print(output)
  
  guideLevelStats = findGuideHitsAllScreens(experiments = unique(binReadMat["expt"]),
                                            countDataFrame = binReadMat, binStats = binStats,
                                            sortBins = c("low", "hi"), unsortedBin = "background",
                                            negativeControl = "isNontargeting")
  
  print(head(guideLevelStats ,2))
  
  elementLevelStats <- getElementwiseStats(experiments = unique(binReadMat["expt"]),
                                           normNBSummaries = guideLevelStats,
                                           negativeControl = "isNontargeting",
                                           elementIDs=c("guide", "name"),tails="both")
  
  print(head(elementLevelStats ,2))
  
  write.table(elementLevelStats[order(elementLevelStats$FDR), ] , file = paste0(output,"_guideLevelStats.tsv"), quote = FALSE, 
              sep = "\t", row.names = FALSE) 
  
  geneLevelStats = getElementwiseStats(experiments = unique(binReadMat["expt"]),
                                           normNBSummaries = guideLevelStats,
                                           negativeControl = "isNontargeting",
                                           elementIDs=c("name"),tails="both")
  
  print(head(geneLevelStats ,2))
  
  write.table(geneLevelStats[order(geneLevelStats$FDR), ] , file = paste0(output,"_geneLevelStats.tsv"), quote = FALSE, 
              sep = "\t", row.names = FALSE) 
  
  cat("\n\n------------------------------------------------------------------------------------------------\n\n\n")

}









