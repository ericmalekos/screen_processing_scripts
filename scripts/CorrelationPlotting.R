#!/usr/bin/env Rscript

# Loading necessary libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(tidyr)

process_pair <- function(file1, file2, n, specific_names) {
  data1 <- read.table(file1, header = TRUE, sep = "\t")
  data2 <- read.table(file2, header = TRUE, sep = "\t")
  
  dir1 <- dirname(file1)
  dir2 <- dirname(file2)
  dir1_name <- basename(dir1)
  dir2_name <- basename(dir2)
  
  # Sorting and selecting the top 'n' rows if 'n' is specified
  if (!is.null(n)) {
    data1_top_n <- head(data1[order(-abs(data1$significanceZ)), ], n)
    data2_top_n <- head(data2[order(-abs(data2$significanceZ)), ], n)
  }
  
  # Select rows that match the specific names in the vector
  data1_specific_names <- data1[data1$name %in% specific_names, ]
  data2_specific_names <- data2[data2$name %in% specific_names, ]
  
  # Combine the top 'n' entries with the specific name entries
  combined_top_n <- unique(rbind(data1_top_n, data2_top_n, data1_specific_names, data2_specific_names))
  
  # Merging the combined top 'n' entries with the original data
  merged_data <- combined_top_n %>%
    merge(data1, by = "name", suffixes = c("", "_file1"), all.x = TRUE) %>%
    merge(data2, by = "name", suffixes = c("", "_file2"), all.x = TRUE)
  
  # Removing duplicate columns
  merged_data <- merged_data[!duplicated(names(merged_data))]
  
  # Renaming the significanceZ columns for clarity
  names(merged_data)[names(merged_data) == "significanceZ.x"] <- "file1"
  names(merged_data)[names(merged_data) == "significanceZ.y"] <- "file2"
  
  merged_data <- merged_data %>% 
    group_by(name) %>% 
    summarize(across(starts_with("significanceZ"), mean), .groups = 'drop')
  
  
  # Computing the correlation coefficient
  correlation_coefficient <- cor(merged_data$significanceZ_file1, merged_data$significanceZ_file2,
                                 use = "complete.obs", method = "spearman")
  
  # Create a new column to indicate the panel label
  # file1_label <- basename(sub(".*(CS_[^/]+)_qb3_reads.*", "\\1", file1))
  # file2_label <- basename(sub(".*(CS_[^/]+)_qb3_reads.*", "\\1", file2))
  panel_label <- paste(dir1_name, "vs", dir2_name, "top", n, "genes\nCorr:", round(correlation_coefficient, 3))
  
  merged_data$panel <- panel_label
  merged_data$color <- ifelse(merged_data$name %in% specific_names, "red", "black")
  
  return(merged_data)
}

# List of input files
# input_files <- c("/home/eric/Projects/Screen/BMDM_CRISPRi_Screens/MAUDE/LPS_A/qb3_reads_geneLevelStats.tsv",
#                  "/home/eric/Projects/Screen/BMDM_CRISPRi_Screens/MAUDE/LPS_B/qb3_reads_geneLevelStats.tsv",
#                  "/home/eric/Projects/Screen/BMDM_CRISPRi_Screens/MAUDE/LPS_C/qb3_reads_geneLevelStats.tsv")


args <- commandArgs(trailingOnly = TRUE)

# Assuming that the command line arguments are as follows:
# args[1]: comma-separated list of input files
# args[2]: value for 'n'
# args[3]: output file name

# Parse arguments
input_files <- strsplit(args[1], ",")[[1]]
n <- as.numeric(args[2])
output_file_name <- args[3]

# Define specific names and 'n'
specific_names <- c("Tlr4", "mP3", "Gfp_sc01", "Gfp_em01", "Myd88","Traf6","Zeb2os","Ptgs2os2","Tug1")

# Generating all pairwise combinations and processing
processed_data <- map_dfr(combn(input_files, 2, simplify = FALSE), 
                          ~process_pair(.x[1], .x[2], n, specific_names))

# Plotting with increased vertical space
p <- ggplot(processed_data, aes(x = significanceZ_file1, y = significanceZ_file2)) +
  geom_point() +
  geom_text_repel(aes(label = name, color = color), 
                  box.padding = 0.35, point.padding = 0.5,
                  segment.color = 'grey50', max.overlaps = 40, size = 3) +
  scale_color_identity() +
  facet_wrap(~ panel, scales = "free", ncol = 1) +  # 'free' scales allow for individual axes
  theme_minimal() +
  theme(legend.position = "none", 
        strip.text = element_text(size = 8),  # Smaller facet titles
        panel.spacing.y = unit(1, "lines"),   # Increase spacing between panels
        panel.background = element_rect(fill = "white", colour = "white"),  # Set panel background to white
        plot.background = element_rect(fill = "white", colour = "white"),   # Set plot background to white
        strip.background = element_rect(fill = "white", colour = "white"),  # Set background for facet titles to white
        strip.placement = "outside")  # Place the strip labels outside of the plot area

# View the plot
if(length(input_files) == 2){
  pdf(paste0(output_file_name,".pdf"))
  # Draw the plot on the PDF device
  print(p)
  # Close the PDF device
  dev.off()
}

# Increase the dimensions of the plot to prevent squishing
ggsave(output_file_name, plot = p, width = 12, height = 20, dpi = 300)



