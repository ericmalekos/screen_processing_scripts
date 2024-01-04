#!/usr/bin/env Rscript

# Loading necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Function to process files and create a plot
process_files_and_plot <- function(directory, pattern, special_genes = NULL) {
  # Recursively find all files in the directory with the given suffix pattern
  file_list <- list.files(path = directory, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  # Convert the comma-separated list of special genes into a vector
  special_genes_vector <- unlist(strsplit(special_genes, ","))
  
  # Loop through files and read data
  for (file in file_list) {
    data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Sort the data frame by significanceZ in ascending order and add the numerical rank
    data <- data %>%
      arrange(significanceZ) %>%
      mutate(Rank = row_number(),
             Label = ifelse(name %in% special_genes_vector | significanceZ > 3 | significanceZ < -3, name, NA),
             Color = ifelse(name %in% special_genes_vector, "blue", ifelse(significanceZ > 3 | significanceZ < -3, "red", "black")))
    
    # Extract the filename without the extension for the plot title
    plot_title <- sub(pattern = "\\.tsv$", replacement = "", x = basename(file))
    
    # Create the plot
    p <- ggplot(data, aes(x = Rank, y = significanceZ, label = Label, color = Color)) +
      geom_point() +
      geom_text_repel(data = subset(data, !is.na(Label)),
                      aes(label = Label), box.padding = 0.35, point.padding = 0.5, size = 3) +
      geom_hline(yintercept = -3, linetype = "dashed", color = "blue") +
      geom_hline(yintercept = 3, linetype = "dashed", color = "blue") +
      scale_color_identity() +
      labs(title = plot_title, x = "Rank", y = "Significance Z-score") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white", colour = "white"),
            plot.background = element_rect(fill = "white", colour = NA),
            axis.title = element_text(colour = "black"),
            axis.text = element_text(colour = "black"),
            axis.line = element_line(colour = "black"),
            strip.background = element_rect(fill = "white", colour = "black"),
            strip.text = element_text(colour = "black"),
            legend.position = "none")
    
    # Define the output file name
    output_file_name <- sub(pattern = "\\.tsv$", replacement = ".png", x = file)
    
    # Save the plot as a PNG file
    ggsave(output_file_name, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  }
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assuming that the command line arguments are as follows:
# args[1]: directory path
# args[2]: file pattern suffix
# args[3]: (optional) comma-separated list of special gene names

# Apply the function with command line arguments
special_genes_input <- if(length(args) > 2) args[3] else ""
process_files_and_plot(args[1], args[2], special_genes_input)