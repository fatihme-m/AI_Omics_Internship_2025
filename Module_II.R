# Assignment 2
# --------------------------
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
#The analysis produces two key measures for each gene:

library(readr)
DEGs_Data_1 <- read_csv("DEGs_Data_1.csv")
View(DEGs_Data_1)
library(readr)
DEGs_Data_2 <- read_csv("DEGs_Data_2.csv")

# Define input and output folders
input_dir <- "Raw_Data" 
output_dir <- "Results"

# create output folder if not already exist
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}



# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.

# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise


# Then:
#   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
#   - Replace missing padj values with 1
#   - Add a new column 'status'
#   - Save processed files into Results folder
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries

#files to process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv") 


#empty list to store results 
result_list <- list()


# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.
classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1   #replace missing values
  if (logFC > 1 && padj < 0.05) {
    return("Upregulated gene")
  } else if (logFC < -1 && padj < 0.05) {
    return("Downregulated gene")
  } else {
    return("Not_Significant")
  }
}

for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  # Import dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  # Replace missing padj with 1 
  if ("padj" %in% names(data)) {
    missing_count <- sum(is.na(data$padj))
    cat("Missing values in 'padj':", missing_count, "\n")
    data$padj[is.na(data$padj)] <- 1
  }
  
  # Classify genes row by row
  data$Gene_Class <- mapply(classify_gene, data$logFC, data$padj)
  cat("Genes have been classified successfully.\n")
  
  # Save results in list
  result_list[[file_names]] <- data 
  
  # Save results in Results folder
  output_file_path <- file.path(output_dir, paste0("classification_", file_names))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
  
  # Print summary counts
  cat("Summary counts for", file_names, ":\n")
  print(table(data$Gene_Class))
}



# Access processed results 
results_1 <- result_list[[files_to_process[1]]]
results_2 <- result_list[[files_to_process[2]]]

save.image(file = "FatihmeMaarawi_Class_2_Assignment.RData")