library("googledrive")
source("~/workspace/bioinformatics_analysis/project_code/r_scripts/download_data.R")
library("digest")
library("dplyr")
library("clusterProfiler")
library("org.Hs.eg.db")
source("~/workspace/bioinformatics_analysis/project_code/r_scripts/enrichment_analysis.R")

GENE_SLICE <- 1000
NUM_BEST_DRUGS <- 50



output_file_path <-download_url(RANDOM_FOREST_RESULTS_URL)
random_forest_results <- read.csv(output_file_path)
sorted_random_forest_results <- random_forest_results %>%
  arrange(desc(pearson_correlation))
gene_symbols <- names(sorted_random_forest_results)[-c(1:4)]
gene_symbols <- sub("\\.\\..*$", "", gene_symbols)  # remove everything after and including '..'

# Create a new row where the first 4 columns are NA and the rest are gene symbols
new_row <- c(rep(NA, 4), gene_symbols)

# Bind the new row at the top of the dataframe
sorted_random_forest_results <- rbind(new_row, sorted_random_forest_results)


random_forest_results_50_best_drugs <- sorted_random_forest_results %>%
  slice_head(n = NUM_BEST_DRUGS)

for (drug_name in random_forest_results_50_best_drugs$drug) {
  if (!is.na(drug_name) && !drug_name == "") {
    message("Running enrichment analysis on drug: ", drug_name)
    GSK_df <- random_forest_results_50_best_drugs[c(1, which(random_forest_results_50_best_drugs$drug == drug_name)), ]
  
    GSK_df_subsetted <- GSK_df[, -c(1:4)]
    transposed_GSK_df_subsetted <- t(GSK_df_subsetted) %>% as.data.frame() 
    transposed_GSK_df_subsetted[[2]] <- as.numeric(transposed_GSK_df_subsetted[[2]])
    # Sort rows by the second column (column index 2) in descending order
    transposed_GSK_df_subsetted_sorted <- transposed_GSK_df_subsetted[order(transposed_GSK_df_subsetted[[2]], decreasing = TRUE), ]
    top_gene_names <- transposed_GSK_df_subsetted_sorted[1:GENE_SLICE, 1]
    
    result <- enrichGO(gene = top_gene_names, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "SYMBOL", 
                       ont = "BP", 
                       pvalueCutoff = 0.2)
    
    result_df <- as.data.frame(result)
    results_df <- run_enrichment_analysis_on_drug(drug_name, sorted_random_forest_results, GENE_SLICE)
    
    # Clean the drug name
    clean_name <- gsub("[^A-Za-z0-9]", "_", drug_name)
    
    # Create full file path
    output_dir <- "~/workspace/bioinformatics_analysis/data"
    full_path <- file.path(output_dir, paste0(clean_name, ".csv"))
    message("Writing results on drug: ", drug_name)
    # Write the data frame
    write.csv(results_df, full_path, row.names = FALSE)
  }
}




