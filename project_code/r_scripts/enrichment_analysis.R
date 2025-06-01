library(dplyr)
library("clusterProfiler")
library("org.Hs.eg.db")


run_enrichment_analysis_on_drug <- function(drug_name, sorted_random_forest_results, gene_slice) {
  GSK_df <- sorted_random_forest_results[c(1, which(sorted_random_forest_results$drug == drug_name)), ]

  GSK_df_subsetted <- GSK_df[, -c(1:4)]
  transposed_GSK_df_subsetted <- t(GSK_df_subsetted) %>% as.data.frame() 
  transposed_GSK_df_subsetted[[2]] <- as.numeric(transposed_GSK_df_subsetted[[2]])
  # Sort rows by the second column (column index 2) in descending order
  transposed_GSK_df_subsetted_sorted <- transposed_GSK_df_subsetted[order(transposed_GSK_df_subsetted[[2]], decreasing = TRUE), ]
  if (gene_slice == -1) {
    top_gene_names <- transposed_GSK_df_subsetted_sorted[, 1]
  } else {
    top_gene_names <- transposed_GSK_df_subsetted_sorted[1:gene_slice, 1]
  }
  result <- enrichGO(gene = top_gene_names, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.05)  # ← Enable multithreading here
  result_df <- as.data.frame(result)
  return(result_df)
}

prepare_results_data <- function(results_df) {
  sorted_random_forest_results <- results_df %>%
    arrange(desc(pearson_correlation))
  
  gene_symbols <- names(sorted_random_forest_results)[-c(1:4)]
  gene_symbols <- sub("\\.\\..*$", "", gene_symbols)  # remove everything after and including '..'
  
  # Create a new row where the first 4 columns are NA and the rest are gene symbols
  new_row <- c(rep(NA, 4), gene_symbols)
  
  # Bind the new row at the top of the dataframe
  sorted_random_forest_results <- rbind(new_row, sorted_random_forest_results)
  
  return(sorted_random_forest_results)
}

filter_to_include_gene_names_with_score_higher_than <- function(prepared_results_df, drug_name, min_score) {
  one_drug_df <- prepared_results_df %>%filter(drug == drug_name | row_number()==1)
  transposed_one_drug_df <- t(one_drug_df) %>% as.data.frame()
  # Get the row index of the drug of interest
  transposed_one_drug_df_filtered <- transposed_one_drug_df %>%
    filter(V2 > min_score | row_number() < 5)  # Assuming V2 is the column with scores
  # Keep first 4 rows as-is
  # First 4 rows (unchanged)
  transposed_one_drug_df_filtered_top <- transposed_one_drug_df_filtered[1:4, ]
  
  # Remaining rows
  transposed_one_drug_df_filtered_rest <- transposed_one_drug_df_filtered[5:nrow(transposed_one_drug_df_filtered), ]
  
  # Sort by V2 in descending order
  transposed_one_drug_df_filtered_rest_sorted <- transposed_one_drug_df_filtered_rest[order(-as.numeric(transposed_one_drug_df_filtered_rest$V2)), ]
  
  # Combine back
  filtered_results_with_score_higher_than_min_score <- rbind(transposed_one_drug_df_filtered_top, transposed_one_drug_df_filtered_rest_sorted)  
  return(t(filtered_results_with_score_higher_than_min_score)  %>% as.data.frame()) 
}