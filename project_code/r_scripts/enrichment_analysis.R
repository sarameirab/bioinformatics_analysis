run_enrichment_analysis_on_drug <- function(drug_name, sorted_random_forest_results, gene_slice) {
  GSK_df <- sorted_random_forest_results[c(1, which(sorted_random_forest_results$drug == drug_name)), ]

  GSK_df_subsetted <- GSK_df[, -c(1:4)]
  transposed_GSK_df_subsetted <- t(GSK_df_subsetted) %>% as.data.frame() 
  transposed_GSK_df_subsetted[[2]] <- as.numeric(transposed_GSK_df_subsetted[[2]])
  # Sort rows by the second column (column index 2) in descending order
  transposed_GSK_df_subsetted_sorted <- transposed_GSK_df_subsetted[order(transposed_GSK_df_subsetted[[2]], decreasing = TRUE), ]
  top_gene_names <- transposed_GSK_df_subsetted_sorted[1:gene_slice, 1]
  result <- enrichGO(gene = top_gene_names, 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "BP", 
                     pvalueCutoff = 0.2)
  result_df <- as.data.frame(result)
  return(result_df)
}