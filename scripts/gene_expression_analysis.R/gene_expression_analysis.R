# Filter for A375 in CCLE
ccle_a375 <- ccle_expression %>%
  filter(str_detect(cell_line, "A375"))

# Filter for A375 in GDSC
gdsc_a375 <- gdsc_response %>%
  filter(str_detect(CELL_LINE_NAME, "A375"))
print(gdsc_a375)


# Analyze gene expression data for A375
gene_expression_a375 <- ccle_a375 %>%
  arrange(desc(rna_expression)) %>%
  select(gene_name, rna_expression)
print("Top 20 highly expressed genes in A375:")
print(head(gene_expression_a375, 20))



# Analyze drug response data for A375
drug_response_a375 <- gdsc_a375 %>%
  arrange(LN_IC50) %>%
  select(DRUG_NAME, PUTATIVE_TARGET, LN_IC50, AUC)
print("Top 20 most effective drugs on A375 (based on lowest LN_IC50):")
print(head(drug_response_a375, 20))
# Save results
output_direct <- "~/Desktop/R_output/Final_data"
dir.create(output_direct, showWarnings = FALSE)
write_csv(gene_expression_a375, file.path(output_direct,"a375_gene_expression.csv"))
write_csv(drug_response_a375, file.path(output_direct,"a375_drug_response.csv"))

