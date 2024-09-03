# Filter for A375 in CCLE
ccle_a375 <- ccle_expression %>%
  filter(str_detect(cell_line, "A375"))



# Analyze gene expression data for A375
gene_expression_a375 <- ccle_a375 %>%
  arrange(desc(rna_expression)) %>%
  select(gene_name, rna_expression)
print("Top 20 highly expressed genes in A375:")
print(head(gene_expression_a375, 20))



# Save results
output_direct <- "~/Desktop/R_output/Final_data"
dir.create(output_direct, showWarnings = FALSE)
write_csv(gene_expression_a375, file.path(output_direct,"a375_gene_expression.csv"))


# Basic statistical summary of gene expression
gene_expression_summary <- summary(ccle_a375$rna_expression)
print("Summary of A375 gene expression:")
print(gene_expression_summary)

# Visualize gene expression distribution
pdf("~/Desktop/R_output/Final_data/a375_gene_expression_histogram.pdf")
hist(ccle_a375$rna_expression, main="Distribution of Gene Expression in A375", xlab="RNA Expression")
dev.off()


