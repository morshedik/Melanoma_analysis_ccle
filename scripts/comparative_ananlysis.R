##############Investigate correlations between gene expression and drug response############
BiocManager::install("biomaRt", force = TRUE)
library(biomaRt)

# Set up biomaRt
ensembl <- tryCatch({
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")
}, error = function(e) {
  message("Error connecting to Ensembl. Using a local approach instead.")
  NULL
})

# Function to get gene synonyms
get_synonyms <- function(gene_name) {
  if (!is.null(ensembl)) {
    results <- getBM(attributes = c("hgnc_symbol", "external_synonym"), 
                     filters = "external_synonym", 
                     values = gene_name, 
                     mart = ensembl)
    unique(c(results$hgnc_symbol, results$external_synonym))
  } else {
    # Fallback: return the gene name itself if biomaRt is unavailable
    return(gene_name)
  }
}

# Create a dictionary of gene synonyms for PUTATIVE_TARGETs
putative_targets <- unique(drug_response$PUTATIVE_TARGET)
gene_dict <- setNames(lapply(putative_targets, get_synonyms), putative_targets)

# Print the first few entries of gene_dict to verify
print(head(gene_dict))

# Create a dictionary of gene synonyms for PUTATIVE_TARGETs
putative_targets <- unique(drug_response$PUTATIVE_TARGET)
gene_dict <- setNames(lapply(putative_targets, get_synonyms), putative_targets)


# Select top 10 effective drugs
top_drugs <- head(drug_response, 15)

# Function to calculate correlation
calc_correlation <- function(drug_name, target) {
  drug_ic50 <- drug_response$LN_IC50[drug_response$DRUG_NAME == drug_name]
  target_expression <- gene_expression$rna_expression[gene_expression$gene_name == target]
  
  if (length(target_expression) > 0 && length(drug_ic50) > 0) {
    mean_expression <- mean(target_expression, na.rm = TRUE)
    return(data.frame(Drug = drug_name, Target = target, 
                      Target_Expression = mean_expression, 
                      Drug_IC50 = mean(drug_ic50, na.rm = TRUE)))
  } else {
    return(data.frame(Drug = drug_name, Target = target, 
                      Target_Expression = NA, Drug_IC50 = NA))
  }
}

print(head(drug_response))
print(head(gene_expression))

# Calculate correlations
correlations <- map2_dfr(top_drugs$DRUG_NAME, top_drugs$PUTATIVE_TARGET, calc_correlation)

# Remove rows with NA values
correlations <- correlations %>% filter(!is.na(Target_Expression) & !is.na(Drug_IC50))

# Print correlations
print("Correlations between drug response and target gene expression:")
print(correlations)

# Calculate overall correlation if there's enough data
if (nrow(correlations) > 2) {
  cor_result <- cor.test(correlations$Target_Expression, correlations$Drug_IC50)
  print(paste("Overall correlation between target gene expression and drug IC50:", 
              round(cor_result$estimate, 3)))
  print(paste("P-value:", format.pval(cor_result$p.value, digits = 3)))
} else {
  print("Not enough data to calculate overall correlation")
}

# Save correlations
write_csv(correlations, "~/Desktop/R_output/Final_data/a375_drug_gene_correlations.csv")


install.packages('pheatmap')
library(ggplot2)
library(pheatmap)

# Create a more informative scatterplot
ggplot(correlations, aes(x = Target_Expression, y = Drug_IC50)) +
  geom_point(aes(color = Drug), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal() +
  labs(title = "Correlation between Target Gene Expression and Drug IC50",
       x = "Target Gene Expression", y = "Drug IC50 (log)")




######Comparison of Pathways
install.packages("rstatix")
library(rstatix)

# Prepare data for analysis
pathway_comparison_data <- prepared_go_results %>%
  select(cell_line, ID, Description, p.adjust) %>%
  group_by(ID, Description) %>%
  filter(n() > 1) %>%  # Keep only pathways present in more than one cell line
  ungroup()

# Perform Kruskal-Wallis test for each pathway
pathway_differences <- pathway_comparison_data %>%
  group_by(ID, Description) %>%
  kruskal_test(p.adjust ~ cell_line) %>%
  adjust_pvalue(method = "BH") %>%
  arrange(p.adj) %>%
  filter(p.adj < 0.05)  # Keep only significantly different pathways

# Print results
print("Top 10 pathways with significant differences between cell lines:")
print(head(pathway_differences, 10))

# Save results
write_csv(pathway_differences, "differentially_enriched_pathways.csv")
