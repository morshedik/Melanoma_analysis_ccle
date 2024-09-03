# Melanoma_analysis_ccle
# Melanoma Cell Line Analysis

This repository contains R code for analyzing gene expression and drug response data in melanoma cell lines using the Cancer Cell Line Encyclopedia (CCLE) and Genomics of Drug Sensitivity in Cancer (GDSC) datasets.

## Table of Contents
1. [Setup and Data Loading](#setup-and-data-loading)
2. [Data Processing](#data-processing)
3. [Gene Expression Analysis](#gene-expression-analysis)
4. [Drug Response Analysis](#drug-response-analysis)
5. [Pathway Enrichment Analysis](#pathway-enrichment-analysis)
6. [Comparative Analysis](#comparative-analysis)

## Setup and Data Loading

```R
# Install and load required packages
packages <- c('tidyverse', 'data.table', "devtools", "usethis", "rvest", "stringdist", 'readxl')
install.packages(packages)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("depmap")

# Load libraries
for (pkg in packages) {
  library(pkg, character.only = TRUE)
}
library(depmap)

# Load CCLE expression data
ccle_expression <- depmap_TPM()
write_csv(ccle_expression, "CCLE_expression.csv")

# Load GDSC drug response data
gdsc_url <- 'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx'
download.file(gdsc_url, destfile = "GDSC2_fitted_dose_response_latest.csv")
gdsc_response <- fread("GDSC2_fitted_dose_response_latest.csv")
```

This section sets up the R environment by installing and loading necessary packages. It then loads the CCLE gene expression data and GDSC drug response data.

## Data Processing

```R
# Function to get unique cell lines and their associated information
get_cell_line_info <- function(data, cell_line_col) {
  data %>%
    select(!!sym(cell_line_col)) %>%
    unique() %>%
    arrange(!!sym(cell_line_col))
}

# Get CCLE and GDSC cell line information
ccle_info <- get_cell_line_info(ccle_expression, "cell_line")
gdsc_info <- get_cell_line_info(gdsc_response, "CELL_LINE_NAME")

# Filter for A375 cell line in CCLE and GDSC data
ccle_a375 <- ccle_expression %>%
  filter(str_detect(cell_line, "A375"))

gdsc_a375 <- gdsc_response %>%
  filter(str_detect(CELL_LINE_NAME, "A375"))
```

This section defines a function to extract unique cell line information and filters the data for the A375 melanoma cell line.

## Gene Expression Analysis

```R
# Analyze gene expression data for A375
gene_expression_a375 <- ccle_a375 %>%
  arrange(desc(rna_expression)) %>%
  select(gene_name, rna_expression)

print("Top 20 highly expressed genes in A375:")
print(head(gene_expression_a375, 20))

# Basic statistical summary of gene expression
gene_expression_summary <- summary(ccle_a375$rna_expression)
print("Summary of A375 gene expression:")
print(gene_expression_summary)

# Visualize gene expression distribution
pdf("a375_gene_expression_histogram.pdf")
hist(ccle_a375$rna_expression, main="Distribution of Gene Expression in A375", xlab="RNA Expression")
dev.off()
```

This section analyzes and visualizes gene expression data for the A375 cell line.

## Drug Response Analysis

```R
# Analyze drug response data for A375
drug_response_a375 <- gdsc_a375 %>%
  arrange(LN_IC50) %>%
  select(DRUG_NAME, PUTATIVE_TARGET, LN_IC50, AUC)

print("Top 20 most effective drugs on A375 (based on lowest LN_IC50):")
print(head(drug_response_a375, 20))

# Basic statistical summary of drug response
drug_response_summary <- summary(gdsc_a375$LN_IC50)
print("Summary of A375 drug response (LN_IC50):")
print(drug_response_summary)

# Visualize drug response distribution
pdf('a375_drug_response_histogram.pdf')
hist(gdsc_a375$LN_IC50, main = "Distribution of Drug Response in A375", xlab="LN_IC50")
dev.off()
```

This section analyzes and visualizes drug response data for the A375 cell line.

## Pathway Enrichment Analysis

```R
library(clusterProfiler)
library(org.Hs.eg.db)

# Select top 100 expressed genes
top_genes <- head(gene_expression$gene_name, 100)

# Perform GO enrichment analysis
go_enrichment <- enrichGO(gene = top_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# Print top enriched pathways
print("Top enriched biological processes:")
print(head(go_enrichment@result, 10))

# Plot top 20 enriched terms
pdf("a375_go_enrichment_plot.pdf", width = 12, height = 8)
barplot(go_enrichment, showCategory = 20)
dev.off()
```

This section performs Gene Ontology (GO) enrichment analysis on the top 100 expressed genes in the A375 cell line.

## Comparative Analysis

```R
# Define melanoma cell lines of interest
melanoma_lines <- c("A375", "M00921", "HMY1", "SKMEL24")

# Filter CCLE and GDSC data for selected melanoma lines
melanoma_expression <- ccle_expression %>%
  filter(str_detect(cell_line, paste(melanoma_lines, collapse = "|")))

melanoma_response <- gdsc_response %>%
  filter(str_detect(CELL_LINE_NAME, paste(melanoma_lines, collapse = "|")))

# Analyze top expressed genes across melanoma lines
top_genes_per_line <- melanoma_expression %>%
  group_by(cell_line) %>%
  top_n(100, rna_expression) %>%
  ungroup()

# Visualize top expressed genes
ggplot(top_genes_per_line, aes(x = cell_line, y = rna_expression, fill = cell_line)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Top 100 Expressed Genes Across Melanoma Cell Lines",
       x = "Cell Line", y = "RNA Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("top_genes_melanoma_lines.pdf", width = 10, height = 6)

# Analyze drug response across melanoma lines
drug_response_summary <- melanoma_response %>%
  group_by(CELL_LINE_NAME, DRUG_NAME) %>%
  summarise(mean_LN_IC50 = mean(LN_IC50, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(CELL_LINE_NAME) %>%
  top_n(-20, mean_LN_IC50)

# Visualize drug response
ggplot(drug_response_summary, aes(x = CELL_LINE_NAME, y = mean_LN_IC50, fill = DRUG_NAME)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Top 20 Most Effective Drugs Across Melanoma Cell Lines",
       x = "Cell Line", y = "Mean LN_IC50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("drug_response_melanoma_lines.pdf", width = 12, height = 8)
```

This section performs a comparative analysis of gene expression and drug response across multiple melanoma cell lines.

For more detailed explanations of specific functions or sections, please refer to the comments in the code or create an issue in this repository.
