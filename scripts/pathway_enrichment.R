########Pathway enrichment analysis##########################


# Load the results
setwd("~/Desktop/R_output/Final_data")
gene_expression <- read_csv("a375_gene_expression.csv")
drug_response <- read_csv("a375_drug_response.csv")

BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("org.Hs.eg.db", force=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)

# Select top 100 expressed genes
top_genes <- head(gene_expression$gene_name, 100)
print(top_genes)


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



# Save full results
output_direct <- "~/Desktop/R_output/Final_data"
dir.create(output_direct, showWarnings = FALSE)
write_csv(go_enrichment@result, file.path(output_direct,"a375_go_enrichment.csv"))



# Plot top 20 enriched terms
pdf("~/Desktop/R_output/Final_data/a375_go_enrichment_plot.pdf", width = 12, height = 8)
barplot(go_enrichment, showCategory = 20)
dev.off()