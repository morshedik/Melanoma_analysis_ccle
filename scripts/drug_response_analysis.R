# Filter for A375 in GDSC
gdsc_a375 <- gdsc_response %>%
  filter(str_detect(CELL_LINE_NAME, "A375"))
print(gdsc_a375)


# Analyze drug response data for A375
drug_response_a375 <- gdsc_a375 %>%
  arrange(LN_IC50) %>%
  select(DRUG_NAME, PUTATIVE_TARGET, LN_IC50, AUC)
print("Top 20 most effective drugs on A375 (based on lowest LN_IC50):")
print(head(drug_response_a375, 20))



# Save results
output_direct <- "~/Desktop/R_output/Final_data"
dir.create(output_direct, showWarnings = FALSE)
write_csv(drug_response_a375, file.path(output_direct,"a375_drug_response.csv"))


# Basic statistical summary of drug response
drug_response_summary <- summary(gdsc_a375$LN_IC50)
print("Summary of A375 drug response (LN_IC50):")
print(drug_response_summary)


pdf('~/Desktop/R_output/Final_data/a375_drug_response_histogram.pdf')
hist(gdsc_a375$LN_IC50, main = "Distribution of Drug Response in A375", xlab="LN_IC50")
dev.off()
