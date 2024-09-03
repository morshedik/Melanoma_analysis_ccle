packages <- c('tidyverse', 'data.table', "devtools", "usethis", "rvest", "stringdist", 'readxl')

install.packages(packages) 
### for conflicted codes, ask to choose which?
##by using dplyr::filter() or plyr::filter() 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("depmap")

for (pkg in packages) {
  library(pkg, character.only = TRUE)
}
library(depmap)


# Load CCLE expression data
ccle_expression <- depmap_TPM()

# Save the data locally
write_csv(ccle_expression, "CCLE_expression.csv")

# Load CCLE expression data
ccle_expression <- fread("CCLE_expression.csv")
print(colnames(ccle_expression))


# Get the latest GDSC data URL
gdsc_url <- 'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx'

# Download the latest GDSC data
if (!file.exists("GDSC2_fitted_dose_response_latest.csv")) {
  download.file(gdsc_url, destfile = "GDSC2_fitted_dose_response_latest.csv")
}

# Load GDSC drug response data
gdsc_response <- fread("GDSC2_fitted_dose_response_latest.csv")
