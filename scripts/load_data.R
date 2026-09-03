packages <- c("tidyverse", "data.table", "readxl", "depmap")
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install dependencies before running the analysis; see requirements.R."
  )
}

for (pkg in packages) {
  library(pkg, character.only = TRUE)
}
library(depmap)

data_dir <- Sys.getenv("MELANOMA_DATA_DIR", "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)


# Load CCLE expression data
ccle_expression <- depmap_TPM()

# Save the data locally
ccle_path <- file.path(data_dir, "CCLE_expression.csv")
write_csv(ccle_expression, ccle_path)
print(colnames(ccle_expression))


# Get the latest GDSC data URL
gdsc_url <- 'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx'

# Download the latest GDSC data
gdsc_path <- file.path(data_dir, "GDSC2_fitted_dose_response_27Oct23.xlsx")
if (!file.exists(gdsc_path)) {
  download.file(gdsc_url, destfile = gdsc_path, mode = "wb")
}

# Load the Excel workbook with the matching reader. The earlier script saved this
# workbook with a .csv suffix and passed it to fread(), which is not a valid import.
gdsc_response <- readxl::read_excel(gdsc_path)

required_gdsc_columns <- c(
  "CELL_LINE_NAME", "DRUG_NAME", "PUTATIVE_TARGET", "LN_IC50", "AUC"
)
missing_gdsc_columns <- setdiff(required_gdsc_columns, names(gdsc_response))
if (length(missing_gdsc_columns) > 0) {
  stop(
    "Unexpected GDSC schema; missing columns: ",
    paste(missing_gdsc_columns, collapse = ", ")
  )
}
