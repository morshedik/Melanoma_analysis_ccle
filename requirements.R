# Run this file once when creating a controlled project environment.
# For publication work, replace unconstrained installation with renv::snapshot().

cran_packages <- c(
  "tidyverse",
  "data.table",
  "readxl",
  "rstatix",
  "pheatmap",
  "conflicted"
)

bioconductor_packages <- c(
  "depmap",
  "biomaRt",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

message("CRAN packages: ", paste(cran_packages, collapse = ", "))
message(
  "Bioconductor packages: ",
  paste(bioconductor_packages, collapse = ", ")
)
message("Install these in a project-local renv environment before running scripts.")
