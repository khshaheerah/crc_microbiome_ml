# scripts/01_load_data_FINAL.R
# Load libraries
library(tidyverse)

# Set your path
setwd("/Users/shaheerah/R_analysis/crc_microbiome_ml")

# Create data directory if it doesn't exist
if(!dir.exists("data")) dir.create("data")

# ============================================
# Download files (if not already downloaded)
# ============================================
cat("Checking for existing files...\n")

if(!file.exists("data/specI_Zeller.tsv")) {
  cat("Downloading Zeller et al. (France) dataset...\n")
  download.file(
    url = "https://zenodo.org/record/815875/files/specI_Zeller.tsv?download=1",
    destfile = "data/specI_Zeller.tsv",
    mode = "wb"
  )
}

if(!file.exists("data/meta_Zeller.tsv")) {
  download.file(
    url = "https://zenodo.org/record/815875/files/meta_Zeller.tsv?download=1",
    destfile = "data/meta_Zeller.tsv",
    mode = "wb"
  )
}

if(!file.exists("data/specI_Yu.tsv")) {
  cat("Downloading Yu et al. (China) dataset...\n")
  download.file(
    url = "https://zenodo.org/record/815875/files/specI_Yu.tsv?download=1",
    destfile = "data/specI_Yu.tsv",
    mode = "wb"
  )
}

if(!file.exists("data/meta_Yu.tsv")) {
  download.file(
    url = "https://zenodo.org/record/815875/files/meta_Yu.tsv?download=1",
    destfile = "data/meta_Yu.tsv",
    mode = "wb"
  )
}

# ============================================
# Read the files correctly
# ============================================
cat("\n📖 Reading files...\n")

# Read feature tables (species abundances)
# Format: rows = species, columns = samples
# The first column contains species names (row names)
feat_fr <- read.table("data/specI_Zeller.tsv", 
                      sep = "\t", 
                      row.names = 1,           # First column is row names
                      header = TRUE,            # First row is column names
                      check.names = FALSE,      # Keep original sample names
                      fill = TRUE,               # Fill missing values if any
                      quote = "")                # No quotes

feat_cn <- read.table("data/specI_Yu.tsv", 
                      sep = "\t", 
                      row.names = 1, 
                      header = TRUE,
                      check.names = FALSE,
                      fill = TRUE,
                      quote = "")

# Read metadata
meta_fr <- read.table("data/meta_Zeller.tsv", 
                      sep = "\t", 
                      row.names = 1, 
                      header = TRUE,
                      check.names = FALSE,
                      quote = "")

meta_cn <- read.table("data/meta_Yu.tsv", 
                      sep = "\t", 
                      row.names = 1, 
                      header = TRUE,
                      check.names = FALSE,
                      quote = "")

# ============================================
# Clean up the data
# ============================================
cat("\n🧹 Cleaning data...\n")

# Remove any empty or NA columns
feat_fr <- feat_fr[, colSums(is.na(feat_fr)) == 0]
feat_cn <- feat_cn[, colSums(is.na(feat_cn)) == 0]

# Remove any rows with all zeros (optional)
feat_fr <- feat_fr[rowSums(feat_fr) > 0, ]
feat_cn <- feat_cn[rowSums(feat_cn) > 0, ]

# ============================================
# Check the data structure
# ============================================
cat("\n=== Data Summary ===\n")
cat("\n🇫🇷 French dataset:\n")
cat("  Features (species):", nrow(feat_fr), "\n")
cat("  Samples:", ncol(feat_fr), "\n")
cat("  Metadata samples:", nrow(meta_fr), "\n")
cat("  Class distribution:\n")
print(table(meta_fr$Group))

cat("\n🇨🇳 Chinese dataset:\n")
cat("  Features (species):", nrow(feat_cn), "\n")
cat("  Samples:", ncol(feat_cn), "\n")
cat("  Metadata samples:", nrow(meta_cn), "\n")
cat("  Class distribution:\n")
print(table(meta_cn$Group))

# Check if sample names match between features and metadata
fr_common <- intersect(colnames(feat_fr), rownames(meta_fr))
cn_common <- intersect(colnames(feat_cn), rownames(meta_cn))

cat("\n✅ Sample matching:\n")
cat("  French - matching samples:", length(fr_common), "/", ncol(feat_fr), "\n")
cat("  Chinese - matching samples:", length(cn_common), "/", ncol(feat_cn), "\n")

# ============================================
# Save as RDS for easier loading later
# ============================================
cat("\n💾 Saving processed data...\n")

saveRDS(list(features = feat_fr, 
             metadata = meta_fr, 
             common_samples = fr_common), 
        "data/french_data.rds")

saveRDS(list(features = feat_cn, 
             metadata = meta_cn,
             common_samples = cn_common), 
        "data/chinese_data.rds")

cat("\n✅ Data downloaded, cleaned, and saved!\n")
cat("   Files saved as RDS for faster loading:\n")
cat("   - data/french_data.rds\n")
cat("   - data/chinese_data.rds\n")

# ============================================
# Quick preview
# ============================================
cat("\n👀 Quick preview of the data:\n")
cat("\nFirst 5 sample names (French):\n")
print(head(colnames(feat_fr), 5))

cat("\nFirst 5 species (French):\n")
print(head(rownames(feat_fr), 5))

cat("\nFirst 5 metadata rows (French):\n")
print(head(meta_fr, 5))
