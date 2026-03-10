# install_packages.R - Run this first, working directory doesn't matter
install.packages(c("tidyverse", "remotes", "vegan", "pheatmap", "ggpubr"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SIAMCAT", "phyloseq"))

# Verify installations
library(SIAMCAT)
library(phyloseq)
print("Packages installed successfully!")
