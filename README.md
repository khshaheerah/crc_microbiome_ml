# Gut Microbiome Signatures in Colorectal Cancer: A Metagenomics & Machine Learning Project

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-≥4.0-blue)](https://www.r-project.org/)
[![GitHub last commit](https://img.shields.io/github/last-commit/khshaheerah/crc_microbiome_ml)](https://github.com/YOUR_USERNAME/crc_microbiome_ml)

## Overview

This project investigates the relationship between the gut microbiome and colorectal cancer (CRC) using public metagenomic sequencing data. We analyze microbial community differences between CRC patients and healthy controls, build a machine learning classifier to predict disease status from stool samples, and validate our findings on an independent cohort. The entire workflow runs on a standard Mac without requiring high-performance computing.

**Key questions addressed:**
- Does the gut microbiome differ between CRC patients and healthy individuals?
- Can we predict CRC status based on microbial composition?
- Which specific bacteria are most associated with CRC?
- Do these patterns generalize across different populations?

---

## Repository Structure

crc_microbiome_ml/
├── scripts/ # All R scripts (run in order)
│ ├── 01_load_data.R # Download and process datasets
│ ├── 02_exploratory_analysis.R # Alpha/beta diversity, differential abundance
│ ├── 03_ml_model.R # Lasso logistic regression model
│ └── 04_visualization.R # Publication-ready figures
├── README.md # This file
└── .gitignore # Files excluded from version control

**Note:** Data files and outputs are excluded from GitHub due to size limitations but will be generated when you run the scripts.

---

## Datasets

We use two publicly available colorectal cancer microbiome studies:

| Cohort | Study | Samples | Purpose |
|--------|-------|---------|---------|
| **French** | Zeller et al. (2014), *Mol Syst Biol* | 156 (CRC + CTR) | Training set |
| **Chinese** | Yu et al. (2017), *Gut* | 120 (CRC + CTR) | Independent test set |

**Data type:** Species-level taxonomic profiles (relative abundances) generated from shotgun metagenomic sequencing of stool samples.

*CTR = healthy control, CRC = colorectal cancer patient*

---

## Methods

### Data Processing
- Sample alignment between feature tables and metadata
- Relative abundance transformation (percentages)
- Filtering of low-abundance species (<0.01% in <10% of samples)

### Exploratory Analysis
- **Alpha diversity**: Shannon index, species richness (Wilcoxon tests)
- **Beta diversity**: Bray-Curtis dissimilarity, PCoA visualization, PERMANOVA
- **Differential abundance**: Wilcoxon tests with FDR correction, volcano plots

### Machine Learning
- **Algorithm**: Lasso logistic regression (L1 regularization)
- **Why Lasso?** Handles high-dimensional data (1000+ features), performs automatic feature selection, produces interpretable sparse models
- **Cross-validation**: 5-fold cross-validation on French cohort
- **Validation**: Independent testing on Chinese cohort
- **Evaluation metrics**: AUC-ROC, accuracy, sensitivity, specificity

### Visualization
- ROC curves for model performance
- Feature importance barplots
- Heatmaps of top discriminative species
- Correlation networks of key bacteria
- Combined publication-ready figures

---

## How to Reproduce

### Prerequisites
- **R ≥ 4.0** installed on your system
- Internet connection for data download

### Step-by-Step Instructions

1. **Clone this repository**
   ```bash
   git clone https://github.com/khshaheerah/crc_microbiome_ml.git
   cd crc_microbiome_ml

2. **Load and install packages**
   ```r
   packages <- c("tidyverse", "vegan", "ggpubr", "glmnet", "pROC", 
              "pheatmap", "gridExtra", "corrplot", "grid", "png")
   install.packages(packages[!packages %in% installed.packages()])

3. **Run the analysis pipeline**
   ```r
   setwd("/path/to/crc_microbiome_ml")  # Change to your path
   source("scripts/01_load_data.R")      # Downloads and processes data
   source("scripts/02_exploratory_analysis.R")  # Creates exploratory plots
   source("scripts/03_ml_model.R")       # Trains ML model
   source("scripts/04_visualization.R")  # Creates final visualizations


## Output Files

### Figures
- `roc_curves.png`: ROC curves for training and test sets
- `cv_curve.png`: Lasso cross-validation curve
- `feature_weights.png`: Barplot of feature weights
- `feature_weights_enhanced.png`: Enhanced feature weight plot with colors
- `heatmap_top_features.png`: Heatmap of top discriminative features
- `top_features_boxplots.png`: Boxplots of abundance for top features
- `feature_correlation.png`: Correlation between top features
- `performance_comparison.png`: Barplot comparing training/test performance
- `combined_figure.png`: Combined ROC and feature weights for publication

### Data Files
- `feature_weights.csv`: All features with non-zero coefficients
- `performance_metrics.csv`: Model performance metrics
- `analysis_summary.txt`: Text summary of all results

### Models
- `lasso_model_complete.rds`: Complete Lasso model object
