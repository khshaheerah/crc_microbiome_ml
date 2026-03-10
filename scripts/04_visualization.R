# scripts/04_visualization_FINAL.R
# Visualization of results from ML model

# Load required libraries
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(gridExtra)
library(corrplot)

# Set your path (CHANGE THIS TO YOUR PATH)
setwd("/Users/shaheerah/R_analysis/crc_microbiome_ml")

# Create output directory if it doesn't exist
if(!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# ============================================
# 1. Load data and results
# ============================================
cat("\n📊 STEP 1: Loading data and results...\n")

# Load original data
french_data <- readRDS("data/french_data.rds")
feat_fr <- french_data$features
meta_fr <- french_data$metadata

# Load model results
if(file.exists("output/feature_weights.csv")) {
  feature_weights <- read.csv("output/feature_weights.csv")
  cat("✅ Loaded feature weights with", nrow(feature_weights), "features\n")
} else {
  stop("❌ feature_weights.csv not found! Run Script 3 first.")
}

if(file.exists("output/performance_metrics.csv")) {
  performance <- read.csv("output/performance_metrics.csv")
  cat("✅ Loaded performance metrics\n")
}

# Align samples
common_samples <- intersect(colnames(feat_fr), rownames(meta_fr))
feat_fr <- feat_fr[, common_samples]
meta_fr <- meta_fr[common_samples, , drop = FALSE]

# Convert to relative abundance
feat_fr_rel <- prop.table(as.matrix(feat_fr), 2) * 100

# ============================================
# 2. Heatmap of top features
# ============================================
cat("\n🔥 STEP 2: Creating heatmap of top features...\n")

# Get top 25 features (mix of positive and negative)
top25 <- feature_weights %>%
  arrange(desc(abs(coefficient))) %>%
  head(25)

top25_features <- top25$feature

# Check which features are in our data
available_features <- intersect(top25_features, rownames(feat_fr_rel))
cat("Features available for heatmap:", length(available_features), "/ 25\n")

if(length(available_features) < 5) {
  cat("⚠️  Not enough top features. Using all available features from model...\n")
  available_features <- intersect(feature_weights$feature[1:50], rownames(feat_fr_rel))
  available_features <- head(available_features, 25)
}

# Extract heatmap data
heatmap_data <- feat_fr_rel[available_features, ]

# Log transform for better visualization
heatmap_data_log <- log1p(heatmap_data)

# Z-score normalize for heatmap (scale across samples)
heatmap_data_scaled <- t(scale(t(heatmap_data_log)))

# Create annotation for samples
annotation_col <- data.frame(
  Group = meta_fr$Group,
  row.names = colnames(heatmap_data)
)

# Create annotation colors
ann_colors <- list(
  Group = c(CTR = "#377EB8", CRC = "#E41A1C")
)

# Create heatmap
png("output/figures/heatmap_top_features.png", width = 1200, height = 800, res = 100)
pheatmap(
  heatmap_data_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  main = "Top CRC-Associated Features - Abundance Patterns (Z-score)",
  fontsize_row = 10,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100)
)
dev.off()
cat("✅ Saved: output/figures/heatmap_top_features.png\n")

# ============================================
# 3. Enhanced feature weight plot
# ============================================
cat("\n📊 STEP 3: Creating enhanced feature weight plot...\n")

# Separate positive and negative
positive_features <- feature_weights %>%
  filter(coefficient > 0) %>%
  arrange(desc(coefficient)) %>%
  mutate(type = "Enriched in CRC")

negative_features <- feature_weights %>%
  filter(coefficient < 0) %>%
  arrange(coefficient) %>%
  mutate(type = "Enriched in Controls")

# Take top 15 from each
top_positive <- head(positive_features, 15)
top_negative <- head(negative_features, 15)
plot_features <- bind_rows(top_positive, top_negative)

# Create plot
p_weights <- ggplot(plot_features, 
                    aes(x = reorder(feature, coefficient), 
                        y = coefficient, 
                        fill = type)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Top Features Associated with Colorectal Cancer",
       subtitle = paste("Based on Lasso model (", nrow(feature_weights), 
                        " total features with non-zero weights)"),
       x = "Bacterial Features",
       y = "Model Coefficient (weight)") +
  scale_fill_manual(values = c("Enriched in CRC" = "#E41A1C", 
                               "Enriched in Controls" = "#377EB8"),
                    name = "Association") +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10, face = "italic"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray50")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

ggsave("output/figures/feature_weights_enhanced.png", p_weights, 
       width = 12, height = 8, dpi = 300)
cat("✅ Saved: output/figures/feature_weights_enhanced.png\n")

# ============================================
# 4. Performance comparison plot
# ============================================
cat("\n📈 STEP 4: Creating performance comparison plot...\n")

if(exists("performance")) {
  # Reshape data for plotting
  perf_long <- performance %>%
    pivot_longer(cols = c(Training, Test), 
                 names_to = "Dataset", 
                 values_to = "Value")
  
  p_perf <- ggplot(perf_long, aes(x = Metric, y = Value, fill = Dataset)) +
    geom_col(position = "dodge", width = 0.7) +
    theme_bw(base_size = 14) +
    labs(title = "Model Performance Comparison",
         subtitle = "Training (French) vs Test (Chinese) cohorts",
         y = "Score",
         x = "") +
    scale_fill_manual(values = c("Training" = "#377EB8", "Test" = "#E41A1C")) +
    geom_text(aes(label = round(Value, 3)), 
              position = position_dodge(0.7), 
              vjust = -0.5,
              size = 4) +
    ylim(0, 1.1) +
    theme(legend.position = "bottom")
  
  ggsave("output/figures/performance_comparison.png", p_perf, 
         width = 10, height = 6, dpi = 300)
  cat("✅ Saved: output/figures/performance_comparison.png\n")
}

# ============================================
# 5. Abundance boxplots for top features
# ============================================
cat("\n📦 STEP 5: Creating abundance boxplots...\n")

# Select top 6 features from each group
top6_positive <- head(positive_features$feature, 6)
top6_negative <- head(negative_features$feature, 6)
top12_features <- c(top6_positive, top6_negative)

# Prepare data for boxplots
boxplot_data <- data.frame()
for(feature in top12_features) {
  if(feature %in% rownames(feat_fr_rel)) {
    temp_df <- data.frame(
      Feature = feature,
      Abundance = as.numeric(feat_fr_rel[feature, ]),
      Group = meta_fr$Group,
      stringsAsFactors = FALSE
    )
    boxplot_data <- rbind(boxplot_data, temp_df)
  }
}

# Shorten feature names for plotting
boxplot_data$Feature_short <- gsub(".*\\[|\\].*", "", boxplot_data$Feature)
boxplot_data$Feature_short <- substr(boxplot_data$Feature_short, 1, 20)

if(nrow(boxplot_data) > 0) {
  # Create boxplots
  p_box <- ggplot(boxplot_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    facet_wrap(~Feature_short, scales = "free_y", ncol = 3) +
    theme_bw(base_size = 10) +
    labs(title = "Abundance Patterns of Top Discriminative Features",
         subtitle = "CRC vs Control groups",
         y = "Relative Abundance (%)") +
    scale_fill_manual(values = c("CTR" = "#377EB8", "CRC" = "#E41A1C")) +
    theme(
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif", 
                       hide.ns = TRUE,
                       size = 3)
  
  ggsave("output/figures/top_features_boxplots.png", p_box, 
         width = 14, height = 10, dpi = 300)
  cat("✅ Saved: output/figures/top_features_boxplots.png\n")
}

# ============================================
# 6. Correlation heatmap of top features
# ============================================
cat("\n🔗 STEP 6: Creating correlation heatmap...\n")

if(length(available_features) >= 5) {
  # Calculate correlation between top features
  top_feat_data <- t(feat_fr_rel[available_features, ])
  cor_matrix <- cor(top_feat_data, method = "spearman")
  
  # Create correlation heatmap
  png("output/figures/feature_correlation.png", width = 1000, height = 1000, res = 100)
  corrplot(cor_matrix, 
           method = "color",
           type = "upper",
           order = "hclust",
           tl.cex = 0.8,
           tl.col = "black",
           title = "Correlation Between Top Features",
           mar = c(0,0,2,0))
  dev.off()
  cat("✅ Saved: output/figures/feature_correlation.png\n")
}

# ============================================
# 7. Create summary report
# ============================================
cat("\n📝 STEP 7: Creating summary report...\n")

sink("output/analysis_summary.txt")
cat("========================================\n")
cat("MICROBIOME CRC ANALYSIS SUMMARY\n")
cat("========================================\n\n")
cat("Date:", date(), "\n\n")
cat("DATASETS\n")
cat("--------\n")
cat("Training cohort: French (Zeller et al.)\n")
cat("  - Samples:", ncol(feat_fr), "\n")
cat("  - Controls:", sum(meta_fr$Group == "CTR"), "\n")
cat("  - CRC cases:", sum(meta_fr$Group == "CRC"), "\n")
cat("  - Features:", nrow(feat_fr), "\n\n")
cat("Test cohort: Chinese (Yu et al.)\n\n")

cat("MODEL PERFORMANCE\n")
cat("----------------\n")
if(exists("performance")) {
  print(performance)
}

cat("\n\nTOP 10 CRC-ENRICHED FEATURES\n")
cat("---------------------------\n")
print(head(positive_features, 10))

cat("\n\nTOP 10 CONTROL-ENRICHED FEATURES\n")
cat("-------------------------------\n")
print(head(negative_features, 10))

cat("\n\nINTERPRETATION\n")
cat("--------------\n")
cat("• The Lasso model selected", nrow(feature_weights), "features with non-zero weights\n")
cat("• Training AUC:", round(performance$Training[performance$Metric == "AUC"], 3), "\n")
cat("• Test AUC:", round(performance$Test[performance$Metric == "AUC"], 3), "\n")
cat("• The drop in performance from training to test is expected due to:\n")
cat("  - Different populations (French vs Chinese)\n")
cat("  - Different sequencing platforms\n")
cat("  - Technical variation between studies\n")
cat("\n• Known CRC-associated genera identified:\n")
if(any(grepl("Fusobacterium", feature_weights$feature))) {
  cat("  - Fusobacterium nucleatum (well-established CRC marker)\n")
}
if(any(grepl("Peptostreptococcus", feature_weights$feature))) {
  cat("  - Peptostreptococcus spp.\n")
}
if(any(grepl("Porphyromonas", feature_weights$feature))) {
  cat("  - Porphyromonas spp.\n")
}

sink()
cat("✅ Saved: output/analysis_summary.txt\n")

# ============================================
# 8. Create combined figure for publication
# ============================================
cat("\n🎯 STEP 8: Creating combined figure...\n")

if(file.exists("output/figures/roc_curves.png") & 
   file.exists("output/figures/feature_weights_enhanced.png")) {
  
  # Load the PNG files
  library(png)
  library(grid)
  
  img_roc <- readPNG("output/figures/roc_curves.png")
  img_weights <- readPNG("output/figures/feature_weights_enhanced.png")
  
  # Create combined figure
  png("output/figures/combined_figure.png", width = 2400, height = 1200, res = 200)
  grid.newpage()
  
  # Layout: 1 row, 2 columns
  pushViewport(viewport(layout = grid.layout(1, 2)))
  
  # ROC plot
  pushViewport(viewport(layout.pos.col = 1))
  grid.raster(img_roc)
  popViewport()
  
  # Feature weights plot
  pushViewport(viewport(layout.pos.col = 2))
  grid.raster(img_weights)
  popViewport()
  
  dev.off()
  cat("✅ Saved: output/figures/combined_figure.png\n")
}

# ============================================
# 9. Create README for output folder
# ============================================
cat("\n📖 STEP 9: Creating README for output folder...\n")

readme_text <- c(
  "# Output Files Description",
  "",
  "## Figures",
  "- `roc_curves.png`: ROC curves for training and test sets",
  "- `cv_curve.png`: Lasso cross-validation curve",
  "- `feature_weights.png`: Barplot of feature weights",
  "- `feature_weights_enhanced.png`: Enhanced feature weight plot with colors",
  "- `heatmap_top_features.png`: Heatmap of top discriminative features",
  "- `top_features_boxplots.png`: Boxplots of abundance for top features",
  "- `feature_correlation.png`: Correlation between top features",
  "- `performance_comparison.png`: Barplot comparing training/test performance",
  "- `combined_figure.png`: Combined ROC and feature weights for publication",
  "",
  "## Data Files",
  "- `feature_weights.csv`: All features with non-zero coefficients",
  "- `performance_metrics.csv`: Model performance metrics",
  "- `analysis_summary.txt`: Text summary of all results",
  "",
  "## Models",
  "- `lasso_model_complete.rds`: Complete Lasso model object",
  "",
  "Generated on: ", date()
)

writeLines(readme_text, "output/README_output.txt")
cat("✅ Saved: output/README_output.txt\n")