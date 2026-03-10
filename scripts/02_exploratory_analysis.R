# scripts/02_exploratory_analysis_FINAL.R
# Exploratory analysis of microbiome data

# Load required libraries
library(tidyverse)
library(vegan)
library(ggpubr)

# Set your path (CHANGE THIS TO YOUR PATH)
setwd("/Users/shaheerah/R_analysis/crc_microbiome_ml")

# Create output directory if it doesn't exist
if(!dir.exists("output/figures")) {
  dir.create("output/figures", recursive = TRUE)
}

# ============================================
# 1. Load data from RDS files
# ============================================
cat("\n📊 STEP 1: Loading data...\n")

french_data <- readRDS("data/french_data.rds")
feat_fr <- french_data$features
meta_fr <- french_data$metadata

cat("✅ Data loaded successfully\n")
cat("   Features:", nrow(feat_fr), "\n")
cat("   Samples:", ncol(feat_fr), "\n")

# ============================================
# 2. Align samples between features and metadata
# ============================================
cat("\n🔧 STEP 2: Aligning samples...\n")

common_samples <- intersect(colnames(feat_fr), rownames(meta_fr))
feat_fr <- feat_fr[, common_samples]
meta_fr <- meta_fr[common_samples, , drop = FALSE]

cat("✅ Aligned", length(common_samples), "samples\n")
cat("   Class distribution:\n")
print(table(meta_fr$Group))

# ============================================
# 3. Convert to relative abundances
# ============================================
cat("\n📈 STEP 3: Calculating relative abundances...\n")

feat_fr_mat <- as.matrix(feat_fr)
feat_fr_rel <- prop.table(feat_fr_mat, 2) * 100  # As percentages

# ============================================
# 4. Alpha Diversity Analysis
# ============================================
cat("\n📊 STEP 4: Calculating alpha diversity...\n")

# Calculate Shannon diversity
alpha_div <- diversity(t(feat_fr_mat), index = "shannon")
alpha_df <- data.frame(
  Sample = names(alpha_div),
  Shannon = alpha_div,
  Group = meta_fr[names(alpha_div), "Group"]
)

# Calculate additional diversity metrics
alpha_df$Richness <- specnumber(t(feat_fr_mat))  # Species richness
alpha_df$Simpson <- diversity(t(feat_fr_mat), index = "simpson")  # Simpson index

# Save diversity metrics
write.csv(alpha_df, "output/alpha_diversity_metrics.csv", row.names = FALSE)

# Plot Shannon diversity
p_alpha_shannon <- ggplot(alpha_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  theme_bw(base_size = 14) +
  labs(title = "Alpha Diversity (Shannon Index)",
       subtitle = paste("French Cohort (n =", nrow(alpha_df), "samples)"),
       y = "Shannon Diversity", x = "") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("CTR" = "#377EB8", "CRC" = "#E41A1C")) +
  theme(legend.position = "none")

ggsave("output/figures/alpha_diversity_shannon.png", 
       p_alpha_shannon, width = 6, height = 5, dpi = 300)

# Plot species richness
p_alpha_richness <- ggplot(alpha_df, aes(x = Group, y = Richness, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  theme_bw(base_size = 14) +
  labs(title = "Species Richness",
       y = "Number of Species", x = "") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("CTR" = "#377EB8", "CRC" = "#E41A1C")) +
  theme(legend.position = "none")

ggsave("output/figures/alpha_diversity_richness.png", 
       p_alpha_richness, width = 6, height = 5, dpi = 300)

# ============================================
# 5. Top Abundant Species
# ============================================
cat("\n📊 STEP 5: Finding top abundant species...\n")

# Calculate mean abundance across all samples
mean_abundance <- rowMeans(feat_fr_rel)
top15 <- names(sort(mean_abundance, decreasing = TRUE))[1:15]

cat("Top 15 most abundant species:\n")
print(data.frame(
  Species = top15,
  `Mean Abundance (%)` = round(sort(mean_abundance, decreasing = TRUE)[1:15], 3)
))

# Prepare data for boxplot
top15_data <- feat_fr_rel[top15, ] %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "Abundance") %>%
  left_join(meta_fr %>% rownames_to_column("Sample"), by = "Sample")

# Plot top 15 species
p_top <- ggplot(top15_data, aes(x = reorder(Species, -Abundance, median), 
                                y = Abundance, 
                                fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 15 Most Abundant Species",
       subtitle = "By group (CRC vs Control)",
       y = "Relative Abundance (%)", 
       x = "") +
  scale_fill_manual(values = c("CTR" = "#377EB8", "CRC" = "#E41A1C")) +
  theme(legend.position = "bottom")

ggsave("output/figures/top_species_boxplot.png", 
       p_top, width = 12, height = 7, dpi = 300)

# ============================================
# 6. Beta Diversity (PCoA)
# ============================================
cat("\n📊 STEP 6: Calculating beta diversity...\n")

# Calculate Bray-Curtis distance
dist_mat <- vegdist(t(feat_fr_rel), method = "bray")

# Perform PCoA
pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  PC1 = pcoa$points[,1],
  PC2 = pcoa$points[,2],
  Group = meta_fr$Group
)

# Calculate variance explained
var_exp <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

# Perform PERMANOVA to test for group differences
adonis_result <- adonis2(dist_mat ~ Group, data = meta_fr, permutations = 999)
cat("\nPERMANOVA results:\n")
print(adonis_result)

# Save PERMANOVA results
sink("output/permanova_results.txt")
cat("PERMANOVA (Bray-Curtis) - Group Differences\n")
cat("============================================\n\n")
print(adonis_result)
sink()

# Plot PCoA
p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  theme_bw(base_size = 14) +
  labs(title = "Beta Diversity (Bray-Curtis PCoA)",
       subtitle = paste("PERMANOVA p =", format(adonis_result$`Pr(>F)`[1], digits = 3)),
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  scale_color_manual(values = c("CTR" = "#377EB8", "CRC" = "#E41A1C")) +
  theme(legend.position = "bottom")

ggsave("output/figures/beta_diversity_pcoa.png", 
       p_pcoa, width = 8, height = 6, dpi = 300)

# ============================================
# 7. Differential Abundance (Wilcoxon test)
# ============================================
cat("\n📊 STEP 7: Testing differential abundance...\n")

# Perform Wilcoxon test for each species
species_names <- rownames(feat_fr_rel)
p_values <- numeric(length(species_names))
fold_changes <- numeric(length(species_names))

for(i in 1:length(species_names)) {
  species <- species_names[i]
  ctr_abund <- feat_fr_rel[species, meta_fr$Group == "CTR"]
  crc_abund <- feat_fr_rel[species, meta_fr$Group == "CRC"]
  
  # Wilcoxon test
  test_result <- wilcox.test(ctr_abund, crc_abund)
  p_values[i] <- test_result$p.value
  
  # Fold change (log2 ratio)
  fold_changes[i] <- log2(mean(crc_abund) + 0.01) - log2(mean(ctr_abund) + 0.01)
}

# Adjust p-values for multiple testing
p_adjusted <- p.adjust(p_values, method = "fdr")

# Create results dataframe
diff_abund <- data.frame(
  Species = species_names,
  p_value = p_values,
  p_adjusted = p_adjusted,
  log2FC = fold_changes,
  mean_abundance = mean_abundance
) %>%
  arrange(p_adjusted)

# Save results
write.csv(diff_abund, "output/differential_abundance.csv", row.names = FALSE)

# Show significant species (FDR < 0.05)
sig_species <- diff_abund %>% filter(p_adjusted < 0.05)
cat("\nSignificant species (FDR < 0.05):", nrow(sig_species), "\n")
if(nrow(sig_species) > 0) {
  print(head(sig_species, 10))
}

# Plot volcano plot
p_volcano <- ggplot(diff_abund, aes(x = log2FC, y = -log10(p_adjusted), 
                                    color = p_adjusted < 0.05)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_bw(base_size = 12) +
  labs(title = "Volcano Plot: Differential Abundance",
       x = "Log2 Fold Change (CRC/CTR)",
       y = "-Log10(Adjusted p-value)") +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "gray50"),
                     name = "Significant (FDR < 0.05)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme(legend.position = "bottom")

ggsave("output/figures/volcano_plot.png", 
       p_volcano, width = 8, height = 6, dpi = 300)

# ============================================
# 8. Summary Statistics
# ============================================
cat("\n📝 STEP 8: Generating summary statistics...\n")

summary_stats <- data.frame(
  Metric = c("Total Samples", "Controls (CTR)", "Cases (CRC)",
             "Total Species", "Mean Species per Sample",
             "Shannon Diversity (CTR)", "Shannon Diversity (CRC)",
             "Significant Species (FDR < 0.05)"),
  Value = c(
    nrow(meta_fr),
    sum(meta_fr$Group == "CTR"),
    sum(meta_fr$Group == "CRC"),
    nrow(feat_fr),
    round(mean(alpha_df$Richness), 1),
    round(mean(alpha_df$Shannon[alpha_df$Group == "CTR"]), 2),
    round(mean(alpha_df$Shannon[alpha_df$Group == "CRC"]), 2),
    nrow(sig_species)
  )
)

write.csv(summary_stats, "output/exploratory_summary.csv", row.names = FALSE)
