# scripts/03_ml_FINAL.R
# Machine Learning: Lasso model for CRC prediction

# Load required libraries
library(tidyverse)
library(glmnet)
library(pROC)

# Set your path (CHANGE THIS TO YOUR PATH)
setwd("/Users/<USERNAME>/Path to/crc_microbiome_ml")

# Create output directories if they don't exist
if(!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
if(!dir.exists("output/models")) dir.create("output/models", recursive = TRUE)

# ============================================
# 1. Load data from RDS files
# ============================================
cat("\n📊 STEP 1: Loading data...\n")

french_data <- readRDS("data/french_data.rds")
chinese_data <- readRDS("data/chinese_data.rds")

feat_fr <- french_data$features
meta_fr <- french_data$metadata
feat_cn <- chinese_data$features
meta_cn <- chinese_data$metadata

cat("French data:", ncol(feat_fr), "samples,", nrow(feat_fr), "features\n")
cat("Chinese data:", ncol(feat_cn), "samples,", nrow(feat_cn), "features\n")

# ============================================
# 2. Preprocess training data (French)
# ============================================
cat("\n🔧 STEP 2: Preprocessing training data...\n")

# Align French samples
fr_common <- intersect(colnames(feat_fr), rownames(meta_fr))
feat_fr <- feat_fr[, fr_common]
meta_fr <- meta_fr[fr_common, , drop = FALSE]

# Convert to relative abundance (percentage)
feat_fr_rel <- prop.table(as.matrix(feat_fr), 2) * 100

# Transpose to samples x features
X_train <- t(feat_fr_rel)

# Create outcome (CTR = 0, CRC = 1)
y_train <- ifelse(meta_fr$Group == "CRC", 1, 0)
names(y_train) <- rownames(meta_fr)

# Remove low abundance features (< 0.01% in at least 10% of samples)
feat_present <- colSums(X_train > 0.01) / nrow(X_train)
X_train <- X_train[, feat_present > 0.1, drop = FALSE]

cat("Training set:", nrow(X_train), "samples,", ncol(X_train), "features\n")
cat("Class distribution:", table(y_train), "\n")

# ============================================
# 3. Preprocess test data (Chinese)
# ============================================
cat("\n STEP 3: Preprocessing test data...\n")

# Align Chinese samples
cn_common <- intersect(colnames(feat_cn), rownames(meta_cn))
feat_cn <- feat_cn[, cn_common]
meta_cn <- meta_cn[cn_common, , drop = FALSE]

# Convert to relative abundance
feat_cn_rel <- prop.table(as.matrix(feat_cn), 2) * 100

# Transpose
X_test_raw <- t(feat_cn_rel)

# Create outcome
y_test <- ifelse(meta_cn$Group == "CRC", 1, 0)
names(y_test) <- rownames(meta_cn)

# Find common features with training set
common_features <- intersect(colnames(X_train), colnames(X_test_raw))
X_train <- X_train[, common_features]
X_test <- X_test_raw[, common_features]

cat("Test set:", nrow(X_test), "samples,", ncol(X_test), "features\n")
cat("Class distribution:", table(y_test), "\n")
cat("Common features:", length(common_features), "\n")

# ============================================
# 4. Train Lasso model with cross-validation
# ============================================
cat("\n STEP 4: Training Lasso model...\n")

set.seed(123)  # For reproducibility

# Cross-validated Lasso
cv_lasso <- cv.glmnet(
  x = X_train,
  y = y_train,
  family = "binomial",
  alpha = 1,  # 1 = Lasso
  nfolds = 5,
  type.measure = "auc",
  parallel = FALSE  # Set to TRUE if you want parallel processing
)

# Best lambda values
lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se

cat("Best lambda (min):", round(lambda_min, 4), "\n")
cat("Best lambda (1se):", round(lambda_1se, 4), "\n")

# Number of non-zero features at best lambda
n_features <- cv_lasso$nzero[which(cv_lasso$lambda == lambda_min)]
cat("Non-zero features:", n_features, "\n")

# ============================================
# 5. Make predictions
# ============================================
cat("\nSTEP 5: Making predictions...\n")

# Predictions on training set
pred_train <- predict(cv_lasso, newx = X_train, s = "lambda.min", type = "response")
pred_train_class <- ifelse(pred_train > 0.5, 1, 0)

# Predictions on test set
pred_test <- predict(cv_lasso, newx = X_test, s = "lambda.min", type = "response")
pred_test_class <- ifelse(pred_test > 0.5, 1, 0)

# ============================================
# 6. Calculate performance metrics
# ============================================
cat("\nSTEP 6: Calculating performance...\n")

# AUC
train_roc <- roc(y_train, as.vector(pred_train))
test_roc <- roc(y_test, as.vector(pred_test))

# Accuracy
train_acc <- mean(pred_train_class == y_train)
test_acc <- mean(pred_test_class == y_test)

# Sensitivity and Specificity
train_sens <- sum(pred_train_class == 1 & y_train == 1) / sum(y_train == 1)
train_spec <- sum(pred_train_class == 0 & y_train == 0) / sum(y_train == 0)
test_sens <- sum(pred_test_class == 1 & y_test == 1) / sum(y_test == 1)
test_spec <- sum(pred_test_class == 0 & y_test == 0) / sum(y_test == 0)

# Create performance dataframe
performance <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity"),
  Training = c(auc(train_roc), train_acc, train_sens, train_spec),
  Test = c(auc(test_roc), test_acc, test_sens, test_spec)
)

cat("\n=== Performance Metrics ===\n")
print(performance)

# Save performance
write.csv(performance, "output/performance_metrics.csv", row.names = FALSE)

# ============================================
# 7. Extract important features
# ============================================
cat("\nSTEP 7: Extracting important features...\n")

# Get coefficients at best lambda
coef_matrix <- as.matrix(coef(cv_lasso, s = "lambda.min"))
feature_weights <- data.frame(
  feature = rownames(coef_matrix)[-1],  # Remove intercept
  coefficient = coef_matrix[-1, 1]
) %>%
  filter(coefficient != 0) %>%
  arrange(desc(abs(coefficient)))

cat("Found", nrow(feature_weights), "features with non-zero coefficients\n")

# Show top features
cat("\nTop 10 features associated with CRC:\n")
top_positive <- feature_weights %>% filter(coefficient > 0) %>% head(10)
print(top_positive)

cat("\nTop 10 features associated with Controls:\n")
top_negative <- feature_weights %>% filter(coefficient < 0) %>% head(10)
print(top_negative)

# Save feature weights
write.csv(feature_weights, "output/feature_weights.csv", row.names = FALSE)

# ============================================
# 8. Create plots
# ============================================
cat("\nSTEP 8: Creating plots...\n")

# Plot 1: ROC Curves
png("output/figures/roc_curves.png", width = 800, height = 600, res = 100)
plot(train_roc, col = "#377EB8", lwd = 3, 
     main = "ROC Curves: Lasso Model Performance",
     cex.main = 1.5, cex.lab = 1.2)
plot(test_roc, col = "#E41A1C", lwd = 3, add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray50")
legend("bottomright", 
       legend = c(paste("Training (AUC =", round(auc(train_roc), 3), ")"),
                  paste("Test (AUC =", round(auc(test_roc), 3), ")")),
       col = c("#377EB8", "#E41A1C"), lwd = 3, cex = 1.2)
grid()
dev.off()
cat("Saved: output/figures/roc_curves.png\n")

# Plot 2: Cross-validation curve
png("output/figures/cv_curve.png", width = 800, height = 600, res = 100)
plot(cv_lasso, main = "Lasso Cross-Validation Curve",
     cex.main = 1.5, cex.lab = 1.2)
abline(v = log(lambda_min), col = "#E41A1C", lty = 2, lwd = 2)
abline(v = log(lambda_1se), col = "#377EB8", lty = 2, lwd = 2)
legend("topright", 
       legend = c("Lambda.min", "Lambda.1se"),
       col = c("#E41A1C", "#377EB8"), lty = 2, lwd = 2)
dev.off()
cat("Saved: output/figures/cv_curve.png\n")

# Plot 3: Feature weights (top 20)
top20 <- head(feature_weights, 20)
p_weights <- ggplot(top20, aes(x = reorder(feature, coefficient), 
                               y = coefficient, 
                               fill = coefficient > 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Top 20 Features Associated with CRC",
       subtitle = paste("Lasso model with", nrow(feature_weights), 
                        "non-zero features"),
       x = "Bacterial Features",
       y = "Model Coefficient") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    name = "Association",
                    labels = c("Enriched in CRC", "Enriched in Controls")) +
  theme(legend.position = "bottom")

ggsave("output/figures/feature_weights.png", p_weights, 
       width = 10, height = 7, dpi = 300)
cat("Saved: output/figures/feature_weights.png\n")

# ============================================
# 9. Save model and session info
# ============================================
cat("\nSTEP 9: Saving model...\n")

# Save model
saveRDS(list(
  model = cv_lasso,
  features = common_features,
  performance = performance,
  feature_weights = feature_weights
), "output/models/lasso_model_complete.rds")

# Save session info
sink("output/session_info.txt")
sessionInfo()
sink()
