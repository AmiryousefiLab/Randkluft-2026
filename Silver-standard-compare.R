# Load required packages (install once)
#install.packages(c("mixtools", "ggplot2", "reshape2", "devtools"))
#devtools::install_github("JiangmeiRubyXiong/GammaGateR")  # For GammaGateR
#library(GammaGateR)
library(mixtools)
library(ggplot2)
library(reshape2)

# Randkluft: Find cutoff where skewness of truncated data <= c is ~0
randkluft_cutoff <- function(data, log_transform = TRUE) {
  if (log_transform) data <- log(data + 1)
  sorted <- sort(data)
  n <- length(sorted)
  best_sk <- Inf
  best_c <- sorted[floor(n * 0.95)]
  step <- max(10, floor(n / 1000))
  for (i in seq(500, n - 100, by = step)) {
    trunc <- sorted[1:i]
    m <- mean(trunc)
    third <- mean((trunc - m)^3)
    sd3 <- sd(trunc)^3
    if (sd3 == 0) next
    sk <- third / sd3
    if (abs(sk) < abs(best_sk)) {
      best_sk <- sk
      best_c <- sorted[i]
    }
  }
  best_c
}

gmm_cutoff <- function(data, log_transform = TRUE) {
  if (log_transform) data <- log(data + 1)
  
  # If data is too uniform or has too few unique values, return median
  if (length(unique(data)) < 10 || sd(data) < 1e-6) {
    return(median(data))
  }
  
  # Try EM with better control and multiple restarts
  best_cutoff <- median(data)  # fallback
  best_loglik <- -Inf
  
  for (attempt in 1:5) {  # Try up to 5 random initializations
    suppressWarnings({
      gmm <- try(mixtools::normalmixEM(
        data, 
        k = 2, 
        maxit = 1000,
        epsilon = 1e-4,   # looser tolerance
        verb = FALSE
      ), silent = TRUE)
    })
    
    if (!inherits(gmm, "try-error") && !is.null(gmm$loglik)) {
      if (gmm$loglik > best_loglik) {
        best_loglik <- gmm$loglik
        
        means <- gmm$mu
        sds <- gmm$sigma
        props <- gmm$lambda
        
        # Approximate valley between components
        f1 <- function(x) props[1] * dnorm(x, means[1], sds[1])
        f2 <- function(x) props[2] * dnorm(x, means[2], sds[2])
        opt <- optim(mean(means), function(c) abs(f1(c) - f2(c)), 
                     method = "Brent", lower = min(means)-2*max(sds), upper = max(means)+2*max(sds))
        best_cutoff <- opt$par
      }
    }
  }
  
  # Final fallback: use 95th percentile (common for rare positive tails)
  if (best_loglik == -Inf) {
    best_cutoff <- quantile(data, 0.95)
  }
  
  return(best_cutoff)
}
# GammaGateR cutoff (posterior probability > 0.5 for positive component)
gammagater_cutoff <- function(data, log_transform = TRUE) {
  if (log_transform) data <- log(data + 1)
  df <- data.frame(x = data)
  fit <- try(fit_gamma_mixture(df, channel = "x", lower_bound = min(data) + 0.1 * sd(data)), silent = TRUE)
  if (inherits(fit, "try-error")) return(median(data))
  probs <- fit$Posterior_Positive
  sort(data)[min(which(probs > 0.5))]
}

# Compute confusion matrix metrics
compute_metrics <- function(y_true, y_pred) {
  tn <- sum(y_true == 0 & y_pred == 0)
  fp <- sum(y_true == 0 & y_pred == 1)
  fn <- sum(y_true == 1 & y_pred == 0)
  tp <- sum(y_true == 1 & y_pred == 1)
  acc <- (tp + tn) / length(y_true)
  prec <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  rec <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1 <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
  spec <- ifelse(tn + fp > 0, tn / (tn + fp), 0)
  list(TP = tp, TN = tn, FP = fp, FN = fn, Accuracy = acc, Precision = prec, Recall = rec, F1 = f1, Specificity = spec)
}

# Cutoff difference metrics
cutoff_diffs <- function(gt_gate, rand_c, gmm_c, gamma_c) {
  list(
    Randkluft_diff = abs(rand_c - gt_gate),
    GMM_diff = abs(gmm_c - gt_gate),
    GammaGateR_diff = abs(gamma_c - gt_gate)
  )
}


######main processing functions. 



#Load your data (adjust paths if needed)
raw <- read.csv("/Users/ali/Desktop/Randkluft/raw.csv", stringsAsFactors = FALSE)
gate <- read.csv("/Users/ali/Desktop/Randkluft/tuulia_data_GT.csv", stringsAsFactors = FALSE)

# Initialize results
results <- data.frame(
  Patient = character(), Marker = character(),
  Randkluft_cutoff = numeric(), GMM_cutoff = numeric(), GammaGateR_cutoff = numeric(),
  Randkluft_diff = numeric(), GMM_diff = numeric(), GammaGateR_diff = numeric(),
  Randkluft_Acc = numeric(), GMM_Acc = numeric(), GammaGateR_Acc = numeric(),
  Randkluft_F1 = numeric(), GMM_F1 = numeric(), GammaGateR_F1 = numeric(),
  Randkluft_Prec = numeric(), GMM_Prec = numeric(), GammaGateR_Prec = numeric(),
  Randkluft_Rec = numeric(), GMM_Rec = numeric(), GammaGateR_Rec = numeric(),
  Randkluft_Spec = numeric(), GMM_Spec = numeric(), GammaGateR_Spec = numeric(),
  stringsAsFactors = FALSE
)

# Create plots directory
dir.create("plots", showWarnings = FALSE)

# Process each GT pair
for (i in 1:nrow(gate)) {
  patient <- gate$Patient[i]
  marker <- gate$Marker[i]
  gt_gate <- gate$Gate[i]  # Already log-scaled
  
  # Extract marker data for patient
  if (!(marker %in% colnames(raw))) next
  expr <- raw[raw$imageid == patient, marker]
  if (length(expr) < 100) next
  
  # Compute cutoffs (log transform inside functions)
  rand_c <- randkluft_cutoff(expr)
  gmm_c <- gmm_cutoff(expr)
  gamma_c <- gammagater_cutoff(expr)
  
  # GT labels (silver standard)
  gt_labels <- as.integer(log(expr + 1) > gt_gate)
  
  # Predictions
  rand_pred <- as.integer(log(expr + 1) > rand_c)
  gmm_pred <- as.integer(log(expr + 1) > gmm_c)
  gamma_pred <- as.integer(log(expr + 1) > gamma_c)
  
  # Metrics
  rand_met <- compute_metrics(gt_labels, rand_pred)
  gmm_met <- compute_metrics(gt_labels, gmm_pred)
  gamma_met <- compute_metrics(gt_labels, gamma_pred)
  
  # Differences
  diffs <- cutoff_diffs(gt_gate, rand_c, gmm_c, gamma_c)
  
  # Append to results
  results <- rbind(results, data.frame(
    Patient = patient, Marker = marker,
    Randkluft_cutoff = rand_c, GMM_cutoff = gmm_c, GammaGateR_cutoff = gamma_c,
    Randkluft_diff = diffs$Randkluft_diff, GMM_diff = diffs$GMM_diff, GammaGateR_diff = diffs$GammaGateR_diff,
    Randkluft_Acc = rand_met$Accuracy, GMM_Acc = gmm_met$Accuracy, GammaGateR_Acc = gamma_met$Accuracy,
    Randkluft_F1 = rand_met$F1, GMM_F1 = gmm_met$F1, GammaGateR_F1 = gamma_met$F1,
    Randkluft_Prec = rand_met$Precision, GMM_Prec = gmm_met$Precision, GammaGateR_Prec = gamma_met$Precision,
    Randkluft_Rec = rand_met$Recall, GMM_Rec = gmm_met$Recall, GammaGateR_Rec = gamma_met$Recall,
    Randkluft_Spec = rand_met$Specificity, GMM_Spec = gmm_met$Specificity, GammaGateR_Spec = gamma_met$Specificity
  ))
  
  # Histogram plot
  pdf(file = paste0("plots/hist_", patient, "_", marker, ".pdf"))
  hist(log(expr + 1), main = paste(patient, marker), xlab = "log Expression")
  abline(v = gt_gate, col = "black", lwd = 2, lty = 2)
  abline(v = rand_c, col = "red", lwd = 2)
  abline(v = gmm_c, col = "green", lwd = 2)
  abline(v = gamma_c, col = "blue", lwd = 2)
  legend("topright", c("GT", "Randkluft", "GMM", "GammaGateR"), col = c("black", "red", "green", "blue"), lty = c(2,1,1,1))
  dev.off()
}

# Save results
write.csv(results, "/Users/ali/Desktop/Randkluft/gating_comparison_results.csv", row.names = FALSE)


####aggregatiom and visualization 

# Per-marker averages
per_marker <- aggregate(. ~ Marker, data = results[, -1], mean)  # Exclude Patient column

# Overall averages
overall <- colMeans(results[, 5:ncol(results)], na.rm = TRUE)

# 2. Re-compute overall averages safely
library(dplyr)

overall_avg <- results %>%
  summarise(
    Randkluft_Acc = mean(Randkluft_Acc, na.rm = TRUE),
    Randkluft_F1 = mean(Randkluft_F1, na.rm = TRUE),
    Randkluft_Prec = mean(Randkluft_Prec, na.rm = TRUE),
    Randkluft_Rec = mean(Randkluft_Rec, na.rm = TRUE),
    Randkluft_Spec = mean(Randkluft_Spec, na.rm = TRUE),
    
    GMM_Acc = mean(GMM_Acc, na.rm = TRUE),
    GMM_F1 = mean(GMM_F1, na.rm = TRUE),
    GMM_Prec = mean(GMM_Prec, na.rm = TRUE),
    GMM_Rec = mean(GMM_Rec, na.rm = TRUE),
    GMM_Spec = mean(GMM_Spec, na.rm = TRUE),
    
    GammaGateR_Acc = mean(GammaGateR_Acc, na.rm = TRUE),
    GammaGateR_F1 = mean(GammaGateR_F1, na.rm = TRUE),
    GammaGateR_Prec = mean(GammaGateR_Prec, na.rm = TRUE),
    GammaGateR_Rec = mean(GammaGateR_Rec, na.rm = TRUE),
    GammaGateR_Spec = mean(GammaGateR_Spec, na.rm = TRUE)
  )

# Check it worked
print(overall_avg)

# Print
print("Per-Marker Averages:")
print(per_marker)
print("Overall Averages:")
print(overall)

# 3. Now build metrics_long correctly
metrics_long <- data.frame(
  Method = rep(c("Randkluft", "GMM", "GammaGateR"), each = 5),
  Metric = rep(c("Accuracy", "F1", "Precision", "Recall", "Specificity"), 3),
  Value = c(
    overall_avg$Randkluft_Acc, overall_avg$Randkluft_F1, overall_avg$Randkluft_Prec, overall_avg$Randkluft_Rec, overall_avg$Randkluft_Spec,
    overall_avg$GMM_Acc, overall_avg$GMM_F1, overall_avg$GMM_Prec, overall_avg$GMM_Rec, overall_avg$GMM_Spec,
    overall_avg$GammaGateR_Acc, overall_avg$GammaGateR_F1, overall_avg$GammaGateR_Prec, overall_avg$GammaGateR_Rec, overall_avg$GammaGateR_Spec
  )
)

# Verify
str(metrics_long)
head(metrics_long)

# 4. Now plot and save (this will work!)
library(ggplot2)

p_overall <- ggplot(metrics_long, aes(x = Metric, y = Value, fill = Method)) +
  geom_col(position = "dodge", width = 0.7) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Randkluft" = "#e63946", "GMM" = "#2a9d8f", "GammaGateR" = "#457b9d")) +
  labs(title = "Overall Gating Performance Comparison",
       subtitle = paste("Based on", nrow(results), "patient-marker pairs"),
       y = "Average Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# Save and display
p_overall
ggsave("plots/overall_metrics_bar.png", plot = p_overall, width = 10, height = 6, dpi = 150)
print(p_overall)

cat("Plot saved successfully to plots/overall_metrics_bar.png\n")

# ===========================================================================
# NEW SECTION: Visualization of cutoff distance to silver standard (log space)
# ===========================================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# Ensure results has the diff columns (absolute differences)
# They are already in results as Randkluft_diff, GMM_diff, GammaGateR_diff

if (!all(c("Randkluft_diff", "GMM_diff", "GammaGateR_diff") %in% names(results))) {
  stop("Difference columns missing. Make sure cutoff_diffs() is used in the loop.")
}

# -----------------------------
# 1. Per-marker mean absolute error
# -----------------------------
per_marker_diff <- results %>%
  group_by(Marker) %>%
  summarise(
    Randkluft_MAE = mean(Randkluft_diff, na.rm = TRUE),
    GMM_MAE       = mean(GMM_diff, na.rm = TRUE),
    GammaGateR_MAE = mean(GammaGateR_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = ends_with("_MAE"),
    names_to = "Method",
    values_to = "MAE"
  ) %>%
  mutate(
    Method = gsub("_MAE", "", Method),
    Marker = factor(Marker)
  )

# Sort markers by overall MAE for better visual order (optional)
marker_order <- per_marker_diff %>%
  group_by(Marker) %>%
  summarise(mean_MAE = mean(MAE), .groups = "drop") %>%
  arrange(mean_MAE) %>%
  pull(Marker)

per_marker_diff$Marker <- factor(per_marker_diff$Marker, levels = marker_order)

# Compute per-marker MAE (if not already done)
per_marker_diff <- results %>%
  group_by(Marker) %>%
  summarise(
    Randkluft_MAE = mean(Randkluft_diff, na.rm = TRUE),
    GMM_MAE       = mean(GMM_diff, na.rm = TRUE),
    GammaGateR_MAE = mean(GammaGateR_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = ends_with("_MAE"),
    names_to = "Method",
    values_to = "MAE"
  ) %>%
  mutate(Method = gsub("_MAE", "", Method))

# Sort markers by increasing Randkluft MAE
marker_order <- per_marker_diff %>%
  filter(Method == "Randkluft") %>%
  arrange(MAE) %>%
  pull(Marker)

per_marker_diff$Marker <- factor(per_marker_diff$Marker, levels = marker_order)

# Updated plot
p_per_marker_diff_sorted <- ggplot(per_marker_diff, aes(x = Marker, y = MAE, fill = Method)) +
  geom_col(position = "dodge", width = 0.7) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(
    values = c("Randkluft" = "#e63946", "GMM" = "#2a9d8f", "GammaGateR" = "#457b9d")
  ) +
  labs(
    title = "Mean Absolute Error of Predicted Cutoffs vs Silver Standard",
    subtitle = "Markers sorted by increasing Randkluft MAE (log space)",
    y = "Mean |Predicted - GT| (log units)",
    x = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold")
  )

# Save and display
ggsave("plots/per_marker_cutoff_MAE_sorted_by_randkluft.png", p_per_marker_diff_sorted, width = 16, height = 8, dpi = 150)
print(p_per_marker_diff_sorted)
ggsave("plots/per_marker_cutoff_MAE.png", p_per_marker_diff, width = 14, height = 8, dpi = 150)
print(p_per_marker_diff)
cat("Per-marker MAE plot saved to plots/per_marker_cutoff_MAE.png\n")

# Compute differences between Randkluft and GMM MAE per marker
diff_rg <- results %>%
  group_by(Marker) %>%
  summarise(
    Randkluft_MAE = mean(Randkluft_diff, na.rm = TRUE),
    GMM_MAE = mean(GMM_diff, na.rm = TRUE),
    Diff = Randkluft_MAE - GMM_MAE,
    Abs_Diff = abs(Randkluft_MAE - GMM_MAE),
    .groups = "drop"
  ) %>%
  arrange(Abs_Diff)

# 5 most similar (smallest absolute difference)
most_similar <- diff_rg %>% slice_head(n = 5)

# 5 most apart (largest absolute difference)
most_apart <- diff_rg %>% slice_tail(n = 5)

cat("5 markers where Randkluft & GMM are MOST SIMILAR (smallest |MAE diff|):\n")
print(most_similar[, c("Marker", "Diff", "Abs_Diff")])

cat("\n5 markers where Randkluft & GMM are MOST APART (largest |MAE diff|):\n")
print(most_apart[, c("Marker", "Diff", "Abs_Diff")])

# ===========================================================================
# UPDATED 2x5 Histogram Panel: Most Similar (top) vs Most Different (bottom)
# Changes:
# - Narrower line widths (lwd = 1)
# - Added silver standard (GT) as thin silver dashed line
# - Used natural log (log()) instead of log10() for all expressions and cutoffs
# ===========================================================================

library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)

# Define markers
similar_markers <- c("CD57", "PCNA", "IRF1", "cPARP", "TIGIT")
different_markers <- c("CD3e", "CD4", "CD11c", "CD206", "KRT14")
all_markers <- c(similar_markers, different_markers)

# Function to create one histogram (now using natural log)
make_hist_plot <- function(marker_name) {
  # Data for this marker across patients with GT
  marker_data <- raw %>%
    filter(imageid %in% gate$Patient[gate$Marker == marker_name]) %>%
    select(all_of(marker_name)) %>%
    pull()
  
  if (length(marker_data) == 0) {
    return(ggplot() + annotate("text", x = 1, y = 1, label = "No data") + theme_void())
  }
  
  # Use natural log
  ln_expr <- log(marker_data + 1)
  
  # Average cutoffs and GT for this marker
  marker_results <- results %>% filter(Marker == marker_name)
  marker_gt <- gate %>% filter(Marker == marker_name) %>% pull(Gate) %>% mean(na.rm = TRUE)
  
  avg_rand <- mean(marker_results$Randkluft_cutoff, na.rm = TRUE)
  avg_gmm  <- mean(marker_results$GMM_cutoff, na.rm = TRUE)
  
  p <- ggplot(data.frame(ln_expr = ln_expr), aes(x = ln_expr)) +
    geom_histogram(bins = 50, fill = "lightgray", color = "black") +
    theme_minimal(base_size = 12) +
    labs(title = marker_name, x = "ln(Expression + 1)", y = "Count") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10)
    )
  
  # Silver standard (GT) - thin silver dashed line
  if (!is.na(marker_gt)) {
    p <- p + geom_vline(xintercept = marker_gt, color = "grey", lwd = 0.7, linetype = "dashed") +
      annotate("text", x = avg_rand, y = Inf, label = "SG", color = "darkgrey", vjust = 4.5, hjust = -0.1, size = 3)
  }
  
  # Randkluft - solid red, narrow
  if (!is.na(avg_rand)) {
    p <- p + geom_vline(xintercept = avg_rand, color = "#e63946", lwd = 0.7, linetype = "solid") +
      annotate("text", x = avg_rand, y = Inf, label = "Randkluft", color = "#e63946", vjust = 1.5, hjust = -0.1, size = 3)
  }
  
  # GMM - dashed green, narrow
  if (!is.na(avg_gmm)) {
    p <- p + geom_vline(xintercept = avg_gmm, color = "#2a9d8f", lwd = 0.7, linetype = "dashed") +
      annotate("text", x = avg_gmm, y = Inf, label = "GMM", color = "#2a9d8f", vjust = 3, hjust = -0.1, size = 3)
  }
  
  return(p)
}

# Generate all 10 plots
plot_list <- lapply(all_markers, make_hist_plot)

# Row titles
title_top <- textGrob("Most Similar Markers (Randkluft vs GMM)", gp = gpar(fontsize = 16, fontface = "bold"))
title_bottom <- textGrob("Most Different Markers (Randkluft vs GMM)", gp = gpar(fontsize = 16, fontface = "bold"))

# Combine into 2 rows with titles
final_plot <- grid.arrange(
  title_top,
  arrangeGrob(grobs = plot_list[1:5], ncol = 5),
  title_bottom,
  arrangeGrob(grobs = plot_list[6:10], ncol = 5),
  ncol = 1,
  heights = c(0.1, 1, 0.1, 1)
)

# Save
ggsave("plots/2x5_histogram_panel_natural_log_with_GT.png", 
       final_plot, width = 20, height = 12, dpi = 150)

print(final_plot)

cat("Updated 2x5 panel (natural log, narrower lines, silver GT line) saved to plots/2x5_histogram_panel_natural_log_with_GT.png\n")
cat("2x5 histogram panel saved to plots/2x5_histogram_panel_similar_vs_different.png\n")



# -----------------------------
# 2. Overall mean absolute error (across all pairs)
# -----------------------------
overall_diff <- results %>%
  summarise(
    Randkluft_MAE = mean(Randkluft_diff, na.rm = TRUE),
    GMM_MAE       = mean(GMM_diff, na.rm = TRUE),
    GammaGateR_MAE = mean(GammaGateR_diff, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "Method", values_to = "MAE") %>%
  mutate(Method = gsub("_MAE", "", Method))

p_overall_diff <- ggplot(overall_diff, aes(x = Method, y = MAE, fill = Method)) +
  geom_col(width = 0.6) +
  theme_minimal(base_size = 16) +
  scale_fill_manual(
    values = c("Randkluft" = "#e63946", "GMM" = "#2a9d8f", "GammaGateR" = "#457b9d"),
    guide = "none"  # hide legend since x-axis is the method
  ) +
  labs(
    title = "Overall Mean Absolute Error vs Silver Standard",
    subtitle = paste("Across all", nrow(results), "patient-marker pairs (log space)"),
    y = "Mean |Predicted - GT| (log units)",
    x = "Method"
  ) +
  geom_text(aes(label = round(MAE, 3)), vjust = -0.5, size = 6) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 14)
  )
p_overall_diff
ggsave("plots/overall_cutoff_MAE.png", p_overall_diff, width = 9, height = 7, dpi = 150)
print(p_overall_diff)
cat("Overall MAE plot saved to plots/overall_cutoff_MAE.png\n")

# Optional: also compute RMSE if you prefer (uncomment)
# per_marker_rmse <- results %>%
#   group_by(Marker) %>%
#   summarise(
#     Randkluft_RMSE = sqrt(mean(Randkluft_diff^2, na.rm = TRUE)),
#     GMM_RMSE       = sqrt(mean(GMM_diff^2, na.rm = TRUE)),
#     GammaGateR_RMSE = sqrt(mean(GammaGateR_diff^2, na.rm = TRUE))
#   )
# ... similar plotting



