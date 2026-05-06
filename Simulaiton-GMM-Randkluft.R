library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(mixtools)

# Custom Laplace generator
rlaplace <- function(n, location = 0, scale = 1) {
  location + scale * (rexp(n) - rexp(n))
}

# Randkluft cutoff
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

# Gammagater placeholder (since fit_gamma_mixture is unavailable)
gammagater_cutoff <- function(data, log_transform = TRUE) {
  if (log_transform) data <- log(pmax(data + 1, 1))
  warning("fit_gamma_mixture not available; using 95th percentile fallback")
  quantile(data, 0.95)
}

# GMM cutoff with safe log transform
gmm_cutoff <- function(data, log_transform = TRUE) {
  orig_data <- data
  shift <- if (min(data) < 0) abs(min(data)) + 1 else 1
  if (log_transform) data <- log(data + shift)
  
  if (length(unique(data)) < 10 || sd(data) < 1e-6 || is.na(sd(data))) {
    cutoff <- median(data)
  } else {
    best_cutoff <- median(data)
    best_loglik <- -Inf
    
    for (attempt in 1:5) {
      suppressWarnings({
        gmm <- try(mixtools::normalmixEM(data, k = 2, maxit = 1000,
                                         epsilon = 1e-4, verb = FALSE), silent = TRUE)
      })
      if (!inherits(gmm, "try-error") && !is.null(gmm$loglik) && gmm$loglik > best_loglik) {
        best_loglik <- gmm$loglik
        means <- gmm$mu
        sds <- gmm$sigma
        props <- gmm$lambda
        f1 <- function(x) props[1] * dnorm(x, means[1], sds[1])
        f2 <- function(x) props[2] * dnorm(x, means[2], sds[2])
        opt <- optim(mean(means), function(c) abs(f1(c) - f2(c)),
                     method = "Brent", lower = min(means) - 2*max(sds), upper = max(means) + 2*max(sds))
        best_cutoff <- opt$par
      }
    }
    if (best_loglik == -Inf) best_cutoff <- quantile(data, 0.95)
    cutoff <- best_cutoff
  }
  
  if (log_transform) cutoff <- exp(cutoff) - shift
  cutoff
}

# Confusion printer
my_confusion <- function(y_true, y_pred, method) {
  tp <- sum(y_true == 1 & y_pred == 1)
  fp <- sum(y_true == 0 & y_pred == 1)
  tn <- sum(y_true == 0 & y_pred == 0)
  fn <- sum(y_true == 1 & y_pred == 0)
  cat(sprintf("%s: TP=%d, FP=%d, TN=%d, FN=%d\n", method, tp, fp, tn, fn))
}

# Metrics DF
compute_metrics_df <- function(y_true, pred_rand, pred_gamma, pred_gmm, scenario_name) {
  calc <- function(y_pred, method) {
    tp <- sum(y_true == 1 & y_pred == 1)
    tn <- sum(y_true == 0 & y_pred == 0)
    fp <- sum(y_true == 0 & y_pred == 1)
    fn <- sum(y_true == 1 & y_pred == 0)
    acc <- (tp + tn) / length(y_true)
    prec <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    rec <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    f1 <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
    spec <- ifelse(tn + fp > 0, tn / (tn + fp), 0)
    data.frame(Method = method, Accuracy = acc, Precision = prec,
               Recall = rec, F1 = f1, Specificity = spec)
  }
  rbind(calc(pred_rand, "Randkluft"),
        calc(pred_gamma, "Gammagater"),
        calc(pred_gmm, "GMM")) %>%
    mutate(Scenario = scenario_name) %>%
    pivot_longer(-c(Method, Scenario), names_to = "Metric", values_to = "Value") %>%
    mutate(Method = factor(Method, levels = c("Randkluft", "Gammagater", "GMM")),
           Metric = factor(Metric, levels = c("Accuracy", "Specificity", "Precision", "Recall", "F1")))
}

set.seed(42)

n <- 10000
p <- 0.25
nn <- n - n*p
ns <- n*p

scenarios <- list(
  list(name = "Scenario 1: Skewed Rare Positive",
       noise = rnorm(nn, 0, 1) + 3,
       signal = rlnorm(ns, 2, 0.7) + 2 + 3,
       n_noise = nn, n_signal = ns),
  list(name = "Scenario 2: Laplace Noise",
       noise = rlaplace(nn, 3, 1) + 5,
       signal = rexp(ns, 0.2) + 4 + 5,
       n_noise = nn, n_signal = ns),
  list(name = "Scenario 3: Heavy-tailed Noise",
       noise = rt(nn, df = 3) * 1.2 + 40,
       signal =rgamma(ns, 1.8, 0.2) + 6 + 40,
       n_noise = nn, n_signal = ns)
)

results <- lapply(scenarios, function(sc) {
  data <- c(sc$noise, sc$signal)
  labels <- c(rep(0, sc$n_noise), rep(1, sc$n_signal))
  
  c_rand <- randkluft_cutoff(data, log_transform = FALSE)
  c_gamma <- gammagater_cutoff(data, log_transform = FALSE)
  c_gmm <- gmm_cutoff(data, log_transform = FALSE)
  
  pred_rand <- as.integer(data > c_rand)
  pred_gamma <- as.integer(data > c_gamma)
  pred_gmm <- as.integer(data > c_gmm)
  
  cat("=== ", sc$name, " ===\n", sep = "")
  my_confusion(labels, pred_rand, "Randkluft")
  my_confusion(labels, pred_gamma, "Gammagater")
  my_confusion(labels, pred_gmm, "GMM")
  
  list(labels = labels, pred_rand = pred_rand, pred_gamma = pred_gamma, pred_gmm = pred_gmm,
       c_rand = c_rand, c_gamma = c_gamma, c_gmm = c_gmm, data = data, 
       noise = sc$noise, signal = sc$signal, name = sc$name)
})

metrics_all <- do.call(rbind, lapply(results, function(r) {
  compute_metrics_df(r$labels, r$pred_rand, r$pred_gamma, r$pred_gmm, r$name)
}))

# === INTEGRATED 3x2 PLOT: Histograms (left) + Metrics (right) ===
h_list <- lapply(results, function(r) {
  df <- data.frame(
    value = c(r$noise, r$signal),
    group = c(rep("Noise", length(r$noise)), rep("Signal", length(r$signal)))
  )
  cut_df <- data.frame(
    cutoff = c(r$c_rand, r$c_gamma, r$c_gmm),
    Method = c("Randkluft", "Gammagater", "GMM")
  )
  
  p <- ggplot(df, aes(x = value, fill = group)) +
    geom_histogram(data = subset(df, group == "Noise"), bins = 80, alpha = 0.6, linewidth = 0.3) +
    geom_histogram(data = subset(df, group == "Signal"), bins = 80, alpha = 0.6, linewidth = 0.3) +
    geom_vline(data = cut_df, aes(xintercept = cutoff, color = Method),
               linetype = "dashed", size = 1.0) +
    scale_fill_manual(values = c("Noise" = "#a8a8a8", "Signal" = "#e76f51"),
                      guide = "none") +   # Hide fill legend (Noise/Signal)
    scale_color_manual(values = c("Randkluft" = "#e63946", 
                                  "Gammagater" = "#457b9d", 
                                  "GMM" = "#2a9d8f")) +
    labs(title = r$name, x = "Value", y = "Count") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10),
          legend.position = "none")
  
  # Add legend only for cutoffs on the first plot
  if (r$name == results[[1]]$name) {
    p <- p + theme(legend.position = "top") +
      guides(color = guide_legend(title = "Cutoffs", 
                                  override.aes = list(linetype = "dashed")))
  }
  p
})

# Prepare metrics plots
m_list <- lapply(unique(metrics_all$Scenario), function(scen) {
  df <- metrics_all %>% filter(Scenario == scen)
  ggplot(df, aes(x = Metric, y = Value, fill = Method)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("Randkluft" = "#e63946",
                                 "Gammagater" = "#457b9d",
                                 "GMM" = "#2a9d8f")) +
    labs(title = "", y = "Score", x = "") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.margin = margin(10, 10, 10, 10),
          legend.position = "none") +
    ylim(0, 1)
})

# Combine into 3x2 layout
integrated_plot <- grid.arrange(
  grobs = c(h_list[1], m_list[1],
            h_list[2], m_list[2],
            h_list[3], m_list[3]),
  ncol = 2,
  nrow = 3,
  widths = c(1.3, 1),
  top = textGrob("Integrated View: Histograms with Noise/Signal Overlay (left) and Performance Metrics (right)",
                 gp = gpar(fontsize = 18, fontface = "bold")),
  bottom = textGrob("Histograms: Gray = Noise, Orange = Signal (transparent overlay); Dashed lines = detected cutoffs",
                    gp = gpar(fontsize = 11, fontface = "italic"))
)

print(integrated_plot)
cat("Integrated 3x2 figure displayed.\n")