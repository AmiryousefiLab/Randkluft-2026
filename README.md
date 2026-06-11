# Randkluft

**Randkluft** is an R-based method for automated marker intensity gating in multiplexed tissue imaging (e.g., CyCIF). It identifies positive-cell cutoffs by finding the truncation threshold at which the skewness of the lower distribution approaches zero — making no assumptions about the number of mixture components. This repository contains the analysis scripts and reproducibility materials accompanying the paper submitted to *Bioinformatics*.

---

## Method overview

Marker gating (i.e. separating positive from negative cells based on protein expression intensity) is a critical step in single-cell tissue imaging analysis. Existing approaches such as Gaussian Mixture Models (GMM) and GammaGateR require parametric assumptions that can fail for skewed, heavy-tailed, or unimodal distributions commonly encountered in CyCIF data.

Randkluft takes a non-parametric approach: it scans candidate cutoffs on log-transformed intensity data and selects the threshold at which the skewness of the truncated below-cutoff distribution is closest to zero, where the negative population is most symmetric. This is robust to a wide range of distributional shapes and requires no manual initialization.

This Repository contains the code and the outputs relevant to the manuscript for reproducibility. For the main application and the and its complete walk-through guide please consult the publication or visit https://github.com/AmiryousefiLab/Randkluft.

---

## Repository structure

```
Randkluft-2026-main/
├── Simulaiton-GMM-Randkluft.R      # Simulation benchmark (3 scenarios)
├── Silver-standard-compare.R       # Validation against a silver-standard gate dataset
└── Figuers/                        # Pre-generated figures (PDF, PNG, AI formats)
    ├── Fig_1.png / Fig_1.ai        # Main figure
    ├── Fig_1_S.png / ...           # Supplementary figures (S1–S4)
    ├── Simulation.pdf              # Simulation results
    ├── Confusion.pdf               # Confusion matrix summary
    ├── phen.pdf                    # Phenotyping output
    ├── Distance to SG.pdf          # Cutoff distance to silver standard
    ├── all markers.pdf             # Per-marker comparison across methods
    └── 2x5 most similar and least.pdf  # 2×5 histogram panel
```

---

## Scripts

### `Simulaiton-GMM-Randkluft.R`

Benchmarks Randkluft against GMM and GammaGateR under three simulated scenarios designed to stress-test gating methods:

- **Scenario 1: Skewed Rare Positive:** Normal noise + log-normal signal
- **Scenario 2: Laplace Noise:** Laplace-distributed noise + exponential signal
- **Scenario 3: Heavy-tailed Noise:** Student-*t* noise (df = 3) + gamma signal

For each scenario, all three methods compute a cutoff, and classification performance is evaluated using Accuracy, Precision, Recall, F1, and Specificity against known ground-truth labels. Results are displayed as an integrated 3×2 panel (histograms with cutoff overlays on the left; performance bar charts on the right).

**Dependencies:** `dplyr`, `tidyr`, `ggplot2`, `gridExtra`, `grid`, `mixtools`

### `Silver-standard-compare.R`

Validates Randkluft against a hand-annotated silver-standard gate dataset derived from real CyCIF patient data. For each patient–marker pair, the script:

1. Computes Randkluft, GMM, and GammaGateR cutoffs
2. Compares predicted cutoffs to the silver-standard gate (Mean Absolute Error in log space)
3. Computes full classification metrics (Accuracy, Precision, Recall, F1, Specificity) using silver-standard labels as ground truth
4. Produces per-marker and overall MAE plots
5. Generates a 2×5 histogram panel showing the 5 markers where Randkluft and GMM agree most closely and the 5 where they diverge most

**Note:** This script requires two input files not included in this repository due to patient data privacy:
- `raw.csv` — cell-level marker intensity table (columns: `imageid`, marker names)
- `tuulia_data_GT.csv` — silver-standard gate table with columns `Patient`, `Marker`, `Gate` (log-scale)
  (These inputs will be deposited in Zenodo with a permanent DOI upon acceptance of the paper.)

Update the file paths at lines ~118–119 before running.

**Dependencies:** `mixtools`, `ggplot2`, `reshape2`, `dplyr`, `gridExtra`, `grid`  
**Optional:** [`GammaGateR`](https://github.com/JiangmeiRubyXiong/GammaGateR) (install via `devtools::install_github("JiangmeiRubyXiong/GammaGateR")`)

---

## Installation

```r
# Required CRAN packages
install.packages(c("mixtools", "ggplot2", "reshape2", "dplyr",
                   "tidyr", "gridExtra", "grid"))

# Optional: GammaGateR (for silver-standard comparison)
devtools::install_github("JiangmeiRubyXiong/GammaGateR")
```

---

## Usage

### Simulation benchmark

```r
source("Simulaiton-GMM-Randkluft.R")
# Produces an integrated 3×2 figure displayed in the R graphics window
```

### Silver-standard validation

Edit the input paths in `Silver-standard-compare.R`:
```r
raw  <- read.csv("/path/to/raw.csv", stringsAsFactors = FALSE)
gate <- read.csv("/path/to/tuulia_data_GT.csv", stringsAsFactors = FALSE)
```
Then:
```r
source("Silver-standard-compare.R")
# Outputs: gating_comparison_results.csv and plots/ directory
```

---

## Randkluft core function

The gating algorithm is self-contained and can be used independently:

```r
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
    sk <- mean((trunc - m)^3) / sd(trunc)^3
    if (is.finite(sk) && abs(sk) < abs(best_sk)) {
      best_sk <- sk
      best_c <- sorted[i]
    }
  }
  best_c
}

# Example
cutoff <- randkluft_cutoff(my_marker_intensities)
positive_cells <- my_marker_intensities > exp(cutoff) - 1  # back-transform if log_transform = TRUE
```

---

## Interactive app

An interactive Shiny application for exploratory gating with Randkluft is available at:  
🔗 **[irscope.shinyapps.io/Randkluft](https://irscope.shinyapps.io/Randkluft)**

---

## Citation

If you use this code, please cite:

> Amiryousefi A, *et al.* (2026). Randkluft: skewness-minimisation gating for multiplexed tissue imaging. *Bioinformatics* (under review).

---

## License

MIT License. See `LICENSE` for details.

---

## Contact

Ali Amiryousefi, PhD. 
Laboratory of Systems Pharmacology & Ludwig Center, Harvard Medical School  
[ali_amiryousefi@hms.harvard.edu](mailto:ali_amiryousefi@hms.harvard.edu)
