# =============================================================================
# GAM threshold analysis for SHAP-derived feature contributions
#
# Purpose:
#   1. Fit cluster-specific GAM curves for SHAP values against raw feature values.
#   2. Identify zero-crossing points where fitted SHAP contributions change sign.
#   3. Export a publication-ready threshold table and a GAM curve figure.
#
# Input:
#   data/processed/shap_feature_values.csv
#
# Required columns:
#   Cluster, Feature, Raw, SHAP
#
# Output:
#   results/03_gam_threshold/
# =============================================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
})

set.seed(123)

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

INPUT_FILE <- "data/processed/shap_feature_values.csv"
OUTPUT_DIR <- "results/03_gam_threshold"

# Edit this vector if a different set of predictors is used in the manuscript.
FEATURES_TO_ANALYZE <- c("precip", "humidity")

MIN_OBSERVATIONS <- 20
GAM_BASIS_K <- 8
PREDICTION_GRID_N <- 300
OUTLIER_LOWER_Q <- 0.01
OUTLIER_UPPER_Q <- 0.99
CONFIDENCE_LEVEL <- 0.95

# Optional plotting labels.
FEATURE_LABELS <- c(
  precip = "Precipitation",
  tempmin = "Minimum temperature"
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

create_output_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

validate_input <- function(df) {
  required_cols <- c("Cluster", "Feature", "Raw", "SHAP")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Input file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

read_shap_data <- function(input_file) {
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file,
         "\nPlease place a processed SHAP dataset at this path or update INPUT_FILE.")
  }

  df <- readr::read_csv(input_file, show_col_types = FALSE)
  validate_input(df)

  df %>%
    mutate(
      Cluster = as.character(Cluster),
      Feature = as.character(Feature),
      Raw = suppressWarnings(as.numeric(Raw)),
      SHAP = suppressWarnings(as.numeric(SHAP))
    ) %>%
    filter(
      Feature %in% FEATURES_TO_ANALYZE,
      is.finite(Raw),
      is.finite(SHAP)
    )
}

remove_outliers_by_group <- function(df,
                                     lower_q = OUTLIER_LOWER_Q,
                                     upper_q = OUTLIER_UPPER_Q) {
  df %>%
    group_by(Cluster, Feature) %>%
    mutate(
      raw_lower = quantile(Raw, lower_q, na.rm = TRUE),
      raw_upper = quantile(Raw, upper_q, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(Raw >= raw_lower, Raw <= raw_upper) %>%
    select(-raw_lower, -raw_upper)
}

prepare_gam_data <- function(df) {
  df %>%
    remove_outliers_by_group() %>%
    mutate(
      signed_log_shap = sign(SHAP) * log1p(abs(SHAP)),
      Feature_label = ifelse(
        Feature %in% names(FEATURE_LABELS),
        FEATURE_LABELS[Feature],
        Feature
      )
    )
}

fit_one_gam <- function(df) {
  if (nrow(df) < MIN_OBSERVATIONS || length(unique(df$Raw)) < 5) {
    return(NULL)
  }

  model <- mgcv::gam(
    signed_log_shap ~ s(Raw, k = min(GAM_BASIS_K, length(unique(df$Raw)) - 1)),
    data = df,
    method = "REML"
  )

  raw_grid <- seq(min(df$Raw, na.rm = TRUE),
                  max(df$Raw, na.rm = TRUE),
                  length.out = PREDICTION_GRID_N)

  pred <- predict(model, newdata = data.frame(Raw = raw_grid), se.fit = TRUE)
  model_summary <- summary(model)

  tibble(
    Cluster = unique(df$Cluster),
    Feature = unique(df$Feature),
    Feature_label = unique(df$Feature_label),
    Raw = raw_grid,
    fit = as.numeric(pred$fit),
    se = as.numeric(pred$se.fit),
    lower = fit - qnorm(1 - (1 - CONFIDENCE_LEVEL) / 2) * se,
    upper = fit + qnorm(1 - (1 - CONFIDENCE_LEVEL) / 2) * se,
    adjusted_r2 = as.numeric(model_summary$r.sq),
    smooth_p_value = as.numeric(model_summary$s.table[1, "p-value"]),
    n = nrow(df)
  )
}

fit_gam_by_cluster_feature <- function(df) {
  df %>%
    group_by(Cluster, Feature) %>%
    group_split() %>%
    purrr::map_dfr(fit_one_gam)
}

find_zero_crossings <- function(pred_df) {
  z <- qnorm(1 - (1 - CONFIDENCE_LEVEL) / 2)
  pred_df <- pred_df %>% arrange(Raw)

  if (nrow(pred_df) < 2) return(tibble())

  crossing_index <- which(pred_df$fit[-nrow(pred_df)] * pred_df$fit[-1] < 0)
  if (length(crossing_index) == 0) return(tibble())

  purrr::map_dfr(crossing_index, function(i) {
    x1 <- pred_df$Raw[i]
    x2 <- pred_df$Raw[i + 1]
    y1 <- pred_df$fit[i]
    y2 <- pred_df$fit[i + 1]
    se1 <- pred_df$se[i]
    se2 <- pred_df$se[i + 1]

    zero_x <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
    weight <- (zero_x - x1) / (x2 - x1)
    zero_se <- se1 * (1 - weight) + se2 * weight

    tibble(
      zero_crossing = zero_x,
      ci_lower = zero_x - z * zero_se,
      ci_upper = zero_x + z * zero_se
    )
  })
}

summarize_zero_crossings <- function(gam_pred_df) {
  gam_pred_df %>%
    group_by(Feature, Feature_label, Cluster) %>%
    group_modify(~ find_zero_crossings(.x)) %>%
    ungroup() %>%
    mutate(
      zero_crossing = round(zero_crossing, 2),
      ci_lower = round(ci_lower, 2),
      ci_upper = round(ci_upper, 2),
      threshold_95ci = paste0(zero_crossing, " (", ci_lower, ", ", ci_upper, ")")
    ) %>%
    arrange(Feature, Cluster)
}

plot_gam_curves <- function(df_clean, gam_pred_df, output_file) {
  p <- ggplot() +
    geom_point(
      data = df_clean,
      aes(x = Raw, y = signed_log_shap, color = Cluster),
      alpha = 0.35,
      size = 0.8
    ) +
    geom_ribbon(
      data = gam_pred_df,
      aes(x = Raw, ymin = lower, ymax = upper, fill = Cluster),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(
      data = gam_pred_df,
      aes(x = Raw, y = fit, color = Cluster),
      linewidth = 0.9
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ Feature_label, scales = "free_x") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = "Raw feature value",
      y = "Signed log-transformed SHAP value",
      color = "Cluster",
      fill = "Cluster"
    )

  ggsave(output_file, p, width = 9, height = 5)
  return(p)
}

# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------

create_output_dir(OUTPUT_DIR)

shap_raw <- read_shap_data(INPUT_FILE)
shap_clean <- prepare_gam_data(shap_raw)

gam_predictions <- fit_gam_by_cluster_feature(shap_clean)

if (nrow(gam_predictions) == 0) {
  stop("No GAM models were fitted. Check sample size and selected features.")
}

zero_crossings <- summarize_zero_crossings(gam_predictions)

model_summary <- gam_predictions %>%
  distinct(Feature, Feature_label, Cluster, n, adjusted_r2, smooth_p_value) %>%
  mutate(
    adjusted_r2 = round(adjusted_r2, 3),
    smooth_p_value = signif(smooth_p_value, 3)
  ) %>%
  arrange(Feature, Cluster)

readr::write_csv(shap_clean, file.path(OUTPUT_DIR, "gam_input_data.csv"))
readr::write_csv(gam_predictions, file.path(OUTPUT_DIR, "gam_predictions.csv"))
readr::write_csv(zero_crossings, file.path(OUTPUT_DIR, "gam_zero_crossing_thresholds.csv"))
readr::write_csv(model_summary, file.path(OUTPUT_DIR, "gam_model_summary.csv"))

plot_gam_curves(
  df_clean = shap_clean,
  gam_pred_df = gam_predictions,
  output_file = file.path(OUTPUT_DIR, "gam_shap_curves.svg")
)

cat("GAM threshold analysis completed.\n")
cat("Output directory:", normalizePath(OUTPUT_DIR), "\n")
