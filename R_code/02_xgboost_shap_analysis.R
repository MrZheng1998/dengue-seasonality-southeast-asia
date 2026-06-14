# =============================================================================
# 02_xgboost_shap_analysis.R
# Main XGBoost-SHAP analysis for national dengue seasonality phenotypes
#
# Input:
#   data/processed/weekly_dengue_processed.csv
#
# Required columns:
#   ISO_A0, week_period, incidence_rate,
#   and all variables listed in FEATURES below.
#
# Output:
#   results/02_xgboost_shap/
#     - model_performance.csv
#     - feature_shap_contribution.csv
#     - group_shap_contribution.csv
#     - feature_group_dictionary.csv
#     - xgb_model_Cluster1.rds
#     - xgb_model_Cluster2.rds
# =============================================================================

suppressPackageStartupMessages({
  library(xgboost)
  library(caret)
  library(dplyr)
  library(purrr)
  library(tibble)
})

set.seed(123)

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

DATA_PATH <- file.path("data", "processed", "weekly_dengue_processed.csv")
OUTPUT_DIR <- file.path("results", "02_xgboost_shap")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SELECTED_COUNTRIES <- c("KHM", "LAO", "MYS", "VNM", "PHL", "SGP", "THA")
ANALYSIS_YEARS <- 2012:2022

FEATURES <- c(
  "annual_rel_amp", "semiannual_rel_amp",
  "tempmax", "tempmin", "temp", "dew", "humidity", "precip", "windspeed",
  "cloudcover", "sealevelpressure", "visibility",
  "GDP", "HEP", "forest_percent", "population_density", "UPR",
  "TP0_14", "TP15_64", "TP65_", "CO2", "Tourist_Arrivals",
  "Dengue.virus.type.1", "Dengue.virus.type.2",
  "Dengue.virus.type.3", "Dengue.virus.type.4"
)

TUNE_GRID <- expand.grid(
  nrounds = c(50, 100, 150),
  eta = c(0.01, 0.05, 0.10),
  max_depth = 3:5,
  gamma = c(0, 0.1, 0.2),
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

CV_CONTROL <- caret::trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  verboseIter = FALSE
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

extract_year <- function(week_period) {
  digits <- gsub("[^0-9]", "", as.character(week_period))
  suppressWarnings(as.integer(substr(digits, 1, 4)))
}

round2 <- function(x) round(as.numeric(x), 2)

get_feature_group <- function(feature) {
  dplyr::case_when(
    feature %in% c("annual_rel_amp", "semiannual_rel_amp") ~ "Seasonal Amplitude",
    feature %in% c(
      "tempmax", "tempmin", "temp", "dew", "humidity", "precip",
      "windspeed", "cloudcover", "sealevelpressure", "visibility"
    ) ~ "Climate",
    feature %in% c("GDP", "HEP", "Tourist_Arrivals", "UPR", "forest_percent", "CO2") ~ "Social Economy",
    feature %in% c("population_density", "TP0_14", "TP15_64", "TP65_") ~ "Population",
    feature %in% c(
      "Dengue.virus.type.1", "Dengue.virus.type.2",
      "Dengue.virus.type.3", "Dengue.virus.type.4"
    ) ~ "Serotype Distribution",
    TRUE ~ "Other"
  )
}

add_cluster_labels <- function(dat) {
  dat %>%
    mutate(
      Cluster = case_when(
        ISO_A0 %in% c("VNM", "LAO", "KHM", "PHL", "THA") ~ "Cluster1",
        ISO_A0 %in% c("SGP", "MYS") ~ "Cluster2",
        TRUE ~ NA_character_
      ),
      Cluster = factor(Cluster, levels = c("Cluster1", "Cluster2"))
    ) %>%
    filter(!is.na(Cluster))
}

read_analysis_data <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }

  dat <- read.csv(path, stringsAsFactors = FALSE)

  required_cols <- c("ISO_A0", "week_period", "incidence_rate", FEATURES)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      Year = extract_year(week_period),
      incidence_rate = suppressWarnings(as.numeric(incidence_rate))
    ) %>%
    filter(
      ISO_A0 %in% SELECTED_COUNTRIES,
      Year %in% ANALYSIS_YEARS,
      !is.na(incidence_rate)
    ) %>%
    add_cluster_labels()
}

prepare_cluster_data <- function(cluster_data) {
  df <- cluster_data %>%
    select(incidence_rate, all_of(FEATURES), ISO_A0, Year, week_period, Cluster)

  for (feature in FEATURES) {
    df[[feature]] <- suppressWarnings(as.numeric(df[[feature]]))
    median_value <- median(df[[feature]], na.rm = TRUE)
    if (!is.finite(median_value)) median_value <- 0
    df[[feature]][!is.finite(df[[feature]]) | is.na(df[[feature]])] <- median_value
  }

  x <- data.matrix(df[, FEATURES, drop = FALSE])
  y <- as.numeric(df$incidence_rate)

  list(
    df = df,
    x = x,
    y = y,
    dmatrix = xgboost::xgb.DMatrix(data = x, label = y),
    features = colnames(x)
  )
}

fit_xgb_shap <- function(cluster_data, cluster_name) {
  message("Training XGBoost model for ", cluster_name, " ...")

  processed <- prepare_cluster_data(cluster_data)

  set.seed(42)
  tuned_model <- caret::train(
    x = processed$x,
    y = processed$y,
    method = "xgbTree",
    trControl = CV_CONTROL,
    tuneGrid = TUNE_GRID,
    metric = "RMSE",
    verbosity = 0
  )

  best <- as.list(tuned_model$bestTune)

  final_params <- list(
    objective = "reg:squarederror",
    booster = "gbtree",
    eval_metric = "rmse",
    eta = best$eta,
    max_depth = best$max_depth,
    gamma = best$gamma,
    subsample = best$subsample,
    colsample_bytree = best$colsample_bytree,
    min_child_weight = best$min_child_weight
  )

  final_model <- xgboost::xgb.train(
    params = final_params,
    data = processed$dmatrix,
    nrounds = best$nrounds,
    verbose = 0
  )

  predictions <- predict(final_model, processed$dmatrix)
  residuals <- processed$y - predictions
  sse <- sum(residuals^2, na.rm = TRUE)
  sst <- sum((processed$y - mean(processed$y, na.rm = TRUE))^2, na.rm = TRUE)

  performance <- tibble(
    Cluster = cluster_name,
    n = length(processed$y),
    RMSE = sqrt(mean(residuals^2, na.rm = TRUE)),
    MAE = mean(abs(residuals), na.rm = TRUE),
    R2 = ifelse(sst > 0, 1 - sse / sst, NA_real_),
    best_nrounds = best$nrounds,
    best_eta = best$eta,
    best_max_depth = best$max_depth,
    best_gamma = best$gamma
  )

  shap_contrib <- as.data.frame(
    predict(final_model, processed$dmatrix, predcontrib = TRUE)
  )

  bias_cols <- grep("BIAS|bias", names(shap_contrib), value = TRUE)
  shap_matrix <- shap_contrib[, setdiff(names(shap_contrib), bias_cols), drop = FALSE]
  shap_matrix <- shap_matrix[, processed$features, drop = FALSE]

  mean_abs_shap <- colMeans(abs(as.matrix(shap_matrix)), na.rm = TRUE)
  total_abs_shap <- sum(mean_abs_shap, na.rm = TRUE)

  feature_contribution <- tibble(
    Cluster = cluster_name,
    Feature = names(mean_abs_shap),
    Mean_ABS_SHAP = as.numeric(mean_abs_shap),
    Contribution = ifelse(total_abs_shap > 0, Mean_ABS_SHAP / total_abs_shap * 100, NA_real_),
    Group = get_feature_group(Feature)
  ) %>%
    arrange(desc(Contribution))

  group_contribution <- feature_contribution %>%
    group_by(Cluster, Group) %>%
    summarise(Total_Contribution = sum(Contribution, na.rm = TRUE), .groups = "drop") %>%
    group_by(Cluster) %>%
    mutate(Total_Contribution = Total_Contribution / sum(Total_Contribution, na.rm = TRUE) * 100) %>%
    ungroup() %>%
    arrange(Cluster, desc(Total_Contribution))

  saveRDS(final_model, file.path(OUTPUT_DIR, paste0("xgb_model_", cluster_name, ".rds")))

  list(
    performance = performance,
    feature_contribution = feature_contribution,
    group_contribution = group_contribution
  )
}

# -----------------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------------

main <- function() {
  dat <- read_analysis_data(DATA_PATH)

  write.csv(
    dat %>% count(Cluster, ISO_A0, Year),
    file.path(OUTPUT_DIR, "data_count_by_cluster_country_year.csv"),
    row.names = FALSE
  )

  cluster_data <- split(dat, dat$Cluster)
  cluster_data <- cluster_data[c("Cluster1", "Cluster2")]
  cluster_data <- cluster_data[!vapply(cluster_data, is.null, logical(1))]

  results <- purrr::imap(cluster_data, fit_xgb_shap)

  performance <- purrr::map_dfr(results, "performance")
  feature_contribution <- purrr::map_dfr(results, "feature_contribution")
  group_contribution <- purrr::map_dfr(results, "group_contribution")

  write.csv(
    performance %>% mutate(across(c(RMSE, MAE, R2), round2)),
    file.path(OUTPUT_DIR, "model_performance.csv"),
    row.names = FALSE
  )

  write.csv(
    feature_contribution %>% mutate(across(c(Mean_ABS_SHAP, Contribution), round2)),
    file.path(OUTPUT_DIR, "feature_shap_contribution.csv"),
    row.names = FALSE
  )

  write.csv(
    group_contribution %>% mutate(Total_Contribution = round2(Total_Contribution)),
    file.path(OUTPUT_DIR, "group_shap_contribution.csv"),
    row.names = FALSE
  )

  write.csv(
    tibble(Feature = FEATURES, Group = get_feature_group(FEATURES)),
    file.path(OUTPUT_DIR, "feature_group_dictionary.csv"),
    row.names = FALSE
  )

  message("XGBoost-SHAP analysis completed. Output directory: ", normalizePath(OUTPUT_DIR))
}

main()
