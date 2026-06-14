# =============================================================================
# 01_dtw_clustering.R
# Country-level dengue seasonality clustering and COVID-period sensitivity analysis
#
# Purpose:
#   1. Build country-level 52-week seasonal profiles from weekly dengue data.
#   2. Perform DTW distance-based hierarchical clustering.
#   3. Evaluate k = 2 to 6 and report k = 2, k = 3, and k = 4 memberships.
#   4. Assess clustering and summer-autumn peak timing stability after excluding
#      COVID-affected years.
#
# Required input columns:
#   ISO_A0        Country ISO3 code, e.g., SGP, MYS, THA
#   week_period   Epidemiological week identifier containing YYYYWW or YYYY-WW
#   dengue_total  Weekly reported dengue cases
#
# Expected input file:
#   data/processed/merged_data_final.csv
#
# Main outputs:
#   results/dtw_clustering/<analysis_window>/cluster_assignments_k*.csv
#   results/dtw_clustering/<analysis_window>/cluster_evaluation_k2_to_k6.csv
#   results/dtw_clustering/comparison_tables/Table_cluster_number_evaluation_for_reviewer.csv
#   results/dtw_clustering/comparison_tables/cluster_membership_comparison_k2.csv
#   results/dtw_clustering/comparison_tables/country_peak_timing_comparison.csv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(dtwclust)
  library(dtw)
  library(proxy)
  library(cluster)
  library(mclust)
  library(zoo)
})

set.seed(123)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

DATA_PATH <- file.path("data", "processed", "merged_data_final.csv")
OUTPUT_ROOT <- file.path("results", "dtw_clustering")
dir.create(OUTPUT_ROOT, showWarnings = FALSE, recursive = TRUE)

SELECTED_COUNTRIES <- c("KHM", "LAO", "MYS", "VNM", "PHL", "SGP", "THA")

# Approximate country centroid latitudes used only for ordering and peak-latitude summaries.
# These values do not define cluster membership.
COUNTRY_LAT <- tibble::tribble(
  ~ISO_A0, ~Latitude,
  "SGP",   1.35,
  "MYS",   4.21,
  "KHM",  12.57,
  "PHL",  12.88,
  "VNM",  14.06,
  "THA",  15.87,
  "LAO",  19.86
) %>%
  filter(ISO_A0 %in% SELECTED_COUNTRIES)

MAIN_DTW_WINDOW <- 4
SENSITIVITY_WINDOWS <- 3:7
NORMALIZATION_METHODS <- c("proportional", "minmax", "rank")
K_RANGE <- 2:6
N_BOOT <- 500
SUMMER_AUTUMN_WEEKS <- 19:44

ANALYSIS_WINDOWS <- list(
  full_2012_2022 = 2012:2022,
  exclude_2020_2021 = setdiff(2012:2022, 2020:2021),
  pre_covid_2012_2019 = 2012:2019
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

safe_mean <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_impute_series <- function(x, renormalize = TRUE) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_

  if (all(is.na(x))) return(rep(0, length(x)))

  if (sum(!is.na(x)) == 1) {
    x2 <- rep(x[which(!is.na(x))[1]], length(x))
  } else {
    x2 <- zoo::na.approx(x, x = seq_along(x), na.rm = FALSE, rule = 2)
  }

  x2[!is.finite(x2)] <- 0

  if (renormalize) {
    s <- sum(x2, na.rm = TRUE)
    if (is.finite(s) && s > 0) x2 <- x2 / s
  }

  x2
}

week_to_month <- function(week) {
  dplyr::case_when(
    week >= 1  & week <= 4  ~ 1,
    week >= 5  & week <= 8  ~ 2,
    week >= 9  & week <= 13 ~ 3,
    week >= 14 & week <= 17 ~ 4,
    week >= 18 & week <= 21 ~ 5,
    week >= 22 & week <= 26 ~ 6,
    week >= 27 & week <= 30 ~ 7,
    week >= 31 & week <= 35 ~ 8,
    week >= 36 & week <= 39 ~ 9,
    week >= 40 & week <= 44 ~ 10,
    week >= 45 & week <= 47 ~ 11,
    week >= 48 & week <= 52 ~ 12,
    TRUE ~ NA_real_
  )
}

normalize_matrix <- function(mat, method = c("proportional", "minmax", "rank")) {
  method <- match.arg(method)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  if (method == "proportional") return(mat)

  out <- t(apply(mat, 1, function(x) {
    if (method == "minmax") {
      rng <- range(x, na.rm = TRUE)
      if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
      return((x - rng[1]) / diff(rng))
    }

    if (all(x == x[1])) return(rep(0.5, length(x)))
    rank(x, ties.method = "average") / length(x)
  }))

  rownames(out) <- rownames(mat)
  out
}

safe_norm01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  if (!is.finite(xmin) || !is.finite(xmax) || xmax == xmin) return(rep(1, length(x)))
  (x - xmin) / (xmax - xmin)
}

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

read_weekly_data <- function(data_path = DATA_PATH) {
  if (!file.exists(data_path)) {
    stop("Input file not found: ", data_path)
  }

  dat <- readr::read_csv(data_path, show_col_types = FALSE)
  required_cols <- c("ISO_A0", "week_period", "dengue_total")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Input data are missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  week_period_digits <- gsub("[^0-9]", "", as.character(dat$week_period))

  dat %>%
    mutate(
      year = suppressWarnings(as.integer(substr(week_period_digits, 1, 4))),
      week = suppressWarnings(as.integer(substr(week_period_digits, 5, 6))),
      dengue_total = suppressWarnings(as.numeric(dengue_total))
    ) %>%
    filter(
      ISO_A0 %in% SELECTED_COUNTRIES,
      !is.na(year),
      !is.na(week),
      week >= 1,
      week <= 52
    ) %>%
    arrange(ISO_A0, year, week)
}

prepare_country_seasonal_matrix <- function(dat, keep_years) {
  d <- dat %>%
    filter(year %in% keep_years) %>%
    mutate(dengue_total = replace_na(as.numeric(dengue_total), 0))

  if (nrow(d) == 0) {
    stop("No records found for the requested analysis window.")
  }

  # Country-year weekly proportions remove differences in annual incidence magnitude
  # and retain the within-year seasonal profile.
  d_prop <- d %>%
    group_by(ISO_A0, year) %>%
    mutate(
      annual_total = sum(dengue_total, na.rm = TRUE),
      incidence_prop = if_else(annual_total > 0, dengue_total / annual_total, NA_real_)
    ) %>%
    ungroup()

  country_week <- d_prop %>%
    group_by(ISO_A0, week) %>%
    summarise(mean_incidence = safe_mean(incidence_prop), .groups = "drop") %>%
    complete(ISO_A0 = SELECTED_COUNTRIES, week = 1:52) %>%
    arrange(ISO_A0, week) %>%
    group_by(ISO_A0) %>%
    mutate(mean_incidence = safe_impute_series(mean_incidence, renormalize = TRUE)) %>%
    ungroup()

  mat_df <- country_week %>%
    pivot_wider(names_from = week, values_from = mean_incidence, names_prefix = "W") %>%
    arrange(match(ISO_A0, SELECTED_COUNTRIES))

  mat <- as.matrix(select(mat_df, -ISO_A0))
  storage.mode(mat) <- "numeric"
  rownames(mat) <- mat_df$ISO_A0

  if (any(!is.finite(mat))) {
    stop("The DTW input matrix contains non-finite values.")
  }

  list(matrix = mat, long = country_week, mat_df = mat_df)
}

# -----------------------------------------------------------------------------
# DTW clustering
# -----------------------------------------------------------------------------

weighted_dtw <- function(x, y, window_size = MAIN_DTW_WINDOW) {
  dtwclust::dtw_basic(
    x, y,
    window.size = window_size,
    window.type = "slantedband",
    step.pattern = dtw::symmetric2
  )
}

calc_dtw_dist <- function(mat, window_size = MAIN_DTW_WINDOW) {
  proxy::dist(
    as.matrix(mat),
    method = function(x, y) weighted_dtw(x, y, window_size)
  )
}

evaluate_k_values <- function(mat, dtw_dist, keep_years, dat, n_boot = N_BOOT) {
  hc <- hclust(dtw_dist, method = "ward.D2")

  eval_tbl <- purrr::map_dfr(K_RANGE, function(k) {
    cl <- cutree(hc, k = k)
    sil <- cluster::silhouette(cl, dtw_dist)
    cluster_sizes <- table(cl)

    tibble(
      k = k,
      silhouette = mean(sil[, "sil_width"], na.rm = TRUE),
      min_cluster_size = min(as.integer(cluster_sizes)),
      n_singleton_clusters = sum(as.integer(cluster_sizes) == 1),
      cluster_size_string = paste(
        paste0("C", names(cluster_sizes), "=", as.integer(cluster_sizes)),
        collapse = "; "
      )
    )
  })

  # Gap statistic is used as a supplementary diagnostic.
  gap_statistic <- rep(NA_real_, length(K_RANGE))
  gap_try <- try({
    gap_obj <- cluster::clusGap(
      as.matrix(dtw_dist),
      FUN = function(x, k) {
        hc_x <- hclust(as.dist(x), method = "ward.D2")
        list(cluster = cutree(hc_x, k = k))
      },
      K.max = max(K_RANGE),
      B = 100
    )
    gap_statistic <- gap_obj$Tab[K_RANGE, "gap"]
  }, silent = TRUE)
  eval_tbl$gap_statistic <- gap_statistic

  # Year-bootstrap stability assesses whether cluster solutions are driven by a small
  # number of years. This is separate from cross-window COVID sensitivity below.
  all_years <- sort(unique(keep_years))
  boot_tbl <- purrr::map_dfr(K_RANGE, function(k) {
    original_cl <- cutree(hc, k = k)
    ari_vec <- rep(NA_real_, n_boot)

    for (b in seq_len(n_boot)) {
      sampled_years <- sample(all_years, size = length(all_years), replace = TRUE)
      prep_b <- tryCatch(prepare_country_seasonal_matrix(dat, sampled_years), error = function(e) NULL)
      if (is.null(prep_b)) next

      dist_b <- tryCatch(calc_dtw_dist(prep_b$matrix, MAIN_DTW_WINDOW), error = function(e) NULL)
      if (is.null(dist_b)) next

      cl_b <- tryCatch(cutree(hclust(dist_b, method = "ward.D2"), k = k), error = function(e) NULL)
      if (!is.null(cl_b)) ari_vec[b] <- mclust::adjustedRandIndex(original_cl, cl_b)
    }

    tibble(
      k = k,
      bootstrap_ARI_mean = mean(ari_vec, na.rm = TRUE),
      bootstrap_ARI_sd = sd(ari_vec, na.rm = TRUE),
      bootstrap_valid_n = sum(!is.na(ari_vec))
    )
  })

  eval_tbl %>%
    left_join(boot_tbl, by = "k") %>%
    mutate(
      silhouette_norm = safe_norm01(silhouette),
      gap_norm = safe_norm01(gap_statistic),
      bootstrap_norm = safe_norm01(bootstrap_ARI_mean),
      combined_score = rowMeans(
        cbind(silhouette_norm, gap_norm, bootstrap_norm),
        na.rm = TRUE
      )
    )
}

save_cluster_assignments <- function(dtw_dist, outdir, analysis_name) {
  hc <- hclust(dtw_dist, method = "ward.D2")
  hc$labels <- labels(dtw_dist)

  cluster_outputs <- list()
  for (k in K_RANGE) {
    cl <- cutree(hc, k = k)
    tbl <- tibble(
      analysis_name = analysis_name,
      ISO_A0 = names(cl),
      k = k,
      Cluster = as.integer(cl)
    )

    readr::write_csv(tbl, file.path(outdir, paste0("cluster_assignments_k", k, ".csv")))
    cluster_outputs[[paste0("k", k)]] <- tbl
  }
  cluster_outputs
}

run_parameter_sensitivity <- function(mat, main_clusters_k2) {
  expand.grid(
    window_size = SENSITIVITY_WINDOWS,
    normalization = NORMALIZATION_METHODS,
    stringsAsFactors = FALSE
  ) %>%
    purrr::pmap_dfr(function(window_size, normalization) {
      test_mat <- normalize_matrix(mat, method = normalization)
      test_dist <- tryCatch(calc_dtw_dist(test_mat, window_size), error = function(e) NULL)

      if (is.null(test_dist)) {
        return(tibble(
          window_size = window_size,
          normalization = normalization,
          avg_silhouette = NA_real_,
          ARI_with_main_k2 = NA_real_
        ))
      }

      test_cl <- cutree(hclust(test_dist, method = "ward.D2"), k = 2)
      sil <- cluster::silhouette(test_cl, test_dist)

      tibble(
        window_size = window_size,
        normalization = normalization,
        avg_silhouette = mean(sil[, "sil_width"], na.rm = TRUE),
        ARI_with_main_k2 = mclust::adjustedRandIndex(main_clusters_k2, test_cl)
      )
    })
}

# -----------------------------------------------------------------------------
# Peak timing analysis
# -----------------------------------------------------------------------------

detect_country_peak_timing <- function(country_week, cluster_k2, analysis_name) {
  country_week %>%
    filter(week %in% SUMMER_AUTUMN_WEEKS) %>%
    group_by(ISO_A0) %>%
    slice_max(order_by = mean_incidence, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      analysis_name = analysis_name,
      ISO_A0,
      peak_week = week,
      peak_month = week_to_month(week),
      peak_intensity = mean_incidence
    ) %>%
    left_join(select(cluster_k2, ISO_A0, Cluster), by = "ISO_A0") %>%
    left_join(COUNTRY_LAT, by = "ISO_A0")
}

summarize_cluster_peak_timing <- function(peak_tbl) {
  peak_tbl %>%
    group_by(analysis_name, Cluster) %>%
    summarise(
      n_countries = n(),
      mean_peak_week = mean(peak_week, na.rm = TRUE),
      sd_peak_week = sd(peak_week, na.rm = TRUE),
      median_peak_week = median(peak_week, na.rm = TRUE),
      mean_peak_month = mean(peak_month, na.rm = TRUE),
      countries = paste(ISO_A0, collapse = ", "),
      .groups = "drop"
    )
}

# -----------------------------------------------------------------------------
# Main analysis functions
# -----------------------------------------------------------------------------

run_one_analysis <- function(dat, analysis_name, keep_years) {
  message("Running DTW analysis: ", analysis_name)

  outdir <- file.path(OUTPUT_ROOT, analysis_name)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  prep <- prepare_country_seasonal_matrix(dat, keep_years)
  dtw_dist <- calc_dtw_dist(prep$matrix, MAIN_DTW_WINDOW)

  readr::write_csv(prep$mat_df, file.path(outdir, "country_seasonal_matrix_52weeks.csv"))
  readr::write_csv(prep$long, file.path(outdir, "country_weekly_mean_profile_long.csv"))
  write.csv(as.matrix(dtw_dist), file.path(outdir, "dtw_distance_matrix.csv"))

  eval_tbl <- evaluate_k_values(prep$matrix, dtw_dist, keep_years, dat, N_BOOT)
  readr::write_csv(eval_tbl, file.path(outdir, "cluster_evaluation_k2_to_k6.csv"))

  cluster_outputs <- save_cluster_assignments(dtw_dist, outdir, analysis_name)

  main_k2 <- cluster_outputs$k2$Cluster
  names(main_k2) <- cluster_outputs$k2$ISO_A0

  sens_tbl <- run_parameter_sensitivity(prep$matrix, main_k2)
  readr::write_csv(sens_tbl, file.path(outdir, "dtw_parameter_sensitivity_k2.csv"))

  peak_tbl <- detect_country_peak_timing(prep$long, cluster_outputs$k2, analysis_name)
  cluster_peak_summary <- summarize_cluster_peak_timing(peak_tbl)

  readr::write_csv(peak_tbl, file.path(outdir, "country_summer_autumn_peak_timing.csv"))
  readr::write_csv(cluster_peak_summary, file.path(outdir, "cluster_peak_timing_summary_k2.csv"))

  list(
    analysis_name = analysis_name,
    keep_years = keep_years,
    eval_tbl = eval_tbl,
    sens_tbl = sens_tbl,
    cluster_outputs = cluster_outputs,
    peak_tbl = peak_tbl,
    cluster_peak_summary = cluster_peak_summary
  )
}

get_cluster_vector <- function(result_obj, k) {
  tbl <- result_obj$cluster_outputs[[paste0("k", k)]]
  v <- tbl$Cluster
  names(v) <- tbl$ISO_A0
  v
}

calc_cross_window_ari <- function(all_results, main_name = "full_2012_2022") {
  purrr::map_dfr(names(all_results), function(aname) {
    purrr::map_dfr(K_RANGE, function(k) {
      main_vec <- get_cluster_vector(all_results[[main_name]], k)
      test_vec <- get_cluster_vector(all_results[[aname]], k)
      common <- intersect(names(main_vec), names(test_vec))

      tibble(
        analysis_name = aname,
        k = k,
        cross_window_ARI_vs_full = mclust::adjustedRandIndex(main_vec[common], test_vec[common]),
        SGP_MYS_grouped_together = unname(test_vec["SGP"] == test_vec["MYS"])
      )
    })
  })
}

write_comparison_outputs <- function(all_results) {
  comparison_dir <- file.path(OUTPUT_ROOT, "comparison_tables")
  dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)

  for (k in c(2, 3, 4)) {
    membership_tbl <- purrr::map_dfr(all_results, ~ .x$cluster_outputs[[paste0("k", k)]]) %>%
      select(analysis_name, ISO_A0, Cluster) %>%
      pivot_wider(names_from = analysis_name, values_from = Cluster)

    readr::write_csv(
      membership_tbl,
      file.path(comparison_dir, paste0("cluster_membership_comparison_k", k, ".csv"))
    )
  }

  cross_window_tbl <- calc_cross_window_ari(all_results)

  eval_compare <- purrr::map_dfr(all_results, function(x) {
    x$eval_tbl %>% mutate(analysis_name = x$analysis_name)
  }) %>%
    left_join(cross_window_tbl, by = c("analysis_name", "k")) %>%
    arrange(analysis_name, k)

  readr::write_csv(
    eval_compare,
    file.path(comparison_dir, "cluster_evaluation_comparison_k2_to_k6.csv")
  )

  reviewer_tbl <- eval_compare %>%
    transmute(
      Analysis_window = analysis_name,
      k = k,
      Average_silhouette_width = round(silhouette, 3),
      Gap_statistic = round(gap_statistic, 3),
      Bootstrap_ARI_mean = round(bootstrap_ARI_mean, 3),
      Bootstrap_ARI_SD = round(bootstrap_ARI_sd, 3),
      Cross_window_ARI_vs_full = round(cross_window_ARI_vs_full, 3),
      Minimum_cluster_size = min_cluster_size,
      Singleton_cluster_count = n_singleton_clusters,
      SGP_MYS_grouped_together = if_else(SGP_MYS_grouped_together, "Yes", "No"),
      Cluster_size_distribution = cluster_size_string
    )

  readr::write_csv(
    reviewer_tbl,
    file.path(comparison_dir, "Table_cluster_number_evaluation_for_reviewer.csv")
  )

  peak_compare <- purrr::map_dfr(all_results, "peak_tbl")
  cluster_peak_compare <- purrr::map_dfr(all_results, "cluster_peak_summary")

  readr::write_csv(peak_compare, file.path(comparison_dir, "country_peak_timing_comparison.csv"))
  readr::write_csv(cluster_peak_compare, file.path(comparison_dir, "cluster_peak_timing_comparison.csv"))
}

# -----------------------------------------------------------------------------
# Run analysis
# -----------------------------------------------------------------------------

dat <- read_weekly_data(DATA_PATH)

all_results <- purrr::imap(
  ANALYSIS_WINDOWS,
  ~ run_one_analysis(dat = dat, analysis_name = .y, keep_years = .x)
)

write_comparison_outputs(all_results)

message("DTW clustering analysis completed.")
message("Outputs saved to: ", normalizePath(OUTPUT_ROOT))
