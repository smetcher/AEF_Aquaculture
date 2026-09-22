# =============================================================================
# FINAL AEF 2024 DEVELOPMENT / REFERENCE PIPELINE
# Tasmania, Chile and Norway
# =============================================================================
#
# INPUTS (created by the final GEE reference-extraction workflow):
#
#   Analysis/Exports/AEF_Tasmania_2024_ReferencePixels.csv
#   Analysis/Exports/AEF_Chile_2024_ReferencePixels.csv
#   Analysis/Exports/AEF_Norway_2024_ReferencePixels.csv
#
# PURPOSE
#
#   1. QA the final 2024 reference data.
#   2. Build full-development class representatives.
#   3. Build K = 5 aquaculture prototypes.
#   4. Generate cosine similarities and class-separation margins.
#   5. Derive the 90/95/99% target-retention cosine thresholds.
#   6. Select production prototype variables using the pre-specified
#      >=80% cumulative prototype-importance rule.
#   7. Run leakage-free, stratified 5-fold SITE-level internal CV diagnostics.
#      All centroids, prototypes and prototype selection are rebuilt inside
#      each training fold before the held-out fold is predicted.
#   8. Run a SITE-level bootstrap for centroid and threshold stability.
#   9. Archive all vectors, diagnostics and provenance needed for the blind
#      2025 GEE evaluation and later manuscript/supplementary analyses.
#
# IMPORTANT
#
#   - This script uses ONLY 2024 development/reference data.
#   - The 2025 evaluation data must not influence any choice below.
#   - Production quantities use all 20 + 20 + 20 development sites.
#   - Internal CV is diagnostic only; the blind 2025 evaluation remains the
#     primary model assessment.
#   - Pixels are classifier observations, but resampling/splitting is always
#     performed by SITE.
#   - All reported similarity variables use explicit cosine normalisation of
#     both observation and reference vectors, matching the GEE workflow.
#
# REQUIRED PACKAGES
#
#   tidyverse
#   randomForest
#   cluster
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(cluster)
})


# =============================================================================
# 1. CONFIGURATION
# =============================================================================

YEAR <- 2024

# Run all regions in one pass with the following setting.
# To rerun only one region, e.g. Tasmania:
# REGIONS_TO_RUN <- c("Tasmania")
REGIONS_TO_RUN <- c("Tasmania", "Chile", "Norway")

INPUT_DIR <- "Analysis/Exports"
OUTPUT_DIR <- "Analysis/R_Outputs"

EXPECTED_SITES_PER_CLASS <- 20
EXPECTED_RADIUS_M <- 50
EXPECTED_SAMPLING_SCALE_M <- 10

CLASS_TARGET <- "target"
CLASS_OCEAN <- "ocean"
CLASS_CONF <- "confounding_marine"
EXPECTED_CLASSES <- c(CLASS_TARGET, CLASS_OCEAN, CLASS_CONF)

K <- 5
KMEANS_NSTART <- 50
KMEANS_ITER_MAX <- 100
CUMULATIVE_IMPORTANCE_TARGET <- 0.80

INTERNAL_CV_FOLDS <- 5
RF_TREES_INTERNAL <- 500
RF_TREES_PRODUCTION_DIAGNOSTIC <- 500

BOOTSTRAP_REPS <- 500

K_DIAGNOSTIC_RANGE <- 2:10
K_DIAGNOSTIC_MAX_PIXELS <- 2500

BASE_SEED <- 123
SEED_KMEANS <- BASE_SEED
SEED_PRODUCTION_SELECTOR <- BASE_SEED + 10
SEED_PRODUCTION_SRF <- BASE_SEED + 20
SEED_PRODUCTION_CRRF <- BASE_SEED + 30
SEED_CV_FOLDS <- BASE_SEED + 100
SEED_CV_BASE <- BASE_SEED + 1000
SEED_BOOTSTRAP <- BASE_SEED + 10000
SEED_K_DIAGNOSTIC <- BASE_SEED + 20000

NORM_WARNING_TOLERANCE <- 0.02

# Final GEE SRF/CR-RF settings are applied in GEE, not here.
# The R random forests below are development diagnostics / prototype selectors.


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 2. GENERAL HELPERS
# =============================================================================

normalise_vector <- function(v) {
  n <- sqrt(sum(v^2))
  
  if (!is.finite(n) || n <= 0) {
    stop("Cannot normalise a zero-length or non-finite vector.")
  }
  
  v / n
}


normalise_matrix_rows <- function(x, vector_label = "input") {
  row_norms <- sqrt(
    rowSums(x^2)
  )
  
  invalid_norms <-
    !is.finite(row_norms) |
    row_norms <= 0
  
  if (any(invalid_norms)) {
    stop(
      paste0(
        "Cannot calculate cosine similarity: one or more ",
        vector_label,
        " vectors have a non-finite or zero norm."
      )
    )
  }
  
  sweep(
    x,
    MARGIN = 1,
    STATS = row_norms,
    FUN = "/"
  )
}


safe_divide <- function(a, b) {
  ifelse(b == 0, NA_real_, a / b)
}


region_prefix <- function(region_name) {
  paste0("AEF_", region_name, "_", YEAR)
}


region_input_file <- function(region_name) {
  file.path(
    INPUT_DIR,
    paste0(region_prefix(region_name), "_ReferencePixels.csv")
  )
}


region_output_file <- function(region_name, suffix) {
  file.path(
    OUTPUT_DIR,
    paste0(region_prefix(region_name), suffix)
  )
}


write_region_output <- function(x, region_name, suffix) {
  path <- region_output_file(region_name, suffix)
  write.csv(x, path, row.names = FALSE)
  path
}


classification_metrics <- function(truth, prediction, classes) {
  truth <- factor(truth, levels = classes)
  prediction <- factor(prediction, levels = classes)
  
  cm <- table(
    Truth = truth,
    Predicted = prediction
  )
  
  per_class <- map_dfr(
    classes,
    function(cl) {
      tp <- as.numeric(cm[cl, cl])
      fn <- as.numeric(sum(cm[cl, ]) - tp)
      fp <- as.numeric(sum(cm[, cl]) - tp)
      tn <- as.numeric(sum(cm) - tp - fn - fp)
      
      precision <- safe_divide(tp, tp + fp)
      recall <- safe_divide(tp, tp + fn)
      specificity <- safe_divide(tn, tn + fp)
      iou <- safe_divide(tp, tp + fp + fn)
      
      f1 <- ifelse(
        is.na(precision) ||
          is.na(recall) ||
          (precision + recall == 0),
        NA_real_,
        2 * precision * recall / (precision + recall)
      )
      
      tibble(
        Class = cl,
        TP = tp,
        TN = tn,
        FP = fp,
        FN = fn,
        Precision = precision,
        Recall = recall,
        Specificity = specificity,
        F1 = f1,
        IoU = iou
      )
    }
  )
  
  overall <- tibble(
    Accuracy = mean(truth == prediction),
    Balanced_Accuracy = mean(per_class$Recall, na.rm = TRUE),
    Macro_F1 = mean(per_class$F1, na.rm = TRUE),
    Macro_IoU = mean(per_class$IoU, na.rm = TRUE)
  )
  
  confusion <- as.data.frame(cm) %>%
    rename(N = Freq)
  
  list(
    overall = overall,
    per_class = per_class,
    confusion = confusion
  )
}


target_vs_rest_metrics <- function(truth, prediction, target_class) {
  truth_positive <- truth == target_class
  pred_positive <- prediction == target_class
  
  tp <- sum(truth_positive & pred_positive)
  tn <- sum(!truth_positive & !pred_positive)
  fp <- sum(!truth_positive & pred_positive)
  fn <- sum(truth_positive & !pred_positive)
  
  precision <- safe_divide(tp, tp + fp)
  recall <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  f1 <- ifelse(
    is.na(precision) ||
      is.na(recall) ||
      (precision + recall == 0),
    NA_real_,
    2 * precision * recall / (precision + recall)
  )
  iou <- safe_divide(tp, tp + fp + fn)
  
  tibble(
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn,
    Precision = precision,
    Recall = recall,
    Specificity = specificity,
    F1 = f1,
    IoU = iou
  )
}


site_majority_predictions <- function(site_key, truth, prediction) {
  tibble(
    site_key = site_key,
    Truth = as.character(truth),
    Predicted = as.character(prediction)
  ) %>%
    count(
      site_key,
      Truth,
      Predicted,
      name = "Votes"
    ) %>%
    arrange(
      site_key,
      desc(Votes),
      Predicted
    ) %>%
    group_by(
      site_key,
      Truth
    ) %>%
    slice(1) %>%
    ungroup()
}


fit_rf <- function(
    train_df,
    test_df,
    feature_names,
    classes,
    ntree,
    seed
) {
  set.seed(seed)
  
  model <- randomForest(
    x = as.data.frame(
      train_df[, feature_names, drop = FALSE]
    ),
    y = factor(
      train_df$class,
      levels = classes
    ),
    ntree = ntree,
    importance = TRUE
  )
  
  prediction <- predict(
    model,
    as.data.frame(
      test_df[, feature_names, drop = FALSE]
    )
  )
  
  list(
    model = model,
    prediction = as.character(prediction)
  )
}


density_overlap <- function(x, y, n = 4096) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }
  
  lower <- min(c(x, y))
  upper <- max(c(x, y))
  
  if (lower == upper) {
    return(1)
  }
  
  range_width <- upper - lower
  fallback_bw <- max(range_width / 1000, .Machine$double.eps)
  
  bw_x <- bw.nrd0(x)
  bw_y <- bw.nrd0(y)
  
  if (!is.finite(bw_x) || bw_x <= 0) {
    bw_x <- fallback_bw
  }
  
  if (!is.finite(bw_y) || bw_y <= 0) {
    bw_y <- fallback_bw
  }
  
  dx <- density(
    x,
    from = lower,
    to = upper,
    n = n,
    bw = bw_x
  )
  
  dy <- density(
    y,
    from = lower,
    to = upper,
    n = n,
    bw = bw_y
  )
  
  minimum_density <- pmin(dx$y, dy$y)
  
  overlap <- sum(
    diff(dx$x) *
      (
        head(minimum_density, -1) +
          tail(minimum_density, -1)
      ) / 2
  )
  
  max(0, min(1, overlap))
}


extract_gini_importance <- function(model) {
  imp <- importance(model)
  
  if (!"MeanDecreaseGini" %in% colnames(imp)) {
    stop("MeanDecreaseGini was not returned by randomForest::importance().")
  }
  
  tibble(
    Feature = rownames(imp),
    MeanDecreaseGini = imp[, "MeanDecreaseGini"]
  ) %>%
    arrange(desc(MeanDecreaseGini))
}


select_prototypes_from_importance <- function(
    importance_df,
    prototype_names,
    cumulative_target
) {
  proto_importance <- importance_df %>%
    filter(Feature %in% prototype_names) %>%
    arrange(desc(MeanDecreaseGini))
  
  importance_sum <- sum(
    proto_importance$MeanDecreaseGini,
    na.rm = TRUE
  )
  
  if (!is.finite(importance_sum) || importance_sum <= 0) {
    stop(
      "Prototype importance sums to zero or less; prototype selection cannot proceed."
    )
  }
  
  proto_importance <- proto_importance %>%
    mutate(
      Relative_Importance = MeanDecreaseGini / importance_sum,
      Cumulative_Importance = cumsum(Relative_Importance)
    )
  
  cutoff <- which(
    proto_importance$Cumulative_Importance >= cumulative_target
  )[1]
  
  selected <- proto_importance[
    seq_len(cutoff),
    ,
    drop = FALSE
  ]
  
  list(
    all = proto_importance,
    selected = selected,
    names = selected$Feature
  )
}


sample_site_indices <- function(index_list, sampled_sites) {
  unlist(
    lapply(
      sampled_sites,
      function(site_name) {
        index_list[[site_name]]
      }
    ),
    use.names = FALSE
  )
}


# =============================================================================
# 3. ENGINEERING HELPERS
# =============================================================================

build_engineering_model <- function(
    training_df,
    embedding_cols,
    k,
    kmeans_nstart,
    kmeans_iter_max,
    seed
) {
  target_mat <- as.matrix(
    training_df[
      training_df$class == CLASS_TARGET,
      embedding_cols,
      drop = FALSE
    ]
  )
  
  ocean_mat <- as.matrix(
    training_df[
      training_df$class == CLASS_OCEAN,
      embedding_cols,
      drop = FALSE
    ]
  )
  
  conf_mat <- as.matrix(
    training_df[
      training_df$class == CLASS_CONF,
      embedding_cols,
      drop = FALSE
    ]
  )
  
  target_centroid_raw <- colMeans(target_mat)
  ocean_centroid_raw <- colMeans(ocean_mat)
  conf_centroid_raw <- colMeans(conf_mat)
  
  target_centroid <- normalise_vector(target_centroid_raw)
  ocean_centroid <- normalise_vector(ocean_centroid_raw)
  conf_centroid <- normalise_vector(conf_centroid_raw)
  
  set.seed(seed)
  
  km <- kmeans(
    target_mat,
    centers = k,
    nstart = kmeans_nstart,
    iter.max = kmeans_iter_max,
    algorithm = "Hartigan-Wong"
  )
  
  prototypes_raw <- km$centers
  prototypes <- prototypes_raw
  
  for (i in seq_len(nrow(prototypes))) {
    prototypes[i, ] <- normalise_vector(prototypes[i, ])
  }
  
  prototype_names <- paste0("proto_", seq_len(k))
  
  rownames(prototypes_raw) <- prototype_names
  rownames(prototypes) <- prototype_names
  
  list(
    target_centroid_raw = target_centroid_raw,
    ocean_centroid_raw = ocean_centroid_raw,
    conf_centroid_raw = conf_centroid_raw,
    target_centroid = target_centroid,
    ocean_centroid = ocean_centroid,
    conf_centroid = conf_centroid,
    kmeans = km,
    prototypes_raw = prototypes_raw,
    prototypes = prototypes,
    prototype_names = prototype_names
  )
}


apply_engineering_model <- function(
    data,
    engineering_model,
    embedding_cols
) {
  out <- data
  
  mat <- as.matrix(
    out[, embedding_cols, drop = FALSE]
  )
  
  # Calculate explicit cosine similarity rather than assuming that every
  # exported AEF pixel vector has an exact L2 norm of one. The class and
  # prototype reference vectors are already normalised in
  # build_engineering_model(); normalising each pixel vector here makes the
  # R-derived similarities mathematically identical to those calculated in
  # the GEE classification workflow.
  mat_unit <- normalise_matrix_rows(
    mat,
    vector_label = "AEF pixel"
  )
  
  out$farm_sim <- as.numeric(
    mat_unit %*% engineering_model$target_centroid
  )
  
  out$ocean_sim <- as.numeric(
    mat_unit %*% engineering_model$ocean_centroid
  )
  
  out$conf_sim <- as.numeric(
    mat_unit %*% engineering_model$conf_centroid
  )
  
  for (i in seq_along(engineering_model$prototype_names)) {
    proto_name <- engineering_model$prototype_names[i]
    
    out[[proto_name]] <- as.numeric(
      mat_unit %*% engineering_model$prototypes[i, ]
    )
  }
  
  out$margin_ocean <- out$farm_sim - out$ocean_sim
  out$margin_conf <- out$farm_sim - out$conf_sim
  out$margin_best_noise <- out$farm_sim - pmax(
    out$ocean_sim,
    out$conf_sim
  )
  
  out
}


# =============================================================================
# 4. REGION PROCESSOR
# =============================================================================

process_region <- function(region_name) {
  cat("\n\n============================================================\n")
  cat("REGION:", region_name, "\n")
  cat("============================================================\n")
  
  input_file <- region_input_file(region_name)
  
  if (!file.exists(input_file)) {
    stop(
      paste0(
        "Input file not found: ",
        input_file
      )
    )
  }
  
  created_files <- character(0)
  
  
  # ==========================================================================
  # 4.1 LOAD DATA
  # ==========================================================================
  
  cat("\nLoading:", input_file, "\n")
  
  df <- read.csv(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  required_columns <- c(
    "class",
    "region",
    "year",
    "sampling_radius_m",
    "sampling_scale_m"
  )
  
  missing_required_columns <- setdiff(
    required_columns,
    names(df)
  )
  
  if (length(missing_required_columns) > 0) {
    stop(
      paste0(
        "Missing required columns: ",
        paste(missing_required_columns, collapse = ", ")
      )
    )
  }
  
  site_id_col <- if (
    "site_id" %in% names(df)
  ) {
    "site_id"
  } else if (
    "poly_id" %in% names(df)
  ) {
    "poly_id"
  } else {
    stop("Neither 'site_id' nor 'poly_id' was found in the input.")
  }
  
  df$class <- trimws(as.character(df$class))
  df[[site_id_col]] <- as.character(df[[site_id_col]])
  
  
  # ==========================================================================
  # 4.2 IDENTIFY AND CHECK AEF DIMENSIONS
  # ==========================================================================
  
  embedding_cols <- grep(
    "^A[0-9]{2}$",
    names(df),
    value = TRUE
  )
  
  embedding_index <- as.integer(
    sub("^A", "", embedding_cols)
  )
  
  embedding_cols <- embedding_cols[
    order(embedding_index)
  ]
  
  embedding_index <- sort(embedding_index)
  
  if (length(embedding_cols) != 64) {
    stop(
      paste0(
        "Expected 64 AEF dimensions but found ",
        length(embedding_cols),
        "."
      )
    )
  }
  
  if (!identical(embedding_index, 0:63)) {
    stop(
      "AEF columns do not correspond exactly to A00-A63."
    )
  }
  
  if (!all(map_lgl(df[, embedding_cols, drop = FALSE], is.numeric))) {
    stop(
      "One or more AEF embedding columns are not numeric."
    )
  }
  
  embedding_matrix_initial <- as.matrix(
    df[, embedding_cols, drop = FALSE]
  )
  
  if (any(!is.finite(embedding_matrix_initial))) {
    stop(
      "Missing or non-finite values detected in AEF dimensions."
    )
  }
  
  cat(
    "Samples:", nrow(df),
    "\nEmbedding dimensions:", length(embedding_cols),
    "\n"
  )
  
  
  # ==========================================================================
  # 4.3 BASIC PROVENANCE QA
  # ==========================================================================
  
  unexpected_classes <- setdiff(
    unique(df$class),
    EXPECTED_CLASSES
  )
  
  missing_classes <- setdiff(
    EXPECTED_CLASSES,
    unique(df$class)
  )
  
  if (length(unexpected_classes) > 0) {
    stop(
      paste0(
        "Unexpected classes found: ",
        paste(unexpected_classes, collapse = ", ")
      )
    )
  }
  
  if (length(missing_classes) > 0) {
    stop(
      paste0(
        "Missing expected classes: ",
        paste(missing_classes, collapse = ", ")
      )
    )
  }
  
  if (any(as.integer(df$year) != YEAR)) {
    stop(
      paste0(
        "Input contains observations outside YEAR = ",
        YEAR,
        "."
      )
    )
  }
  
  region_values <- unique(trimws(as.character(df$region)))
  
  if (
    length(region_values) != 1 ||
    region_values[1] != region_name
  ) {
    stop(
      paste0(
        "Input region field does not match expected region '",
        region_name,
        "'. Found: ",
        paste(region_values, collapse = ", ")
      )
    )
  }
  
  radius_values <- unique(
    as.numeric(df$sampling_radius_m)
  )
  
  radius_values <- radius_values[
    is.finite(radius_values)
  ]
  
  if (
    length(radius_values) != 1 ||
    abs(radius_values[1] - EXPECTED_RADIUS_M) > 1e-8
  ) {
    stop(
      paste0(
        "Sampling radius is not uniformly ",
        EXPECTED_RADIUS_M,
        " m."
      )
    )
  }
  
  scale_values <- unique(
    as.numeric(df$sampling_scale_m)
  )
  
  scale_values <- scale_values[
    is.finite(scale_values)
  ]
  
  if (
    length(scale_values) != 1 ||
    abs(scale_values[1] - EXPECTED_SAMPLING_SCALE_M) > 1e-8
  ) {
    stop(
      paste0(
        "Sampling scale is not uniformly ",
        EXPECTED_SAMPLING_SCALE_M,
        " m."
      )
    )
  }
  
  if (
    "site_id" %in% names(df) &&
    "poly_id" %in% names(df)
  ) {
    if (
      any(
        as.character(df$site_id) !=
        as.character(df$poly_id)
      )
    ) {
      stop(
        "site_id and poly_id are not identical in the final export."
      )
    }
  }
  
  
  # ==========================================================================
  # 4.4 DETERMINISTIC ROW ORDER AND SITE KEY
  # ==========================================================================
  
  df$site_key <- paste(
    df$class,
    df[[site_id_col]],
    sep = "::"
  )
  
  ordering_data <- df %>%
    select(
      class,
      site_key,
      any_of(c("pixel_x", "pixel_y", "pixel_lon", "pixel_lat")),
      all_of(embedding_cols)
    )
  
  row_order <- do.call(
    order,
    ordering_data
  )
  
  df <- df[
    row_order,
    ,
    drop = FALSE
  ]
  
  rownames(df) <- NULL
  df$row_id <- seq_len(nrow(df))
  
  embeddings <- as.matrix(
    df[, embedding_cols, drop = FALSE]
  )
  
  
  # ==========================================================================
  # 4.5 SITE QA
  # ==========================================================================
  
  site_qa <- df %>%
    count(
      class,
      site_key,
      .data[[site_id_col]],
      name = "n_pixels"
    )
  
  names(site_qa)[
    names(site_qa) == site_id_col
  ] <- "site_id"
  
  metadata_candidates <- c(
    "source_id",
    "centroid_lon",
    "centroid_lat",
    "sampling_radius_m",
    "sampling_scale_m",
    "sampling_crs",
    "sampling_utm_zone",
    "source_aef_utm_zones",
    "aef_dataset_versions",
    "aef_model_versions",
    "aef_processing_versions"
  )
  
  metadata_cols <- intersect(
    metadata_candidates,
    names(df)
  )
  
  if (length(metadata_cols) > 0) {
    site_metadata <- df %>%
      group_by(
        class,
        site_key
      ) %>%
      summarise(
        across(
          all_of(metadata_cols),
          ~ first(.x)
        ),
        .groups = "drop"
      )
    
    site_qa <- site_qa %>%
      left_join(
        site_metadata,
        by = c("class", "site_key")
      )
  }
  
  site_counts <- site_qa %>%
    count(
      class,
      name = "n_sites"
    )
  
  cat("\nUnique sites per class:\n")
  print(site_counts)
  
  if (
    nrow(site_counts) != length(EXPECTED_CLASSES) ||
    any(site_counts$n_sites != EXPECTED_SITES_PER_CLASS)
  ) {
    stop(
      paste0(
        "Expected exactly ",
        EXPECTED_SITES_PER_CLASS,
        " sites in each class."
      )
    )
  }
  
  
  # ==========================================================================
  # 4.6 PIXEL SUPPORT QA
  # ==========================================================================
  
  pixel_support_summary <- site_qa %>%
    group_by(class) %>%
    summarise(
      n_sites = n(),
      total_pixels = sum(n_pixels),
      min_pixels_site = min(n_pixels),
      q25_pixels_site = as.numeric(
        quantile(n_pixels, 0.25)
      ),
      median_pixels_site = median(n_pixels),
      mean_pixels_site = mean(n_pixels),
      q75_pixels_site = as.numeric(
        quantile(n_pixels, 0.75)
      ),
      max_pixels_site = max(n_pixels),
      sd_pixels_site = sd(n_pixels),
      cv_pixels_site = safe_divide(
        sd(n_pixels),
        mean(n_pixels)
      ),
      .groups = "drop"
    )
  
  cat("\nPixel support by class:\n")
  print(pixel_support_summary)
  
  
  # ==========================================================================
  # 4.7 EMBEDDING NORM QA
  # ==========================================================================
  
  embedding_norms <- sqrt(
    rowSums(embeddings^2)
  )
  
  embedding_norm_summary <- tibble(
    Metric = c(
      "Minimum",
      "Q025",
      "Median",
      "Mean",
      "Q975",
      "Maximum",
      "Mean_Absolute_Deviation_From_1",
      "Maximum_Absolute_Deviation_From_1"
    ),
    Value = c(
      min(embedding_norms),
      as.numeric(quantile(embedding_norms, 0.025)),
      median(embedding_norms),
      mean(embedding_norms),
      as.numeric(quantile(embedding_norms, 0.975)),
      max(embedding_norms),
      mean(abs(embedding_norms - 1)),
      max(abs(embedding_norms - 1))
    )
  )
  
  if (
    max(abs(embedding_norms - 1)) > NORM_WARNING_TOLERANCE
  ) {
    warning(
      paste0(
        region_name,
        ": some AEF vectors deviate from unit length by more than ",
        NORM_WARNING_TOLERANCE,
        ". Inspect the embedding norm QA output."
      )
    )
  }
  
  
  # ==========================================================================
  # 4.8 SITE-MEAN EMBEDDINGS
  # ==========================================================================
  
  site_mean_embeddings <- df %>%
    group_by(
      class,
      site_key,
      .data[[site_id_col]]
    ) %>%
    summarise(
      across(
        all_of(embedding_cols),
        mean
      ),
      n_pixels = n(),
      .groups = "drop"
    )
  
  names(site_mean_embeddings)[
    names(site_mean_embeddings) == site_id_col
  ] <- "site_id"
  
  
  # ==========================================================================
  # 4.9 FULL-DEVELOPMENT ENGINEERING MODEL
  # ==========================================================================
  
  cat("\nBuilding full-development class representatives and K =", K, "prototypes...\n")
  
  production_engineering <- build_engineering_model(
    training_df = df,
    embedding_cols = embedding_cols,
    k = K,
    kmeans_nstart = KMEANS_NSTART,
    kmeans_iter_max = KMEANS_ITER_MAX,
    seed = SEED_KMEANS
  )
  
  production_df <- apply_engineering_model(
    data = df,
    engineering_model = production_engineering,
    embedding_cols = embedding_cols
  )
  
  prototype_names <- production_engineering$prototype_names
  
  target_rows <- which(
    production_df$class == CLASS_TARGET
  )
  
  production_df$prototype_cluster <- NA_integer_
  production_df$prototype_cluster[target_rows] <-
    production_engineering$kmeans$cluster
  
  
  # ==========================================================================
  # 4.10 CLASS REPRESENTATIVE ARCHIVE
  # ==========================================================================
  
  centroid_archive <- tibble(
    band = embedding_cols,
    farm_raw = production_engineering$target_centroid_raw,
    farm_normalised = production_engineering$target_centroid,
    ocean_raw = production_engineering$ocean_centroid_raw,
    ocean_normalised = production_engineering$ocean_centroid,
    conf_raw = production_engineering$conf_centroid_raw,
    conf_normalised = production_engineering$conf_centroid
  )
  
  centroid_relationship_matrix <- rbind(
    farm = production_engineering$target_centroid,
    ocean = production_engineering$ocean_centroid,
    conf = production_engineering$conf_centroid
  )
  
  centroid_relationship_matrix <-
    centroid_relationship_matrix %*%
    t(centroid_relationship_matrix)
  
  centroid_relationships <- as.data.frame(
    as.table(centroid_relationship_matrix)
  ) %>%
    rename(
      Class_1 = Var1,
      Class_2 = Var2,
      Cosine_Similarity = Freq
    )
  
  
  # ==========================================================================
  # 4.11 COSINE THRESHOLDS
  # ==========================================================================
  
  target_similarity_values <- production_df$farm_sim[
    production_df$class == CLASS_TARGET
  ]
  
  cosine_thresholds <- tibble(
    Retention = c(90, 95, 99),
    Quantile = c(0.10, 0.05, 0.01),
    Threshold = as.numeric(
      quantile(
        target_similarity_values,
        probs = c(0.10, 0.05, 0.01),
        names = FALSE
      )
    )
  )
  
  cat("\nCosine thresholds:\n")
  print(cosine_thresholds)
  
  
  # ==========================================================================
  # 4.12 SIMILARITY / MARGIN SUMMARIES
  # ==========================================================================
  
  class_similarity_cols <- c(
    "farm_sim",
    "ocean_sim",
    "conf_sim"
  )
  
  margin_cols <- c(
    "margin_ocean",
    "margin_conf",
    "margin_best_noise"
  )
  
  diagnostic_variables <- c(
    class_similarity_cols,
    margin_cols,
    prototype_names
  )
  
  similarity_summary <- production_df %>%
    select(
      class,
      all_of(diagnostic_variables)
    ) %>%
    pivot_longer(
      cols = -class,
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    group_by(
      class,
      Variable
    ) %>%
    summarise(
      N = n(),
      Mean = mean(Value),
      SD = sd(Value),
      Q01 = as.numeric(quantile(Value, 0.01)),
      Q05 = as.numeric(quantile(Value, 0.05)),
      Q10 = as.numeric(quantile(Value, 0.10)),
      Q25 = as.numeric(quantile(Value, 0.25)),
      Median = median(Value),
      Q75 = as.numeric(quantile(Value, 0.75)),
      Q90 = as.numeric(quantile(Value, 0.90)),
      Q95 = as.numeric(quantile(Value, 0.95)),
      Q99 = as.numeric(quantile(Value, 0.99)),
      .groups = "drop"
    )
  
  class_pairs <- list(
    c(CLASS_TARGET, CLASS_OCEAN),
    c(CLASS_TARGET, CLASS_CONF),
    c(CLASS_OCEAN, CLASS_CONF)
  )
  
  overlap_diagnostics <- map_dfr(
    diagnostic_variables,
    function(variable_name) {
      map_dfr(
        class_pairs,
        function(pair) {
          x <- production_df[[variable_name]][
            production_df$class == pair[1]
          ]
          
          y <- production_df[[variable_name]][
            production_df$class == pair[2]
          ]
          
          tibble(
            Variable = variable_name,
            Class_1 = pair[1],
            Class_2 = pair[2],
            Kernel_Density_Overlap = density_overlap(x, y)
          )
        }
      )
    }
  )
  
  
  # ==========================================================================
  # 4.13 SITE-MEAN RELATIONSHIPS
  # ==========================================================================
  
  site_mean_matrix <- as.matrix(
    site_mean_embeddings[, embedding_cols, drop = FALSE]
  )
  
  site_mean_matrix_unit <- normalise_matrix_rows(
    site_mean_matrix,
    vector_label = "site-mean AEF"
  )
  
  site_mean_relationships <- site_mean_embeddings %>%
    select(
      class,
      site_key,
      site_id,
      n_pixels
    ) %>%
    mutate(
      farm_sim = as.numeric(
        site_mean_matrix_unit %*%
          production_engineering$target_centroid
      ),
      ocean_sim = as.numeric(
        site_mean_matrix_unit %*%
          production_engineering$ocean_centroid
      ),
      conf_sim = as.numeric(
        site_mean_matrix_unit %*%
          production_engineering$conf_centroid
      ),
      margin_ocean = farm_sim - ocean_sim,
      margin_conf = farm_sim - conf_sim,
      margin_best_noise = farm_sim - pmax(
        ocean_sim,
        conf_sim
      )
    )
  
  
  # ==========================================================================
  # 4.14 PROTOTYPE DIAGNOSTICS
  # ==========================================================================
  
  prototypes <- production_engineering$prototypes
  prototypes_raw <- production_engineering$prototypes_raw
  
  prototype_similarity_matrix <-
    prototypes %*%
    t(prototypes)
  
  dimnames(prototype_similarity_matrix) <- list(
    prototype_names,
    prototype_names
  )
  
  prototype_similarity_long <- as.data.frame(
    as.table(prototype_similarity_matrix)
  ) %>%
    rename(
      Prototype_1 = Var1,
      Prototype_2 = Var2,
      Cosine_Similarity = Freq
    )
  
  prototype_site_membership <- production_df %>%
    filter(
      class == CLASS_TARGET
    ) %>%
    count(
      prototype_cluster,
      site_key,
      .data[[site_id_col]],
      name = "n_pixels"
    )
  
  names(prototype_site_membership)[
    names(prototype_site_membership) == site_id_col
  ] <- "site_id"
  
  prototype_cluster_summary <- prototype_site_membership %>%
    group_by(
      prototype_cluster
    ) %>%
    summarise(
      n_pixels = sum(n_pixels),
      n_sites = n_distinct(site_key),
      .groups = "drop"
    )
  
  prototype_archive <- tibble(
    band = embedding_cols
  )
  
  for (i in seq_len(K)) {
    prototype_archive[[paste0("proto_", i, "_raw")]] <-
      prototypes_raw[i, ]
    
    prototype_archive[[paste0("proto_", i, "_normalised")]] <-
      prototypes[i, ]
  }
  
  
  # ==========================================================================
  # 4.15 K-MEANS DIAGNOSTICS
  # ==========================================================================
  
  cat("\nCalculating K-means diagnostics...\n")
  
  target_matrix <- as.matrix(
    production_df[
      production_df$class == CLASS_TARGET,
      embedding_cols,
      drop = FALSE
    ]
  )
  
  set.seed(SEED_K_DIAGNOSTIC)
  
  if (
    nrow(target_matrix) > K_DIAGNOSTIC_MAX_PIXELS
  ) {
    k_diag_rows <- sort(
      sample(
        seq_len(nrow(target_matrix)),
        K_DIAGNOSTIC_MAX_PIXELS,
        replace = FALSE
      )
    )
    
    target_matrix_kdiag <- target_matrix[
      k_diag_rows,
      ,
      drop = FALSE
    ]
  } else {
    target_matrix_kdiag <- target_matrix
  }
  
  target_distance <- dist(
    target_matrix_kdiag
  )
  
  k_diagnostic_list <- vector(
    "list",
    length(K_DIAGNOSTIC_RANGE)
  )
  
  for (j in seq_along(K_DIAGNOSTIC_RANGE)) {
    k_value <- K_DIAGNOSTIC_RANGE[j]
    
    set.seed(
      SEED_K_DIAGNOSTIC + k_value
    )
    
    k_result <- tryCatch(
      {
        km_temp <- kmeans(
          target_matrix_kdiag,
          centers = k_value,
          nstart = KMEANS_NSTART,
          iter.max = KMEANS_ITER_MAX,
          algorithm = "Hartigan-Wong"
        )
        
        sil <- silhouette(
          km_temp$cluster,
          target_distance
        )
        
        tibble(
          K = k_value,
          N_Pixels_Used = nrow(target_matrix_kdiag),
          Total_WithinSS = km_temp$tot.withinss,
          BetweenSS = km_temp$betweenss,
          TotalSS = km_temp$totss,
          Variance_Explained = safe_divide(
            km_temp$betweenss,
            km_temp$totss
          ),
          Mean_Silhouette = mean(
            sil[, "sil_width"]
          )
        )
      },
      error = function(e) {
        warning(
          paste0(
            region_name,
            ": K diagnostic failed for K = ",
            k_value,
            ": ",
            conditionMessage(e)
          )
        )
        
        tibble(
          K = k_value,
          N_Pixels_Used = nrow(target_matrix_kdiag),
          Total_WithinSS = NA_real_,
          BetweenSS = NA_real_,
          TotalSS = NA_real_,
          Variance_Explained = NA_real_,
          Mean_Silhouette = NA_real_
        )
      }
    )
    
    k_diagnostic_list[[j]] <- k_result
  }
  
  k_diagnostics <- bind_rows(
    k_diagnostic_list
  )
  
  
  # ==========================================================================
  # 4.16 PRODUCTION PROTOTYPE SELECTION
  # ==========================================================================
  # All 2024 reference sites are legitimate development data.
  # The genuinely blind assessment is the spatially + temporally held-out
  # 2025 evaluation AOI. Therefore production prototype selection is fitted
  # using the complete 2024 development set.
  # ==========================================================================
  
  selector_features <- c(
    embedding_cols,
    class_similarity_cols,
    margin_cols,
    prototype_names
  )
  
  set.seed(SEED_PRODUCTION_SELECTOR)
  
  production_selector_model <- randomForest(
    x = as.data.frame(
      production_df[, selector_features, drop = FALSE]
    ),
    y = factor(
      production_df$class,
      levels = EXPECTED_CLASSES
    ),
    ntree = RF_TREES_PRODUCTION_DIAGNOSTIC,
    importance = TRUE
  )
  
  production_selector_importance <- extract_gini_importance(
    production_selector_model
  )
  
  production_prototype_selection <-
    select_prototypes_from_importance(
      importance_df = production_selector_importance,
      prototype_names = prototype_names,
      cumulative_target = CUMULATIVE_IMPORTANCE_TARGET
    )
  
  prototype_importance <-
    production_prototype_selection$all
  
  selected_prototypes <-
    production_prototype_selection$selected
  
  selected_proto_names <-
    production_prototype_selection$names
  
  cat("\nProduction selected prototypes:\n")
  print(selected_proto_names)
  
  
  # ==========================================================================
  # 4.17 PRODUCTION FEATURE SETS
  # ==========================================================================
  
  production_feature_sets <- list(
    Raw64 = embedding_cols,
    Raw64_ClassSimilarity = c(
      embedding_cols,
      class_similarity_cols
    ),
    Raw64_ClassSimilarity_Margins = c(
      embedding_cols,
      class_similarity_cols,
      margin_cols
    ),
    Raw64_ClassSimilarity_Margins_AllPrototypes = c(
      embedding_cols,
      class_similarity_cols,
      margin_cols,
      prototype_names
    ),
    CRRF_Selected = c(
      embedding_cols,
      class_similarity_cols,
      margin_cols,
      selected_proto_names
    ),
    EngineeredOnly_Selected = c(
      class_similarity_cols,
      margin_cols,
      selected_proto_names
    )
  )
  
  production_feature_manifest <- imap_dfr(
    production_feature_sets,
    function(features, variant_name) {
      tibble(
        Variant = variant_name,
        Feature_Order = seq_along(features),
        Feature = features
      )
    }
  )
  
  production_feature_sets_compact <- imap_dfr(
    production_feature_sets,
    function(features, variant_name) {
      tibble(
        Variant = variant_name,
        N_Features = length(features),
        Features = paste(
          features,
          collapse = ";"
        )
      )
    }
  )
  
  
  # ==========================================================================
  # 4.18 PRODUCTION RF IMPORTANCE ARCHIVE
  # ==========================================================================
  # These models are descriptive only. Final mapping RFs are trained in GEE.
  # ==========================================================================
  
  set.seed(SEED_PRODUCTION_SRF)
  
  production_srf_model <- randomForest(
    x = as.data.frame(
      production_df[, embedding_cols, drop = FALSE]
    ),
    y = factor(
      production_df$class,
      levels = EXPECTED_CLASSES
    ),
    ntree = RF_TREES_PRODUCTION_DIAGNOSTIC,
    importance = TRUE
  )
  
  production_srf_importance <- extract_gini_importance(
    production_srf_model
  )
  
  set.seed(SEED_PRODUCTION_CRRF)
  
  production_crrf_model <- randomForest(
    x = as.data.frame(
      production_df[
        ,
        production_feature_sets$CRRF_Selected,
        drop = FALSE
      ]
    ),
    y = factor(
      production_df$class,
      levels = EXPECTED_CLASSES
    ),
    ntree = RF_TREES_PRODUCTION_DIAGNOSTIC,
    importance = TRUE
  )
  
  production_crrf_importance <- extract_gini_importance(
    production_crrf_model
  )
  
  
  # ==========================================================================
  # 4.19 STRATIFIED 5-FOLD SITE ASSIGNMENT
  # ==========================================================================
  # With 20 sites/class and 5 folds, each held-out fold contains exactly
  # 4 independent sites from each class.
  # ==========================================================================
  
  site_table <- production_df %>%
    distinct(
      class,
      site_key,
      .data[[site_id_col]]
    )
  
  names(site_table)[
    names(site_table) == site_id_col
  ] <- "site_id"
  
  if (
    EXPECTED_SITES_PER_CLASS %% INTERNAL_CV_FOLDS != 0
  ) {
    stop(
      "EXPECTED_SITES_PER_CLASS must be divisible by INTERNAL_CV_FOLDS."
    )
  }
  
  fold_assignments <- map_dfr(
    seq_along(EXPECTED_CLASSES),
    function(class_index) {
      cl <- EXPECTED_CLASSES[class_index]
      
      class_sites <- site_table %>%
        filter(class == cl) %>%
        arrange(site_key)
      
      set.seed(
        SEED_CV_FOLDS + class_index
      )
      
      shuffled_indices <- sample(
        seq_len(nrow(class_sites)),
        size = nrow(class_sites),
        replace = FALSE
      )
      
      shuffled_sites <- class_sites[
        shuffled_indices,
        ,
        drop = FALSE
      ]
      
      shuffled_sites$Fold <- rep(
        seq_len(INTERNAL_CV_FOLDS),
        length.out = nrow(shuffled_sites)
      )
      
      shuffled_sites
    }
  )
  
  fold_summary <- fold_assignments %>%
    count(
      class,
      Fold,
      name = "n_sites"
    )
  
  cat("\nInternal CV site allocation:\n")
  print(fold_summary)
  
  
  # ==========================================================================
  # 4.20 LEAKAGE-FREE INTERNAL CV ABLATION
  # ==========================================================================
  # For each fold:
  #   - held-out sites are excluded first;
  #   - class representatives are rebuilt from training sites only;
  #   - K=5 prototypes are rebuilt from training target sites only;
  #   - engineered variables are generated from training-only quantities;
  #   - prototype selection is fitted on training-fold data only;
  #   - the untouched held-out sites are then predicted.
  # ==========================================================================
  
  cv_prediction_list <- list()
  cv_fold_summary_list <- list()
  cv_fold_pixel_class_list <- list()
  cv_fold_site_class_list <- list()
  cv_fold_confusion_list <- list()
  cv_fold_site_votes_list <- list()
  cv_selected_proto_list <- list()
  cv_selector_importance_list <- list()
  
  cv_record_index <- 1
  
  for (fold in seq_len(INTERNAL_CV_FOLDS)) {
    cat(
      "\nInternal CV fold",
      fold,
      "/",
      INTERNAL_CV_FOLDS,
      "\n"
    )
    
    held_out_keys <- fold_assignments %>%
      filter(Fold == fold) %>%
      pull(site_key)
    
    train_raw <- production_df %>%
      filter(!site_key %in% held_out_keys) %>%
      select(-any_of(diagnostic_variables)) %>%
      select(-any_of("prototype_cluster"))
    
    test_raw <- production_df %>%
      filter(site_key %in% held_out_keys) %>%
      select(-any_of(diagnostic_variables)) %>%
      select(-any_of("prototype_cluster"))
    
    fold_engineering <- build_engineering_model(
      training_df = train_raw,
      embedding_cols = embedding_cols,
      k = K,
      kmeans_nstart = KMEANS_NSTART,
      kmeans_iter_max = KMEANS_ITER_MAX,
      seed = SEED_CV_BASE + fold
    )
    
    train_fold <- apply_engineering_model(
      data = train_raw,
      engineering_model = fold_engineering,
      embedding_cols = embedding_cols
    )
    
    test_fold <- apply_engineering_model(
      data = test_raw,
      engineering_model = fold_engineering,
      embedding_cols = embedding_cols
    )
    
    fold_proto_names <- fold_engineering$prototype_names
    
    fold_selector_features <- c(
      embedding_cols,
      class_similarity_cols,
      margin_cols,
      fold_proto_names
    )
    
    set.seed(
      SEED_CV_BASE + 100 + fold
    )
    
    fold_selector_model <- randomForest(
      x = as.data.frame(
        train_fold[
          ,
          fold_selector_features,
          drop = FALSE
        ]
      ),
      y = factor(
        train_fold$class,
        levels = EXPECTED_CLASSES
      ),
      ntree = RF_TREES_INTERNAL,
      importance = TRUE
    )
    
    fold_selector_importance <- extract_gini_importance(
      fold_selector_model
    )
    
    fold_selection <- select_prototypes_from_importance(
      importance_df = fold_selector_importance,
      prototype_names = fold_proto_names,
      cumulative_target = CUMULATIVE_IMPORTANCE_TARGET
    )
    
    fold_selected_proto_names <- fold_selection$names
    
    cv_selected_proto_list[[fold]] <- fold_selection$selected %>%
      mutate(
        Fold = fold,
        Selected_Prototype_Count = length(fold_selected_proto_names)
      )
    
    cv_selector_importance_list[[fold]] <- fold_selector_importance %>%
      mutate(
        Fold = fold
      )
    
    fold_variant_features <- list(
      Raw64 = embedding_cols,
      Raw64_ClassSimilarity = c(
        embedding_cols,
        class_similarity_cols
      ),
      Raw64_ClassSimilarity_Margins = c(
        embedding_cols,
        class_similarity_cols,
        margin_cols
      ),
      Raw64_ClassSimilarity_Margins_AllPrototypes = c(
        embedding_cols,
        class_similarity_cols,
        margin_cols,
        fold_proto_names
      ),
      CRRF_Selected = c(
        embedding_cols,
        class_similarity_cols,
        margin_cols,
        fold_selected_proto_names
      ),
      EngineeredOnly_Selected = c(
        class_similarity_cols,
        margin_cols,
        fold_selected_proto_names
      )
    )
    
    variant_names <- names(
      fold_variant_features
    )
    
    for (variant_index in seq_along(variant_names)) {
      variant_name <- variant_names[variant_index]
      features <- fold_variant_features[[variant_name]]
      
      rf_result <- fit_rf(
        train_df = train_fold,
        test_df = test_fold,
        feature_names = features,
        classes = EXPECTED_CLASSES,
        ntree = RF_TREES_INTERNAL,
        seed = SEED_CV_BASE +
          fold * 100 +
          variant_index
      )
      
      pixel_metrics <- classification_metrics(
        truth = test_fold$class,
        prediction = rf_result$prediction,
        classes = EXPECTED_CLASSES
      )
      
      site_votes <- site_majority_predictions(
        site_key = test_fold$site_key,
        truth = test_fold$class,
        prediction = rf_result$prediction
      )
      
      site_metrics <- classification_metrics(
        truth = site_votes$Truth,
        prediction = site_votes$Predicted,
        classes = EXPECTED_CLASSES
      )
      
      target_metrics <- target_vs_rest_metrics(
        truth = test_fold$class,
        prediction = rf_result$prediction,
        target_class = CLASS_TARGET
      )
      
      cv_fold_summary_list[[cv_record_index]] <- tibble(
        Fold = fold,
        Variant = variant_name,
        N_Features = length(features),
        Selected_Prototype_Count = ifelse(
          variant_name %in% c(
            "CRRF_Selected",
            "EngineeredOnly_Selected"
          ),
          length(fold_selected_proto_names),
          NA_integer_
        ),
        Pixel_Accuracy = pixel_metrics$overall$Accuracy,
        Pixel_Balanced_Accuracy = pixel_metrics$overall$Balanced_Accuracy,
        Pixel_Macro_F1 = pixel_metrics$overall$Macro_F1,
        Pixel_Macro_IoU = pixel_metrics$overall$Macro_IoU,
        Site_Accuracy = site_metrics$overall$Accuracy,
        Site_Balanced_Accuracy = site_metrics$overall$Balanced_Accuracy,
        Site_Macro_F1 = site_metrics$overall$Macro_F1,
        Site_Macro_IoU = site_metrics$overall$Macro_IoU,
        Target_Precision = target_metrics$Precision,
        Target_Recall = target_metrics$Recall,
        Target_F1 = target_metrics$F1,
        Target_IoU = target_metrics$IoU
      )
      
      cv_fold_pixel_class_list[[cv_record_index]] <-
        pixel_metrics$per_class %>%
        mutate(
          Fold = fold,
          Variant = variant_name
        )
      
      cv_fold_site_class_list[[cv_record_index]] <-
        site_metrics$per_class %>%
        mutate(
          Fold = fold,
          Variant = variant_name
        )
      
      cv_fold_confusion_list[[cv_record_index]] <-
        pixel_metrics$confusion %>%
        mutate(
          Fold = fold,
          Variant = variant_name
        )
      
      cv_fold_site_votes_list[[cv_record_index]] <-
        site_votes %>%
        mutate(
          Fold = fold,
          Variant = variant_name
        )
      
      cv_prediction_list[[cv_record_index]] <- tibble(
        row_id = test_fold$row_id,
        site_key = test_fold$site_key,
        site_id = test_fold[[site_id_col]],
        class = test_fold$class,
        Fold = fold,
        Variant = variant_name,
        Prediction = rf_result$prediction
      )
      
      cv_record_index <- cv_record_index + 1
    }
  }
  
  cv_predictions <- bind_rows(
    cv_prediction_list
  )
  
  cv_fold_summary <- bind_rows(
    cv_fold_summary_list
  )
  
  cv_fold_pixel_class_metrics <- bind_rows(
    cv_fold_pixel_class_list
  )
  
  cv_fold_site_class_metrics <- bind_rows(
    cv_fold_site_class_list
  )
  
  cv_fold_confusions <- bind_rows(
    cv_fold_confusion_list
  )
  
  cv_fold_site_votes <- bind_rows(
    cv_fold_site_votes_list
  )
  
  cv_selected_prototypes <- bind_rows(
    cv_selected_proto_list
  )
  
  cv_selector_importance <- bind_rows(
    cv_selector_importance_list
  )
  
  
  # ==========================================================================
  # 4.21 POOLED OUT-OF-FOLD INTERNAL CV METRICS
  # ==========================================================================
  
  pooled_cv_summary_list <- list()
  pooled_cv_pixel_class_list <- list()
  pooled_cv_site_class_list <- list()
  pooled_cv_confusion_list <- list()
  pooled_cv_target_list <- list()
  
  cv_variant_names <- unique(
    cv_predictions$Variant
  )
  
  for (variant_index in seq_along(cv_variant_names)) {
    variant_name <- cv_variant_names[variant_index]
    
    variant_predictions <- cv_predictions %>%
      filter(
        Variant == variant_name
      ) %>%
      arrange(row_id)
    
    if (
      nrow(variant_predictions) != nrow(production_df)
    ) {
      stop(
        paste0(
          "Internal CV did not produce exactly one out-of-fold prediction per row for variant ",
          variant_name,
          "."
        )
      )
    }
    
    pooled_pixel_metrics <- classification_metrics(
      truth = variant_predictions$class,
      prediction = variant_predictions$Prediction,
      classes = EXPECTED_CLASSES
    )
    
    pooled_site_votes <- site_majority_predictions(
      site_key = variant_predictions$site_key,
      truth = variant_predictions$class,
      prediction = variant_predictions$Prediction
    )
    
    pooled_site_metrics <- classification_metrics(
      truth = pooled_site_votes$Truth,
      prediction = pooled_site_votes$Predicted,
      classes = EXPECTED_CLASSES
    )
    
    pooled_target_metrics <- target_vs_rest_metrics(
      truth = variant_predictions$class,
      prediction = variant_predictions$Prediction,
      target_class = CLASS_TARGET
    )
    
    pooled_cv_summary_list[[variant_index]] <- tibble(
      Variant = variant_name,
      Pixel_Accuracy = pooled_pixel_metrics$overall$Accuracy,
      Pixel_Balanced_Accuracy = pooled_pixel_metrics$overall$Balanced_Accuracy,
      Pixel_Macro_F1 = pooled_pixel_metrics$overall$Macro_F1,
      Pixel_Macro_IoU = pooled_pixel_metrics$overall$Macro_IoU,
      Site_Accuracy = pooled_site_metrics$overall$Accuracy,
      Site_Balanced_Accuracy = pooled_site_metrics$overall$Balanced_Accuracy,
      Site_Macro_F1 = pooled_site_metrics$overall$Macro_F1,
      Site_Macro_IoU = pooled_site_metrics$overall$Macro_IoU,
      Target_Precision = pooled_target_metrics$Precision,
      Target_Recall = pooled_target_metrics$Recall,
      Target_F1 = pooled_target_metrics$F1,
      Target_IoU = pooled_target_metrics$IoU
    )
    
    pooled_cv_pixel_class_list[[variant_index]] <-
      pooled_pixel_metrics$per_class %>%
      mutate(
        Variant = variant_name
      )
    
    pooled_cv_site_class_list[[variant_index]] <-
      pooled_site_metrics$per_class %>%
      mutate(
        Variant = variant_name
      )
    
    pooled_cv_confusion_list[[variant_index]] <-
      pooled_pixel_metrics$confusion %>%
      mutate(
        Variant = variant_name
      )
    
    pooled_cv_target_list[[variant_index]] <-
      pooled_target_metrics %>%
      mutate(
        Variant = variant_name
      )
  }
  
  pooled_cv_summary <- bind_rows(
    pooled_cv_summary_list
  )
  
  pooled_cv_pixel_class_metrics <- bind_rows(
    pooled_cv_pixel_class_list
  )
  
  pooled_cv_site_class_metrics <- bind_rows(
    pooled_cv_site_class_list
  )
  
  pooled_cv_confusions <- bind_rows(
    pooled_cv_confusion_list
  )
  
  pooled_cv_target_metrics <- bind_rows(
    pooled_cv_target_list
  )
  
  cat("\nPooled leakage-free internal CV summary:\n")
  print(pooled_cv_summary)
  
  
  # ==========================================================================
  # 4.22 FOLD-TO-FOLD INTERNAL CV SUMMARY
  # ==========================================================================
  
  cv_fold_variability <- cv_fold_summary %>%
    group_by(Variant) %>%
    summarise(
      Folds = n(),
      Pixel_Accuracy_Mean = mean(Pixel_Accuracy),
      Pixel_Accuracy_SD = sd(Pixel_Accuracy),
      Pixel_Macro_F1_Mean = mean(Pixel_Macro_F1),
      Pixel_Macro_F1_SD = sd(Pixel_Macro_F1),
      Target_F1_Mean = mean(Target_F1),
      Target_F1_SD = sd(Target_F1),
      Target_IoU_Mean = mean(Target_IoU),
      Target_IoU_SD = sd(Target_IoU),
      Site_Accuracy_Mean = mean(Site_Accuracy),
      Site_Accuracy_SD = sd(Site_Accuracy),
      .groups = "drop"
    )
  
  
  # ==========================================================================
  # 4.23 SITE-LEVEL BOOTSTRAP
  # ==========================================================================
  # Entire sites are resampled within class. Individual pixels are never
  # resampled independently.
  # ==========================================================================
  
  cat(
    "\nRunning",
    BOOTSTRAP_REPS,
    "site-level bootstrap replicates...\n"
  )
  
  target_df <- production_df %>%
    filter(class == CLASS_TARGET)
  
  ocean_df <- production_df %>%
    filter(class == CLASS_OCEAN)
  
  conf_df <- production_df %>%
    filter(class == CLASS_CONF)
  
  target_mat <- as.matrix(
    target_df[, embedding_cols, drop = FALSE]
  )
  
  target_mat_unit <- normalise_matrix_rows(
    target_mat,
    vector_label = "bootstrap target AEF"
  )
  
  ocean_mat <- as.matrix(
    ocean_df[, embedding_cols, drop = FALSE]
  )
  
  conf_mat <- as.matrix(
    conf_df[, embedding_cols, drop = FALSE]
  )
  
  target_indices_by_site <- split(
    seq_len(nrow(target_df)),
    target_df$site_key
  )
  
  ocean_indices_by_site <- split(
    seq_len(nrow(ocean_df)),
    ocean_df$site_key
  )
  
  conf_indices_by_site <- split(
    seq_len(nrow(conf_df)),
    conf_df$site_key
  )
  
  target_sites <- names(
    target_indices_by_site
  )
  
  ocean_sites <- names(
    ocean_indices_by_site
  )
  
  conf_sites <- names(
    conf_indices_by_site
  )
  
  bootstrap_metrics_list <- vector(
    "list",
    BOOTSTRAP_REPS
  )
  
  bootstrap_target_centroids <- matrix(
    NA_real_,
    nrow = BOOTSTRAP_REPS,
    ncol = length(embedding_cols)
  )
  
  bootstrap_ocean_centroids <- matrix(
    NA_real_,
    nrow = BOOTSTRAP_REPS,
    ncol = length(embedding_cols)
  )
  
  bootstrap_conf_centroids <- matrix(
    NA_real_,
    nrow = BOOTSTRAP_REPS,
    ncol = length(embedding_cols)
  )
  
  set.seed(SEED_BOOTSTRAP)
  
  for (b in seq_len(BOOTSTRAP_REPS)) {
    sampled_target_sites <- sample(
      target_sites,
      size = length(target_sites),
      replace = TRUE
    )
    
    sampled_ocean_sites <- sample(
      ocean_sites,
      size = length(ocean_sites),
      replace = TRUE
    )
    
    sampled_conf_sites <- sample(
      conf_sites,
      size = length(conf_sites),
      replace = TRUE
    )
    
    target_idx <- sample_site_indices(
      target_indices_by_site,
      sampled_target_sites
    )
    
    ocean_idx <- sample_site_indices(
      ocean_indices_by_site,
      sampled_ocean_sites
    )
    
    conf_idx <- sample_site_indices(
      conf_indices_by_site,
      sampled_conf_sites
    )
    
    boot_target_centroid <- normalise_vector(
      colMeans(
        target_mat[target_idx, , drop = FALSE]
      )
    )
    
    boot_ocean_centroid <- normalise_vector(
      colMeans(
        ocean_mat[ocean_idx, , drop = FALSE]
      )
    )
    
    boot_conf_centroid <- normalise_vector(
      colMeans(
        conf_mat[conf_idx, , drop = FALSE]
      )
    )
    
    bootstrap_target_centroids[b, ] <-
      boot_target_centroid
    
    bootstrap_ocean_centroids[b, ] <-
      boot_ocean_centroid
    
    bootstrap_conf_centroids[b, ] <-
      boot_conf_centroid
    
    bootstrap_target_similarity <- as.numeric(
      target_mat_unit[target_idx, , drop = FALSE] %*%
        boot_target_centroid
    )
    
    boot_thresholds <- as.numeric(
      quantile(
        bootstrap_target_similarity,
        probs = c(0.10, 0.05, 0.01),
        names = FALSE
      )
    )
    
    bootstrap_metrics_list[[b]] <- tibble(
      Replicate = b,
      Threshold_Retain90 = boot_thresholds[1],
      Threshold_Retain95 = boot_thresholds[2],
      Threshold_Retain99 = boot_thresholds[3],
      Target_Centroid_Stability = sum(
        boot_target_centroid *
          production_engineering$target_centroid
      ),
      Ocean_Centroid_Stability = sum(
        boot_ocean_centroid *
          production_engineering$ocean_centroid
      ),
      Conf_Centroid_Stability = sum(
        boot_conf_centroid *
          production_engineering$conf_centroid
      ),
      Target_vs_Ocean_Centroid_Similarity = sum(
        boot_target_centroid *
          boot_ocean_centroid
      ),
      Target_vs_Conf_Centroid_Similarity = sum(
        boot_target_centroid *
          boot_conf_centroid
      ),
      Ocean_vs_Conf_Centroid_Similarity = sum(
        boot_ocean_centroid *
          boot_conf_centroid
      )
    )
    
    if (b %% 50 == 0) {
      cat(
        "  Bootstrap",
        b,
        "/",
        BOOTSTRAP_REPS,
        "\n"
      )
    }
  }
  
  bootstrap_metrics <- bind_rows(
    bootstrap_metrics_list
  )
  
  bootstrap_summary <- bootstrap_metrics %>%
    pivot_longer(
      cols = -Replicate,
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    group_by(Metric) %>%
    summarise(
      Mean = mean(Value),
      SD = sd(Value),
      Q025 = as.numeric(quantile(Value, 0.025)),
      Median = median(Value),
      Q975 = as.numeric(quantile(Value, 0.975)),
      .groups = "drop"
    )
  
  target_boot_df <- as.data.frame(
    bootstrap_target_centroids
  )
  
  names(target_boot_df) <- paste0(
    "target_",
    embedding_cols
  )
  
  ocean_boot_df <- as.data.frame(
    bootstrap_ocean_centroids
  )
  
  names(ocean_boot_df) <- paste0(
    "ocean_",
    embedding_cols
  )
  
  conf_boot_df <- as.data.frame(
    bootstrap_conf_centroids
  )
  
  names(conf_boot_df) <- paste0(
    "conf_",
    embedding_cols
  )
  
  bootstrap_centroids <- bind_cols(
    tibble(
      Replicate = seq_len(BOOTSTRAP_REPS)
    ),
    target_boot_df,
    ocean_boot_df,
    conf_boot_df
  )
  
  
  # ==========================================================================
  # 4.24 GEE VECTOR TABLES
  # ==========================================================================
  # Vectors.csv contains only the selected production prototypes.
  # Vectors_AllPrototypes.csv retains all K = 5 prototypes for ablation and
  # reproducibility.
  # ==========================================================================
  
  gee_vectors <- data.frame(
    band = embedding_cols,
    farm = production_engineering$target_centroid,
    ocean = production_engineering$ocean_centroid,
    conf = production_engineering$conf_centroid,
    stringsAsFactors = FALSE
  )
  
  for (proto_name in selected_proto_names) {
    proto_index <- as.numeric(
      sub(
        "proto_",
        "",
        proto_name
      )
    )
    
    gee_vectors[[proto_name]] <-
      production_engineering$prototypes[
        proto_index,
      ]
  }
  
  gee_vectors_all <- data.frame(
    band = embedding_cols,
    farm = production_engineering$target_centroid,
    ocean = production_engineering$ocean_centroid,
    conf = production_engineering$conf_centroid,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(K)) {
    gee_vectors_all[[paste0("proto_", i)]] <-
      production_engineering$prototypes[i, ]
  }
  
  
  # ==========================================================================
  # 4.25 COMPACT DIAGNOSTICS
  # ==========================================================================
  
  pooled_raw64 <- pooled_cv_summary %>%
    filter(
      Variant == "Raw64"
    )
  
  pooled_crrf <- pooled_cv_summary %>%
    filter(
      Variant == "CRRF_Selected"
    )
  
  diagnostics <- tibble(
    Metric = c(
      "Total_Pixels",
      "Target_Pixels",
      "Ocean_Pixels",
      "Confounding_Pixels",
      "Target_Sites",
      "Ocean_Sites",
      "Confounding_Sites",
      "K",
      "Selected_Prototype_Count",
      "Internal_CV_Folds",
      "Internal_CV_Raw64_Target_F1",
      "Internal_CV_CRRF_Target_F1",
      "Internal_CV_Target_F1_Difference",
      "Internal_CV_Raw64_Target_IoU",
      "Internal_CV_CRRF_Target_IoU",
      "Internal_CV_Target_IoU_Difference",
      "Bootstrap_Reps"
    ),
    Value = c(
      nrow(production_df),
      sum(production_df$class == CLASS_TARGET),
      sum(production_df$class == CLASS_OCEAN),
      sum(production_df$class == CLASS_CONF),
      n_distinct(
        production_df$site_key[
          production_df$class == CLASS_TARGET
        ]
      ),
      n_distinct(
        production_df$site_key[
          production_df$class == CLASS_OCEAN
        ]
      ),
      n_distinct(
        production_df$site_key[
          production_df$class == CLASS_CONF
        ]
      ),
      K,
      length(selected_proto_names),
      INTERNAL_CV_FOLDS,
      pooled_raw64$Target_F1,
      pooled_crrf$Target_F1,
      pooled_crrf$Target_F1 - pooled_raw64$Target_F1,
      pooled_raw64$Target_IoU,
      pooled_crrf$Target_IoU,
      pooled_crrf$Target_IoU - pooled_raw64$Target_IoU,
      BOOTSTRAP_REPS
    )
  )
  
  
  # ==========================================================================
  # 4.26 PROVENANCE MANIFEST
  # ==========================================================================
  
  input_md5 <- as.character(
    tools::md5sum(input_file)
  )
  
  unique_or_na <- function(column_name) {
    if (!column_name %in% names(production_df)) {
      return(NA_character_)
    }
    
    values <- unique(
      as.character(production_df[[column_name]])
    )
    
    values <- values[
      !is.na(values) & values != ""
    ]
    
    if (length(values) == 0) {
      NA_character_
    } else {
      paste(values, collapse = "|")
    }
  }
  
  manifest <- tibble(
    Parameter = c(
      "Region",
      "Year",
      "Input_File",
      "Input_MD5",
      "Expected_Sites_Per_Class",
      "Expected_Radius_m",
      "Expected_Sampling_Scale_m",
      "Sampling_CRS",
      "Sampling_UTM_Zone",
      "Source_AEF_UTM_Zones",
      "AEF_Dataset_Versions",
      "AEF_Model_Versions",
      "AEF_Processing_Versions",
      "AEF_Dimensions",
      "K",
      "KMeans_NStart",
      "KMeans_Iter_Max",
      "Prototype_Cumulative_Importance_Target",
      "Production_Selector_RF_Trees",
      "Internal_CV_Folds",
      "Internal_RF_Trees",
      "Bootstrap_Reps",
      "Seed_KMeans",
      "Seed_Production_Selector",
      "Seed_CV_Folds",
      "Seed_CV_Base",
      "Seed_Bootstrap",
      "Selected_Prototypes",
      "R_Version",
      "randomForest_Version",
      "cluster_Version",
      "dplyr_Version"
    ),
    Value = c(
      region_name,
      YEAR,
      input_file,
      input_md5,
      EXPECTED_SITES_PER_CLASS,
      EXPECTED_RADIUS_M,
      EXPECTED_SAMPLING_SCALE_M,
      unique_or_na("sampling_crs"),
      unique_or_na("sampling_utm_zone"),
      unique_or_na("source_aef_utm_zones"),
      unique_or_na("aef_dataset_versions"),
      unique_or_na("aef_model_versions"),
      unique_or_na("aef_processing_versions"),
      length(embedding_cols),
      K,
      KMEANS_NSTART,
      KMEANS_ITER_MAX,
      CUMULATIVE_IMPORTANCE_TARGET,
      RF_TREES_PRODUCTION_DIAGNOSTIC,
      INTERNAL_CV_FOLDS,
      RF_TREES_INTERNAL,
      BOOTSTRAP_REPS,
      SEED_KMEANS,
      SEED_PRODUCTION_SELECTOR,
      SEED_CV_FOLDS,
      SEED_CV_BASE,
      SEED_BOOTSTRAP,
      paste(selected_proto_names, collapse = ";"),
      R.version.string,
      as.character(packageVersion("randomForest")),
      as.character(packageVersion("cluster")),
      as.character(packageVersion("dplyr"))
    )
  )
  
  
  # ==========================================================================
  # 4.27 EXPORT CSV OUTPUTS
  # ==========================================================================
  
  created_files <- c(
    created_files,
    write_region_output(
      production_df,
      region_name,
      "_Augmented.csv"
    ),
    write_region_output(
      gee_vectors,
      region_name,
      "_Vectors.csv"
    ),
    write_region_output(
      gee_vectors_all,
      region_name,
      "_Vectors_AllPrototypes.csv"
    ),
    write_region_output(
      centroid_archive,
      region_name,
      "_Centroids_Raw_Normalised.csv"
    ),
    write_region_output(
      centroid_relationships,
      region_name,
      "_ClassCentroidRelationships.csv"
    ),
    write_region_output(
      prototype_archive,
      region_name,
      "_Prototypes_Raw_Normalised.csv"
    ),
    write_region_output(
      site_qa,
      region_name,
      "_SiteQA.csv"
    ),
    write_region_output(
      pixel_support_summary,
      region_name,
      "_PixelSupportSummary.csv"
    ),
    write_region_output(
      embedding_norm_summary,
      region_name,
      "_EmbeddingNormQA.csv"
    ),
    write_region_output(
      site_mean_embeddings,
      region_name,
      "_SiteMeanEmbeddings.csv"
    ),
    write_region_output(
      site_mean_relationships,
      region_name,
      "_SiteMeanRelationships.csv"
    ),
    write_region_output(
      cosine_thresholds,
      region_name,
      "_CosineThresholds.csv"
    ),
    write_region_output(
      similarity_summary,
      region_name,
      "_SimilaritySummary.csv"
    ),
    write_region_output(
      overlap_diagnostics,
      region_name,
      "_SimilarityOverlap.csv"
    ),
    write_region_output(
      prototype_similarity_long,
      region_name,
      "_PrototypeSimilarities.csv"
    ),
    write_region_output(
      prototype_site_membership,
      region_name,
      "_PrototypeSiteMembership.csv"
    ),
    write_region_output(
      prototype_cluster_summary,
      region_name,
      "_PrototypeClusterSummary.csv"
    ),
    write_region_output(
      k_diagnostics,
      region_name,
      "_KMeansDiagnostics.csv"
    ),
    write_region_output(
      production_selector_importance,
      region_name,
      "_FeatureImportance_AllPrototypeSelector.csv"
    ),
    write_region_output(
      prototype_importance,
      region_name,
      "_PrototypeImportance.csv"
    ),
    write_region_output(
      selected_prototypes,
      region_name,
      "_SelectedPrototypes.csv"
    ),
    write_region_output(
      production_srf_importance,
      region_name,
      "_FeatureImportance_SRF.csv"
    ),
    write_region_output(
      production_crrf_importance,
      region_name,
      "_FeatureImportance_CRRF_Selected.csv"
    ),
    write_region_output(
      production_feature_manifest,
      region_name,
      "_ProductionFeatureManifest.csv"
    ),
    write_region_output(
      production_feature_sets_compact,
      region_name,
      "_ProductionFeatureSets.csv"
    ),
    write_region_output(
      fold_assignments,
      region_name,
      "_InternalCV_FoldAssignments.csv"
    ),
    write_region_output(
      cv_fold_summary,
      region_name,
      "_InternalCV_FoldMetrics.csv"
    ),
    write_region_output(
      cv_fold_variability,
      region_name,
      "_InternalCV_FoldVariability.csv"
    ),
    write_region_output(
      cv_predictions,
      region_name,
      "_InternalCV_OutOfFoldPredictions.csv"
    ),
    write_region_output(
      pooled_cv_summary,
      region_name,
      "_InternalCV_PooledSummary.csv"
    ),
    write_region_output(
      pooled_cv_pixel_class_metrics,
      region_name,
      "_InternalCV_PooledPixelClassMetrics.csv"
    ),
    write_region_output(
      pooled_cv_site_class_metrics,
      region_name,
      "_InternalCV_PooledSiteClassMetrics.csv"
    ),
    write_region_output(
      pooled_cv_target_metrics,
      region_name,
      "_InternalCV_PooledTargetMetrics.csv"
    ),
    write_region_output(
      pooled_cv_confusions,
      region_name,
      "_InternalCV_PooledConfusions.csv"
    ),
    write_region_output(
      cv_fold_pixel_class_metrics,
      region_name,
      "_InternalCV_FoldPixelClassMetrics.csv"
    ),
    write_region_output(
      cv_fold_site_class_metrics,
      region_name,
      "_InternalCV_FoldSiteClassMetrics.csv"
    ),
    write_region_output(
      cv_fold_confusions,
      region_name,
      "_InternalCV_FoldConfusions.csv"
    ),
    write_region_output(
      cv_fold_site_votes,
      region_name,
      "_InternalCV_FoldSiteVotes.csv"
    ),
    write_region_output(
      cv_selected_prototypes,
      region_name,
      "_InternalCV_SelectedPrototypes.csv"
    ),
    write_region_output(
      cv_selector_importance,
      region_name,
      "_InternalCV_SelectorImportance.csv"
    ),
    write_region_output(
      bootstrap_metrics,
      region_name,
      "_BootstrapMetrics.csv"
    ),
    write_region_output(
      bootstrap_summary,
      region_name,
      "_BootstrapSummary.csv"
    ),
    write_region_output(
      bootstrap_centroids,
      region_name,
      "_BootstrapCentroids.csv"
    ),
    write_region_output(
      diagnostics,
      region_name,
      "_Diagnostics.csv"
    ),
    write_region_output(
      manifest,
      region_name,
      "_Manifest.csv"
    )
  )
  
  
  # ==========================================================================
  # 4.28 SAVE COMPLETE ANALYSIS BUNDLE
  # ==========================================================================
  
  analysis_bundle <- list(
    region = region_name,
    year = YEAR,
    input_file = input_file,
    config = list(
      EXPECTED_SITES_PER_CLASS = EXPECTED_SITES_PER_CLASS,
      EXPECTED_RADIUS_M = EXPECTED_RADIUS_M,
      EXPECTED_SAMPLING_SCALE_M = EXPECTED_SAMPLING_SCALE_M,
      K = K,
      KMEANS_NSTART = KMEANS_NSTART,
      KMEANS_ITER_MAX = KMEANS_ITER_MAX,
      CUMULATIVE_IMPORTANCE_TARGET = CUMULATIVE_IMPORTANCE_TARGET,
      INTERNAL_CV_FOLDS = INTERNAL_CV_FOLDS,
      RF_TREES_INTERNAL = RF_TREES_INTERNAL,
      RF_TREES_PRODUCTION_DIAGNOSTIC = RF_TREES_PRODUCTION_DIAGNOSTIC,
      BOOTSTRAP_REPS = BOOTSTRAP_REPS
    ),
    embedding_cols = embedding_cols,
    site_qa = site_qa,
    production_engineering = production_engineering,
    production_selected_prototypes = selected_proto_names,
    production_prototype_importance = prototype_importance,
    cosine_thresholds = cosine_thresholds,
    production_feature_sets = production_feature_sets,
    production_selector_model = production_selector_model,
    production_srf_model = production_srf_model,
    production_crrf_model = production_crrf_model,
    internal_cv_fold_assignments = fold_assignments,
    internal_cv_predictions = cv_predictions,
    internal_cv_pooled_summary = pooled_cv_summary,
    internal_cv_selected_prototypes = cv_selected_prototypes,
    bootstrap_metrics = bootstrap_metrics,
    bootstrap_summary = bootstrap_summary,
    k_diagnostics = k_diagnostics
  )
  
  rds_path <- region_output_file(
    region_name,
    "_AnalysisBundle.rds"
  )
  
  saveRDS(
    analysis_bundle,
    rds_path
  )
  
  created_files <- c(
    created_files,
    rds_path
  )
  
  
  # ==========================================================================
  # 4.29 SAVE SESSION INFO
  # ==========================================================================
  
  session_path <- region_output_file(
    region_name,
    "_SessionInfo.txt"
  )
  
  writeLines(
    capture.output(
      sessionInfo()
    ),
    con = session_path
  )
  
  created_files <- c(
    created_files,
    session_path
  )
  
  
  # ==========================================================================
  # 4.30 REGION COMPLETE
  # ==========================================================================
  
  cat("\n------------------------------------------------------------\n")
  cat("REGION COMPLETE:", region_name, "\n")
  cat("------------------------------------------------------------\n")
  cat(
    "Selected production prototypes:",
    paste(selected_proto_names, collapse = ", "),
    "\n"
  )
  cat(
    "Production files for GEE:\n  ",
    region_output_file(region_name, "_Augmented.csv"),
    "\n  ",
    region_output_file(region_name, "_Vectors.csv"),
    "\n  ",
    region_output_file(region_name, "_Vectors_AllPrototypes.csv"),
    "\n"
  )
  cat(
    "Created",
    length(created_files),
    "files for",
    region_name,
    "\n"
  )
  
  tibble(
    Region = region_name,
    Input_Pixels = nrow(production_df),
    Target_Sites = n_distinct(
      production_df$site_key[
        production_df$class == CLASS_TARGET
      ]
    ),
    Ocean_Sites = n_distinct(
      production_df$site_key[
        production_df$class == CLASS_OCEAN
      ]
    ),
    Confounding_Sites = n_distinct(
      production_df$site_key[
        production_df$class == CLASS_CONF
      ]
    ),
    Selected_Prototypes = paste(
      selected_proto_names,
      collapse = ";"
    ),
    Internal_CV_Raw64_Target_F1 = pooled_raw64$Target_F1,
    Internal_CV_CRRF_Target_F1 = pooled_crrf$Target_F1,
    Internal_CV_Raw64_Target_IoU = pooled_raw64$Target_IoU,
    Internal_CV_CRRF_Target_IoU = pooled_crrf$Target_IoU
  )
}


# =============================================================================
# 5. RUN REQUESTED REGIONS
# =============================================================================

run_summary_list <- vector(
  "list",
  length(REGIONS_TO_RUN)
)

for (i in seq_along(REGIONS_TO_RUN)) {
  run_summary_list[[i]] <- process_region(
    REGIONS_TO_RUN[i]
  )
}

run_summary <- bind_rows(
  run_summary_list
)

run_summary_path <- file.path(
  OUTPUT_DIR,
  paste0(
    "AEF_",
    YEAR,
    "_AllRegions_RunSummary.csv"
  )
)

write.csv(
  run_summary,
  run_summary_path,
  row.names = FALSE
)


# =============================================================================
# 6. COMPLETE
# =============================================================================

cat("\n\n============================================================\n")
cat("ALL REQUESTED REGIONS COMPLETE\n")
cat("============================================================\n")
print(run_summary)
cat("\nCombined run summary:\n  ", run_summary_path, "\n")
cat("\nDone.\n")
