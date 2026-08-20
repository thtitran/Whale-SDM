library(dplyr)
library(readr)
library(tibble)
library(dismo)
library(gbm)
library(ggplot2)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

model_data_dir <- paste0("C:/github/Whale-SDM/", "startcases/CCS/6H-Gapless/Model-data")

bundle_file <- file.path(model_data_dir, "BRT-ready", "brt_model_data_bundle.rds")

brt_output_dir <- paste0("C:/github/Whale-SDM/", "startcases/CCS/6H-Gapless/BRT-output")


# ------------------------------------------------------------
# Reload prepared BRT datasets
# ------------------------------------------------------------

brt_bundle <- readRDS(
  bundle_file
)

# Convert to ordinary data.frames for gbm.step
brt_static_data <- as.data.frame(
  brt_bundle$data$static
)

brt_sst_data <- as.data.frame(
  brt_bundle$data$sst
)

brt_chla_data <- as.data.frame(
  brt_bundle$data$chlorophyll
)

brt_kd490_data <- as.data.frame(
  brt_bundle$data$kd490
)

brt_combined_data <- as.data.frame(
  brt_bundle$data$combined
)

# ------------------------------------------------------------
# Reload predictor lists
# ------------------------------------------------------------
static_predictors <- brt_bundle$predictors$static

sst_predictors <- brt_bundle$predictors$sst

chla_predictors <- brt_bundle$predictors$chlorophyll

kd490_predictors <- brt_bundle$predictors$kd490

combined_predictors <- brt_bundle$predictors$combined

# Recreate static + dynamic predictor combinations
static_sst_predictors <- unique(
  c(
    static_predictors,
    sst_predictors
  )
)

static_chla_predictors <- unique(
  c(
    static_predictors,
    chla_predictors
  )
)

static_kd490_predictors <- unique(
  c(
    static_predictors,
    kd490_predictors
  )
)


c(
  static_rows = nrow(brt_static_data),
  sst_rows = nrow(brt_sst_data),
  chla_rows = nrow(brt_chla_data),
  kd490_rows = nrow(brt_kd490_data),
  combined_rows = nrow(brt_combined_data)
)



sapply(
  list(
    static = brt_static_data,
    sst = brt_sst_data,
    chla = brt_chla_data,
    kd490 = brt_kd490_data,
    combined = brt_combined_data
  ),
  function(x) {
    c(
      rows = nrow(x),
      used = sum(x$PresAbs == 1L),
      whales = dplyr::n_distinct(x$whale_id),
      folds = dplyr::n_distinct(x$fold_whale)
    )
  }
)

brt_output_dir <- "C:/github/Whale-SDM/startcases/CCS/6H-gapless/BRT-output"

saveRDS(
  list(
    static = brt_static_data,
    sst = brt_sst_data,
    chla = brt_chla_data,
    kd490 = brt_kd490_data,
    combined = brt_combined_data
  ),
  file.path(
    brt_output_dir,
    "brt_prepared_datasets.rds"
  )
)


# ------------------------------------------------------------
# Reload fitted BRT models
# ------------------------------------------------------------

brt_static_lr001t1 <- readRDS(
  file.path(
    brt_output_dir,
    "static_tc1_lr001",
    "static_tc1_lr001.rds"
  )
)

brt_static_sst <- readRDS(
  file.path(
    brt_output_dir,
    "static_plus_sst",
    "static_plus_sst.rds"
  )
)

brt_static_chla <- readRDS(
  file.path(
    brt_output_dir,
    "static_plus_chla",
    "static_plus_chla.rds"
  )
)

brt_static_kd490 <- readRDS(
  file.path(
    brt_output_dir,
    "static_plus_kd490",
    "static_plus_kd490.rds"
  )
)

brt_combined <- readRDS(
  file.path(
    brt_output_dir,
    "combined_all_predictors",
    "combined_all_predictors.rds"
  )
)

brt_static_on_sst_sample <- readRDS(
  file.path(
    brt_output_dir,
    "static_on_sst_sample",
    "static_on_sst_sample.rds"
  )
)

brt_static_on_chla_sample <- readRDS(
  file.path(
    brt_output_dir,
    "static_on_chla_sample",
    "static_on_chla_sample.rds"
  )
)

brt_static_on_kd490_sample <- readRDS(
  file.path(
    brt_output_dir,
    "static_on_kd490_sample",
    "static_on_kd490_sample.rds"
  )
)

brt_static_on_combined_sample <- readRDS(
  file.path(
    brt_output_dir,
    "static_on_combined_sample",
    "static_on_combined_sample.rds"
  )
)

best_static_model <- brt_static_lr001t1

model_check <- c(
  static = exists("best_static_model"),
  sst = exists("brt_static_sst"),
  chla = exists("brt_static_chla"),
  kd490 = exists("brt_static_kd490"),
  combined = exists("brt_combined"),
  static_sst_match = exists("brt_static_on_sst_sample"),
  static_chla_match = exists("brt_static_on_chla_sample"),
  static_kd490_match = exists("brt_static_on_kd490_sample"),
  static_combined_match = exists("brt_static_on_combined_sample")
)
print(model_check)

# ------------------------------------------------------------
# per-whale held-out performance
# ------------------------------------------------------------
c(static_has_fold_fit = "fold.fit" %in% names(best_static_model), combined_has_fold_fit = "fold.fit" %in% names(brt_combined),
  matched_static_has_fold_fit =
    "fold.fit" %in% names(brt_static_on_combined_sample))

c(static_has_fold_vector = "fold.vector" %in% names(best_static_model), combined_has_fold_vector = "fold.vector" %in% names(brt_combined),
  matched_static_has_fold_vector =
    "fold.vector" %in% names(brt_static_on_combined_sample))

c(static_predictions = length(best_static_model$fold.fit), static_rows = nrow(brt_static_data), combined_predictions = length(brt_combined$fold.fit),
  combined_rows =
    nrow(brt_combined_data))

#Convert the withheld-whale predictions to probability
static_cv_predictions <- plogis(best_static_model$fold.fit)
combined_cv_predictions <- plogis(brt_combined$fold.fit)
matched_static_cv_predictions <- plogis(brt_static_on_combined_sample$fold.fit)

#sanity check (all should be btween 0-1)
summary(combined_cv_predictions)
range(combined_cv_predictions, na.rm = TRUE)

#gbm.step(), where AUC is already calculated separately for each CV fold.
#this step computes each fold's AUC from predictions made on that withheld fold.
#Since one fold = one biological whale in your setup, it immediately map those AUCs back to whales
combined_fold_map <- brt_combined_data %>%
  dplyr::distinct(fold_whale, whale_id) %>%
  dplyr::arrange(fold_whale)

#combined AUC
combined_auc_by_whale <- combined_fold_map %>%
  dplyr::mutate(cv_auc = brt_combined$cv.roc.matrix)

print(combined_auc_by_whale)

#static AUC
matched_static_auc_by_whale <- combined_fold_map %>%
  dplyr::mutate(cv_auc = brt_static_on_combined_sample$cv.roc.matrix)

print(matched_static_auc_by_whale)

#Combine dynamic and matched-static AUCs by whale
auc_by_whale_comparison <- combined_auc_by_whale %>%
  dplyr::rename(
    combined_cv_auc = cv_auc
  ) %>%
  dplyr::left_join(
    matched_static_auc_by_whale %>%
      dplyr::rename(
        matched_static_cv_auc = cv_auc
      ),
    by = c(
      "fold_whale",
      "whale_id"
    )
  ) %>%
  dplyr::mutate(
    delta_cv_auc =
      combined_cv_auc - matched_static_cv_auc)

# Display result
auc_by_whale_comparison %>%
  tibble::as_tibble() %>%
  print(width = Inf)

# Save result
readr::write_csv(auc_by_whale_comparison, file.path(brt_output_dir, "auc_by_whale_comparison.csv"))

# Check structure
class(auc_by_whale_comparison)

# ------------------------------------------------------------
#choice rank overall
# ------------------------------------------------------------
#Out-of-fold predictions for combined model

combined_choice_eval <- brt_combined_data %>%
  dplyr::mutate(cv_prediction = plogis( brt_combined$fold.fit))
# ------------------------------------------------------------
# Rank observed location within each choice set
# ------------------------------------------------------------

choice_rank_results <- combined_choice_eval %>%
  dplyr::group_by(whale_id, choice_id) %>%
  dplyr::summarise(n_available = sum(PresAbs == 0L), observed_prediction = cv_prediction[ PresAbs == 1L][1],
                   
                   observed_percentile = {
                     used_pred <-
                       cv_prediction[
                         PresAbs == 1L][1]
                     available_pred <-
                       cv_prediction[
                         PresAbs == 0L
                       ]
                     (sum(
                       available_pred < used_pred) + 0.5 *
                         sum(available_pred == used_pred)
                     ) /
                       length(available_pred)
                   },
                   .groups = "drop")

#summarize by whale
choice_rank_by_whale <- choice_rank_results %>%
  dplyr::group_by(whale_id
  ) %>%
  dplyr::summarise(n_choices = dplyr::n(), mean_percentile =
                     mean(observed_percentile, na.rm = TRUE),
                   
                   median_percentile =median(observed_percentile, na.rm = TRUE),
                   
                   prop_above_50 =mean(observed_percentile > 0.50, na.rm = TRUE),
                   
                   prop_above_75 = mean(observed_percentile > 0.75, na.rm = TRUE),
                   
                   prop_above_90 =mean(observed_percentile > 0.90, na.rm = TRUE),
                   
                   .groups = "drop")

print(choice_rank_by_whale)

# Save result
readr::write_csv(choice_rank_by_whale, file.path(brt_output_dir, "choice_rank_by_whale.csv"))


#overall result
choice_rank_overall <- choice_rank_results %>%
  dplyr::summarise(n_choices =  dplyr::n(),
                   mean_percentile = mean(observed_percentile, na.rm = TRUE),
                   median_percentile = median(observed_percentile, na.rm = TRUE),
                   prop_above_50 = mean(observed_percentile > 0.50, na.rm = TRUE),
                   prop_above_75 = mean(observed_percentile > 0.75, na.rm = TRUE),
                   prop_above_90 = mean(observed_percentile > 0.90, na.rm = TRUE)
  )

print(choice_rank_overall)

# Save result
readr::write_csv(choice_rank_overall, file.path(brt_output_dir, "choice_rank_overall.csv"))


# ------------------------------------------------------------
# Dynamic complete-case retention by whale
# ------------------------------------------------------------
static_choices_by_whale <- brt_static_data %>%
  dplyr::filter(PresAbs == 1L) %>%
  dplyr::distinct(whale_id, choice_id) %>%
  dplyr::count(
    whale_id,
    name = "static_choices")

combined_choices_by_whale <- brt_combined_data %>%
  dplyr::filter(PresAbs == 1L) %>%
  dplyr::distinct(whale_id, choice_id) %>%
  dplyr::count(
    whale_id,
    name = "combined_choices")

coverage_by_whale <- static_choices_by_whale %>%
  dplyr::left_join(
    combined_choices_by_whale,
    by = "whale_id"
  ) %>%
  dplyr::mutate(
    combined_choices =
      dplyr::coalesce(
        combined_choices,
        0L
      ),
    
    retained_fraction =
      combined_choices /
      static_choices,
    
    retained_percent =
      retained_fraction * 100)

print(coverage_by_whale)

# Save result
readr::write_csv(coverage_by_whale, file.path(brt_output_dir, "coverage_by_whale.csv"))


#Mark which original choices survive into combined model
combined_retained_choices <- brt_combined_data %>%
  dplyr::filter(PresAbs == 1L) %>%
  dplyr::distinct(whale_id, choice_id)  %>%
  dplyr::mutate(combined_retained = TRUE)


all_choice_retention <- brt_static_data %>%
  dplyr::filter(PresAbs == 1L) %>%
  dplyr::distinct(whale_id, track_id, choice_id, dTime) %>%
  dplyr::left_join(combined_retained_choices, by = c("whale_id", "choice_id")) %>%
  dplyr::mutate(combined_retained = dplyr::coalesce(combined_retained, FALSE), date = as.Date(dTime))

#Full versus retained temporal span/coverage
temporal_coverage_by_whale <- all_choice_retention %>%
  dplyr::group_by(whale_id) %>%
  dplyr::summarise(full_start = min(date, na.rm = TRUE),
                   full_end = max(date, na.rm = TRUE),
                   retained_start =
                     min(date[combined_retained], na.rm = TRUE),
                   
                   retained_end = max(date[combined_retained], na.rm = TRUE),
                   total_choices = dplyr::n(),
                   retained_choices = sum(combined_retained),
                   .groups = "drop")

# Save result
readr::write_csv(temporal_coverage_by_whale, file.path(brt_output_dir, "temporal_coverage_by_whale.csv"))

print(temporal_coverage_by_whale)

#Weekly temporal retention
weekly_retention <- all_choice_retention %>%
  dplyr::mutate(week = lubridate::floor_date(date, unit = "week")) %>%
  dplyr::group_by(whale_id, week) %>%
  dplyr::summarise(total_choices = dplyr::n(),
                   retained_choices =
                     sum(combined_retained),
                   retained_percent = 100 * retained_choices / total_choices,
                   .groups = "drop")

print(weekly_retention,n = Inf)

#Save result
readr::write_csv(weekly_retention, file.path(brt_output_dir, "weekly_retention.csv"))


#spatial extent check
grep("lon|lat|longitude|latitude", names(brt_static_data), value = TRUE, ignore.case = TRUE)

#Spatial retention of observed choices
all_choice_retention_spatial <- brt_static_data %>%
  dplyr::filter(PresAbs == 1L) %>%
  dplyr::distinct(whale_id, track_id,  choice_id,  dTime,  long,  lat) %>%
  dplyr::left_join(combined_retained_choices, by = c("whale_id", "choice_id")) %>%
  dplyr::mutate(combined_retained = dplyr::coalesce(combined_retained, FALSE))

ggplot2::ggplot(
  all_choice_retention_spatial,
  ggplot2::aes(x = long, y = lat)
) +
  ggplot2::geom_path(
    ggplot2::aes(group = track_id),
    linewidth = 0.4,
    alpha = 0.35
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(
      all_choice_retention_spatial,
      !combined_retained
    ),
    color = "grey70",
    size = 0.5,
    alpha = 0.7
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(
      all_choice_retention_spatial,
      combined_retained
    ),
    color = "red",
    size = 0.7
  ) +
  ggplot2::facet_wrap(~ whale_id) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    title = "Spatial distribution of choices retained in combined BRT",
    subtitle = "Grey = lost due to dynamic-data completeness; red = retained"
  ) +
  ggplot2::theme_bw()

ggsave("Spatial_distribution_of retained_choices.png")

#calculate whether the retained points occupy approximately the same geographic extent as the full observations
spatial_extent_summary <- all_choice_retention_spatial %>%
  dplyr::group_by(
    whale_id,
    combined_retained
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    
    lon_min = min(long, na.rm = TRUE),
    lon_max = max(long, na.rm = TRUE),
    lon_median = median(long, na.rm = TRUE),
    
    lat_min = min(lat, na.rm = TRUE),
    lat_max = max(lat, na.rm = TRUE),
    lat_median = median(lat, na.rm = TRUE),
    
    .groups = "drop"
  )


#Save result
readr::write_csv(spatial_extent_summary, file.path(brt_output_dir, "spatial_extent_summary.csv"))



print(
  spatial_extent_summary,
  width = Inf
)




#static-environment comparison
retention_environment_summary <- all_choice_retention_spatial %>%
  dplyr::left_join(
    brt_static_data %>%
      dplyr::filter(PresAbs == 1L) %>%
      dplyr::select(
        whale_id,
        choice_id,
        depth_m,
        slope_deg,
        rugosity_sd_2p5km_m,
        dist_200m_isobath_km
      ),
    by = c("whale_id", "choice_id")
  ) %>%
  dplyr::group_by(
    whale_id,
    combined_retained
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    
    median_depth =
      median(depth_m, na.rm = TRUE),
    
    median_slope =
      median(slope_deg, na.rm = TRUE),
    
    median_rugosity =
      median(rugosity_sd_2p5km_m, na.rm = TRUE),
    
    median_dist_200m =
      median(
        dist_200m_isobath_km,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


#Save result
readr::write_csv(retention_environment_summary, file.path(brt_output_dir, "retention_environment_summary.csv"))

print(
  retention_environment_summary,
  width = Inf
)


