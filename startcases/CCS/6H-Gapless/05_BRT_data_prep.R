library(dplyr)
library(readr)
library(lubridate)

out_dir <- paste0("C:/github/Whale-SDM/", "startcases/CCS/6H-Gapless/Model-data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#load unfiltered master PA + track + var table
master_file <- file.path("C:/github/Whale-SDM/startcases/CCS/6H-Gapless/Model-data", "observed_crw_four_whales_static_dynamic_unfiltered.csv")

master_pts <- readr::read_csv(master_file, col_types = readr::cols(.default = readr::col_guess(), dTime = readr::col_datetime()), show_col_types = FALSE) %>%
mutate(whale_id = as.character(whale_id), track_id = as.character(track_id), choice_id = as.character(choice_id), PresAbs = as.integer(PresAbs))

prepare_brt_data <- function(data, predictors, min_available = 20)
{
  
required_columns <- c("whale_id", "track_id", "choice_id", "PresAbs", predictors)

#check for missing columns
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

#===================================
#1. Retain choice sets whose observed point has every predictor required by the candidate model
#===================================
valid_used_choices <- data %>%
    filter(
    PresAbs == 1L) %>%
    filter(if_all(all_of(predictors), ~ !is.na(.))) %>%
    distinct(whale_id, choice_id)
  
  
model_data <- data %>%
  semi_join(valid_used_choices, by = c("whale_id", "choice_id")) %>%
    
#Remove incomplete available alternatives within otherwise valid choice sets
  filter(if_all(all_of(predictors),
  ~ !is.na(.)
  )
)
  
#===================================
#2. Count valid alternatives remaining in each choice set
#===================================
choice_counts <- model_data %>%
  group_by(whale_id, choice_id) %>%
  summarise(n_used = sum(PresAbs == 1L),
  n_available = sum(PresAbs == 0L),
  .groups = "drop")

#Require at least one used location and adequate number of remaining alternatives
  valid_choices <- choice_counts %>%
  filter(n_used == 1L, n_available >= min_available)
  
  model_data <- model_data %>%
    semi_join(valid_choices, by = c("whale_id", "choice_id")) %>%
    left_join(valid_choices %>%
    select(whale_id, choice_id, n_available),
    by = c("whale_id", "choice_id")
    )

#===================================
#3. Count retained decisions per biological whale
#===================================
whale_counts <- valid_choices %>%
count(whale_id, name = "n_choices_whale")

model_data <- model_data %>%
left_join(whale_counts, by = "whale_id")

#===================================
#4. Calculate choice and whale weights
#===================================
model_data <- model_data %>%
  mutate(
# One total used unit and one total available unit
# within every choice set
      choice_class_weight = if_else(PresAbs == 1L,
        1,
        1 / n_available),

#Equal total contribution from each biological whale
      whale_balance_weight = 1 / n_choices_whale,
      
    model_weight_raw = choice_class_weight * whale_balance_weight,
      
#Model weight scaling to keep weight near their mean
      model_weight = model_weight_raw / mean(model_weight_raw),
#Same biological whale always receives the same fold.
#Both 254026 track segments therefore remain together.
    fold_whale = as.integer(factor(whale_id,
        levels = sort(unique(whale_id)
      )
    )
  )
)

#Maintain whale ID as metadata/fold info environmental predictor
#remember only GAMM can use ID as a random effect variable
model_data$whale_id <- factor(model_data$whale_id)
model_data$track_id <- factor(model_data$track_id)
  
#===================================
#5. Quality check/parsing for NA and anomalous values
#===================================

if (anyNA(model_data$model_weight)) {
stop("Model weights contain NA values.")
}

if (any(model_data$model_weight <= 0)) {
stop("Model weights must all be positive.")
}
  
return(model_data)
}

#prep candidate datasets ensembles/categories
static_predictors <- c("depth_m", "slope_deg", "rugosity_sd_2p5km_m", "dist_200m_isobath_km")
sst_predictors <- c("sst_c", "sst_front_position")
chla_predictors <- c("chla_mg_m3")
kd490_predictors <- c("kd490_m_inverse")
combined_predictors <- c(static_predictors, sst_predictors, chla_predictors, kd490_predictors)

#create datasets
#start with around 20 min_avaiable, try up to 25-50ish
brt_static_data <- prepare_brt_data(data = master_pts, predictors = static_predictors, min_available = 20)
brt_sst_data <- prepare_brt_data(data = master_pts, predictors = sst_predictors, min_available = 20)
brt_chla_data <- prepare_brt_data(data = master_pts, predictors = chla_predictors, min_available = 20)
brt_kd490_data <- prepare_brt_data(data = master_pts, predictors = kd490_predictors, min_available = 20)
brt_combined_data <- prepare_brt_data(data = master_pts, predictors = combined_predictors, min_available = 20)

#weighting check
weight_check <- brt_static_data %>%
  group_by(whale_id, PresAbs) %>%
  summarise(n_rows = n(),
    total_weight = sum(model_weight),
    .groups = "drop")

print(weight_check, width = Inf)

#retained choice test
brt_combined_data %>%
  distinct(whale_id, choice_id) %>%
  count(whale_id, name = "n_choices")

readr::write_rds(brt_static_data, file.path(out_dir, "brt_static_data.rds"))
readr::write_rds(brt_sst_data, file.path(out_dir, "brt_sst_data.rds"))
readr::write_rds(brt_chla_data, file.path(out_dir, "brt_chla_data.rds"))
readr::write_rds(brt_kd490_data, file.path(out_dir, "brt_kd490_data.rds"))
readr::write_rds(brt_combined_data, file.path(out_dir, "brt_combined_data.rds"))


#save BRT data for use in analysis script
#RDS used because preserves weight, datetime, factors
brt_data_dir <- file.path(out_dir, "BRT-ready")
dir.create(brt_data_dir, recursive = TRUE, showWarnings = FALSE)


brt_data_bundle <- list(
  data = list(
    static = brt_static_data,
    sst = brt_sst_data,
    chlorophyll = brt_chla_data,
    kd490 = brt_kd490_data,
    combined = brt_combined_data
  ),
  
  predictors = list(
    static = static_predictors,
    sst = sst_predictors,
    chlorophyll = chla_predictors,
    kd490 = kd490_predictors,
    combined = combined_predictors
  ),
  
settings = list(minimum_available_per_choice = 20,
#min choice subject to change
  weighting =paste0(
              "Observed weight = 1 per choice; ",
              "available weight totals 1 per choice; ",
              "biological whales have equal total weight; ",
              "weights scaled to mean 1."
              ),
  
validation = paste0("Leave-one-biological-whale-out folds; ", "both 254026 track segments remain in the same fold."),
  
omitted_predictors = c("sst_gradient_c_km")))

bundle_file <- file.path(brt_data_dir,  "brt_model_data_bundle.rds")

saveRDS(brt_data_bundle, bundle_file)

cat("\nSaved BRT data bundle to:\n", bundle_file, "\n")

#back up csv's
readr::write_csv(
  brt_static_data,
  file.path(
    brt_data_dir,
    "brt_static_model_data.csv"),
  na = "NA"
)

readr::write_csv(brt_sst_data, file.path(brt_data_dir, "brt_sst_model_data.csv"), na = "NA")
readr::write_csv(brt_chla_data, file.path(brt_data_dir, "brt_chlorophyll_model_data.csv"), na = "NA")
readr::write_csv(brt_kd490_data, file.path(brt_data_dir, "brt_kd490_model_data.csv"), na = "NA")
readr::write_csv(brt_combined_data, file.path(brt_data_dir, "brt_combined_model_data.csv"), na = "NA")

#dataset summary
make_brt_dataset_summary <- function(data, dataset_name)
{
  
choice_table <- data %>%
  group_by(whale_id, choice_id) %>%
  summarise(n_available = sum(PresAbs == 0L),
    .groups = "drop"
)
  
tibble(dataset = dataset_name, n_rows = nrow(data), n_observed = sum(data$PresAbs == 1L),
    n_available = sum(data$PresAbs == 0L),
    n_choice_sets = nrow(choice_table),
    n_whales = n_distinct(data$whale_id),
    min_available_per_choice = min(choice_table$n_available),
    median_available_per_choice = median(choice_table$n_available),
    mean_available_per_choice = mean(choice_table$n_available),
    max_available_per_choice =  max(choice_table$n_available)
)
}

brt_dataset_summary <- bind_rows(make_brt_dataset_summary(brt_static_data, "static"),
  make_brt_dataset_summary(brt_sst_data, "sst"),
  make_brt_dataset_summary(brt_chla_data, "chlorophyll"),
  make_brt_dataset_summary( brt_kd490_data, "kd490"),
  make_brt_dataset_summary(brt_combined_data, "combined")
)

print(brt_dataset_summary, width = Inf)
readr::write_csv(brt_dataset_summary, file.path(brt_data_dir, "brt_dataset_summary.csv"))