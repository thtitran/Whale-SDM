library(dismo)
library(gbm)
library(dplyr)
library(readr)
library(tibble)
library(lubridate)

#load model data
model_data_dir <- paste0("C:/github/Whale-SDM/", "startcases/CCS/6H/Model-data")
bundle_file <- file.path(model_data_dir, "BRT-ready", "brt_model_data_bundle.rds")
brt_output_dir <- paste0("C:/github/Whale-SDM/", "startcases/CCS/6H/BRT-output")

#===================================
#1. Load model datasets + predictor lists from previous script
#===================================
brt_bundle <- readRDS(bundle_file)

brt_static_data <- brt_bundle$data$static

brt_sst_data <- brt_bundle$data$sst

brt_chla_data <- brt_bundle$data$chlorophyll

brt_kd490_data <- brt_bundle$data$kd490

brt_combined_data <- brt_bundle$data$combined

static_predictors <- brt_bundle$predictors$static

sst_predictors <- brt_bundle$predictors$sst

chla_predictors <- brt_bundle$predictors$chlorophyll

kd490_predictors <- brt_bundle$predictors$kd490

combined_predictors <- brt_bundle$predictors$combined

brt_bundle <- readRDS("C:/github/Whale-SDM/startcases/CCS/6H/Model-data/BRT-ready/brt_model_data_bundle.rds")

brt_static_data <- brt_bundle$data$static

static_predictors <- brt_bundle$predictors$static

#convert RDS to dataframe (gbm.step later cannot read tibbles, just dataframes)
brt_static_data <- as.data.frame(brt_bundle$data$static)
brt_sst_data <- as.data.frame(brt_bundle$data$sst)
brt_chla_data <- as.data.frame(brt_bundle$data$chlorophyll)
brt_kd490_data <- as.data.frame(brt_bundle$data$kd490)
brt_combined_data <- as.data.frame(brt_bundle$data$combined)

#retain the static template and add the dynamic predictors to produces these candidate models
#static
#static + SST/front
#static + chlorophyll
#static + Kd490
#static + SST/front + chlorophyll + Kd490
static_sst_predictors <- unique(
  c(static_predictors, sst_predictors))

static_chla_predictors <- unique(
  c(static_predictors, chla_predictors))

static_kd490_predictors <- unique(
  c(static_predictors, kd490_predictors))

#=========================
#2. weighted BRT fit func
#==========================
fit_weighted_brt <- function(
    model_data,
    predictors,
    model_name,
    tree_complexity = 2,
    learning_rate = 0.005,
    bag_fraction = 0.5,
    initial_trees = 500,
    step_size = 500,
    max_trees = 20000,
    seed = 123,
    plot_main = TRUE)
{
  
required_columns <- c("PresAbs", "model_weight", "fold_whale", predictors)

missing_columns <- setdiff(
  required_columns,
  names(model_data))

if (length(missing_columns) > 0) {
  stop(model_name, " is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
  
if (anyNA(model_data[, predictors])) {
  stop(model_name, " contains missing prediction values.")
}

if (anyNA(model_data$model_weight)) {
  stop(model_name, " contains missing model weights.")
}

if (any(model_data$model_weight <= 0)) {
  stop(model_name, " contains negative model weights.")
}

n_folds <- dplyr::n_distinct(model_data$fold_whale)

cat(
  "\n====================================================\n",
  "Fitting model: ", model_name, "\n",
  "Rows: ", nrow(model_data), "\n",
  "Observed points: ", sum(model_data$PresAbs == 1L), "\n",
  "Whale folds: ", n_folds, "\n",
  "Tree complexity: ", tree_complexity, "\n",
  "Learning rate: ", learning_rate, "\n",
  "Bag fraction: ", bag_fraction, "\n",
  "====================================================\n",
  sep = "")

#random seed for BRT
set.seed(seed)

fitted_model <- dismo::gbm.step(data = model_data,

gbm.x = match(predictors, names(model_data)),
gbm.y = match("PresAbs", names(model_data)),

family = "bernoulli",

#equal-whale + equal-choice weights
site.weights = model_data$model_weight,
  
#Leave one biological whale out in each fold validation
fold.vector = model_data$fold_whale,
n.folds = n_folds,

#folds  already manually defined
prev.stratify = FALSE,

tree.complexity = tree_complexity,

learning.rate = learning_rate,

bag.fraction = bag_fraction,

n.trees = initial_trees,

step.size = step_size,

max.trees = max_trees,

tolerance.method = "auto",

plot.main = plot_main,

plot.folds = FALSE,

verbose = TRUE,
silent = FALSE,

keep.fold.models = FALSE,
keep.fold.vector = TRUE,
keep.fold.fit = TRUE
)

if (is.null(fitted_model)) {
stop(model_name, " did not return a fitted model.")
}

fitted_model$model_name <- model_name

return(fitted_model)
}


#model performance extraction
#extract BRT metrics function that:
#saves
#validates output
#saves performance metrics
#generates part dependence plots
extract_brt_metrics <- function(model_object, model_name)
{

training_deviance_explained <- ((model_object$self.statistics$mean.null - model_object$self.statistics$mean.resid)
  /
  model_object$self.statistics$mean.null) * 100
  tibble::tibble(
  model = model_name,
  
  tree_complexity = model_object$gbm.call$tree.complexity,
  
  learning_rate = model_object$gbm.call$learning.rate,
  
  bag_fraction = model_object$gbm.call$bag.fraction,
  
  best_trees = model_object$gbm.call$best.trees,
  
  training_deviance_explained = training_deviance_explained,
  
  cv_deviance = model_object$cv.statistics$deviance.mean,
  
  cv_deviance_se = model_object$cv.statistics$deviance.se,
  
  training_auc =model_object$self.statistics$discrimination,
  
  cv_auc = model_object$cv.statistics$discrimination.mean,
  
  cv_auc_se = model_object$cv.statistics$discrimination.se)
}


#save BRT output
save_brt_outputs <- function(model_object, model_name, output_dir)
  {
#Validate input
  if (is.null(model_object)) {
    stop("Cannot save ", model_name, ": model object is NULL.")
}

if (length(model_name) != 1L || is.na(model_name) || model_name == "")
  {
  stop("model_name must be one non-empty character value.")
  }
  
#create BRT output directory
model_dir <- file.path(output_dir, model_name)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(model_dir)) {
    stop("Could not create model output directory:\n", model_dir)
  }

#Save  fitted model
model_file <- file.path(model_dir, paste0(model_name, ".rds"))
saveRDS(model_object, model_file)
  
#Save performance metrics
performance <- extract_brt_metrics(model_object, model_name)
performance_file <- file.path(model_dir, paste0(model_name, "_performance.csv"))
readr::write_csv(performance, performance_file)

#save variable importance
variable_importance <- summary(model_object, plotit = FALSE) %>%
tibble::as_tibble()
importance_file <- file.path(model_dir,paste0(model_name, "_variable_importance.csv"))
readr::write_csv(variable_importance, importance_file)
  
#Save partial dependence plots
partial_plot_file <- file.path(model_dir, paste0(model_name,"_partial_dependence.png"))
grDevices::png(filename = partial_plot_file, width = 14, height = 10, units = "in", res = 300, bg = "white")
tryCatch(
    {
    dismo::gbm.plot(model_object, smooth = TRUE, write.title = TRUE)
    },
  finally = {
  grDevices::dev.off()
}
)

cat("\nSaved outputs for ", model_name, " to:\n", model_dir, "\n", sep = "")
print(list.files(model_dir, full.names = TRUE))
invisible(model_dir)
}

{
training_deviance_explained <- ((model_object$self.statistics$mean.null - model_object$self.statistics$mean.resid) / model_object$self.statistics$mean.null) * 100
  
tibble(
  model = model_name,
  tree_complexity = model_object$gbm.call$tree.complexity,
  learning_rate = model_object$gbm.call$learning.rate,
  bag_fraction = model_object$gbm.call$bag.fraction,
  best_trees = model_object$gbm.call$best.trees,
  training_deviance_explained = training_deviance_explained,
  cv_deviance = model_object$cv.statistics$deviance.mean,
  cv_deviance_se = model_object$cv.statistics$deviance.se,
  training_auc = model_object$self.statistics$discrimination,
  cv_auc = model_object$cv.statistics$discrimination.mean,
  cv_auc_se = model_object$cv.statistics$discrimination.se)
}



#=================================================================
#STOP HERE, either use the following section to test between different BRT model parameters for the code to automatically
#pick the best performer or skip to the next section around line 393 for a chosen preset with test values
#=================================================================

#Use this model comparison for more/actual extensive model testing

#Two model comparison
#.005 not working
#brt_static_lr005 <- fit_weighted_brt(
#  model_data = brt_static_data,
# predictors = static_predictors,
#  model_name = "static_lr005",
#  tree_complexity = 2,
#  learning_rate = 0.005,
#  bag_fraction = 0.5,
#  initial_trees = 50,
#  step_size = 50,
#  max_trees = 20000,
#  plot_main = TRUE)

#this lr works better than .005
brt_static_lr001 <- fit_weighted_brt(
  model_data = brt_static_data,
  predictors = static_predictors,
  model_name = "static_lr001",
  tree_complexity = 2,
  learning_rate = 0.001,
  bag_fraction = 0.5,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = TRUE)

brt_static_lr010 <- fit_weighted_brt(
  model_data = brt_static_data,
  predictors = static_predictors,
  model_name = "static_lr010",
  tree_complexity = 2,
  learning_rate = 0.01,
  bag_fraction = 0.5,
  plot_main = TRUE)

#comparison test
static_tuning_results <- bind_rows(extract_brt_metrics(brt_static_lr001, "static_lr001"),
                                   extract_brt_metrics(brt_static_lr010, "static_lr010"))%>%
arrange(cv_deviance)

cat("\nStatic tuning comparison:\n")
print(static_tuning_results, width = Inf)
readr::write_csv(static_tuning_results, file.path( brt_output_dir, "static_initial_tuning_results.csv"))

#inspect  number of trees
cat("\nLR 0.005 best trees:", brt_static_lr005$gbm.call$best.trees, "\n")
cat("LR 0.01 best trees:", brt_static_lr010$gbm.call$best.trees, "\n")

#set to retain better of the two models:
if (brt_static_lr001$cv.statistics$deviance.mean <= brt_static_lr010$cv.statistics$deviance.mean)
{
  best_static_model <- brt_static_lr001
} else {
  best_static_model <- brt_static_lr010
}

#automatically selects the best parameters
chosen_tree_complexity <- best_static_model$gbm.call$tree.complexity
chosen_learning_rate <- best_static_model$gbm.call$learning.rate
chosen_bag_fraction <- best_static_model$gbm.call$bag.fraction
cat("\nChosen preliminary settings:\n", "Tree complexity: ", chosen_tree_complexity, "\n", "Learning rate: ", chosen_learning_rate, "\n", "Bag fraction: ", chosen_bag_fraction, "\n", sep = "")

#Save selected static model from automated selection
best_static_model <- brt_static_lr001t1
best_static_model$model_name <- "static_tc1_lr001"
save_brt_outputs(best_static_model, "static_tc1_lr001", brt_output_dir)



#===================================
#Use this manually changed parameter comparison for specific testing purposes
#===================================
#tree complexity models and comparison
brt_static_lr001t1 <- fit_weighted_brt(
  model_data = brt_static_data,
  predictors = static_predictors,
  model_name = "static_tc1_lr001",
  tree_complexity = 1,
  learning_rate = 0.001,
  bag_fraction = 0.5,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = TRUE)

brt_static_lr001t2 <- fit_weighted_brt(
  model_data = brt_static_data,
  predictors = static_predictors,
  model_name = "static_tc2_lr001",
  tree_complexity = 2,
  learning_rate = 0.001,
  bag_fraction = 0.5,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = TRUE)

brt_static_lr001t3 <- fit_weighted_brt(
  model_data = brt_static_data,
  predictors = static_predictors,
  model_name = "static_tc3_lr001",
  tree_complexity = 3,
  learning_rate = 0.001,
  bag_fraction = 0.5,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = TRUE)


#save tree complexity models
save_brt_outputs(brt_static_lr001t1, "static_tc1_lr001", brt_output_dir)
save_brt_outputs(brt_static_lr001t2, "static_tc2_lr001", brt_output_dir)
save_brt_outputs(brt_static_lr001t3, "static_tc3_lr001", brt_output_dir)

#use this instead for testing; these metrics work decent enough
chosen_tree_complexity <- 1
chosen_learning_rate   <- 0.001
chosen_bag_fraction    <- 0.5

#Save selected static model from manual testcase (in this case it was brt_static_lr001t1, the one with the tree complexity = 1)
best_static_model <- brt_static_lr001t1
best_static_model$model_name <- "static_tc1_lr001"
save_brt_outputs(best_static_model, "static_tc1_lr001", brt_output_dir)


# ------------------------------------------------------------
#Fit the dynamic candidate models; do only if static tunings are completed w/o error
# ------------------------------------------------------------
brt_static_sst <- fit_weighted_brt(
  model_data = brt_sst_data,
  predictors = static_sst_predictors,
  model_name = "static_plus_sst",
  tree_complexity = chosen_tree_complexity,
  learning_rate = chosen_learning_rate,
  bag_fraction = chosen_bag_fraction,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = TRUE
)

brt_static_chla <- fit_weighted_brt(
  model_data = brt_chla_data,
  predictors = static_chla_predictors,
  model_name = "static_plus_chla",
  tree_complexity = chosen_tree_complexity,
  learning_rate = chosen_learning_rate,
  bag_fraction = chosen_bag_fraction,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = TRUE
)

brt_static_kd490 <- fit_weighted_brt(
  model_data = brt_kd490_data,
  predictors = static_kd490_predictors,
  model_name = "static_plus_kd490",
  tree_complexity = chosen_tree_complexity,
  learning_rate = chosen_learning_rate,
  bag_fraction = chosen_bag_fraction,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = TRUE
)

brt_combined <- fit_weighted_brt(
  model_data = brt_combined_data,
  predictors = combined_predictors,
  model_name = "combined_all_predictors",
  tree_complexity = chosen_tree_complexity,
  learning_rate = chosen_learning_rate,
  bag_fraction = chosen_bag_fraction,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = TRUE
)

#save dynamic models
save_brt_outputs(brt_static_sst, "static_plus_sst", brt_output_dir)
save_brt_outputs(brt_static_chla, "static_plus_chla", brt_output_dir)
save_brt_outputs(brt_static_kd490, "static_plus_kd490", brt_output_dir)
save_brt_outputs(brt_combined, "combined_all_predictors", brt_output_dir)

#===========================================
#Fit matched static comparison models on each dynamic model's exact dataset if the full combined model has greatly
#reduced choices under observational conditions, ergo performance cant directly tell us if dynamic predictors improved model
#===========================================
brt_static_on_sst_sample <- fit_weighted_brt(
  model_data = brt_sst_data,
  predictors = static_predictors,
  model_name = "static_on_sst_sample",
  tree_complexity = chosen_tree_complexity,
  learning_rate = chosen_learning_rate,
  bag_fraction = chosen_bag_fraction,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = FALSE
)

brt_static_on_chla_sample <- fit_weighted_brt(
  model_data = brt_chla_data,
  predictors = static_predictors,
  model_name = "static_on_chla_sample",
  tree_complexity = chosen_tree_complexity,
  learning_rate = chosen_learning_rate,
  bag_fraction = chosen_bag_fraction,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = FALSE
)

brt_static_on_kd490_sample <- fit_weighted_brt(
  model_data = brt_kd490_data,
  predictors = static_predictors,
  model_name = "static_on_kd490_sample",
  tree_complexity = chosen_tree_complexity,
  learning_rate = chosen_learning_rate,
  bag_fraction = chosen_bag_fraction,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = FALSE
)

brt_static_on_combined_sample <- fit_weighted_brt(
  model_data = brt_combined_data,
  predictors = static_predictors,
  model_name = "static_on_combined_sample",
  tree_complexity = chosen_tree_complexity,
  learning_rate = chosen_learning_rate,
  bag_fraction = chosen_bag_fraction,
  initial_trees = 50,
  step_size = 50,
  max_trees = 50000,
  plot_main = FALSE
)

#save each dynamic-static matched model
save_brt_outputs(brt_static_on_sst_sample, "static_on_sst_sample", brt_output_dir)
save_brt_outputs(brt_static_on_chla_sample, "static_on_chla_sample", brt_output_dir)
save_brt_outputs(brt_static_on_kd490_sample, "static_on_kd490_sample", brt_output_dir)
save_brt_outputs(brt_static_on_combined_sample, "static_on_combined_sample", brt_output_dir)


#full candidate-model comparison
candidate_model_comparison <- dplyr::bind_rows(
extract_brt_metrics(best_static_model, "static_full_sample_reference"),
extract_brt_metrics(brt_static_on_sst_sample, "static_on_sst_sample"),
extract_brt_metrics(brt_static_sst, "static_plus_sst"),
extract_brt_metrics(brt_static_on_chla_sample, "static_on_chla_sample"),
extract_brt_metrics(brt_static_chla, "static_plus_chla"),
extract_brt_metrics(brt_static_on_kd490_sample, "static_on_kd490_sample"),
extract_brt_metrics(brt_static_kd490, "static_plus_kd490"),
extract_brt_metrics(brt_static_on_combined_sample, "static_on_combined_sample"),
extract_brt_metrics(brt_combined, "combined_all_predictors"
 )
)

#print and show model comparisons
print(candidate_model_comparison,width = Inf)
readr::write_csv(candidate_model_comparison, file.path(brt_output_dir, "candidate_model_performance.csv"))
