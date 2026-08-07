library(terra)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)

out_dir <- paste0("C:/github/Whale-SDM/", "startcases/CCS/6H/Model-data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#===================================
#1. load static-complete table
#===================================
static_file <- file.path("C:/github/Whale-SDM/startcases/CCS/6H/Model-data", "observed_crw_four_whales_static_complete.csv")

dynamic_model_pts <- readr::read_csv(static_file, col_types = readr::cols(.default = readr::col_guess(), dTime = readr::col_character()), show_col_types = FALSE) %>%

mutate(dTime = lubridate::ymd_hms(dTime,tz = "UTC", quiet = TRUE), date = as.Date(dTime))

if (!"date" %in% names(dynamic_model_pts)) {
  dynamic_model_pts$date <- as.Date(dynamic_model_pts$dTime)
}

if (anyNA(dynamic_model_pts$date)) {
  stop("At least one point has a missing date.")
}

#===================================
#2. Dynamic netCDF paths
#===================================
#Dynamic file paths
l3s_file <- paste0("C:/github/Whale-SDM/data/L3S-LEO/", "L3S_LEO_20240331_20240731_combined.nc")
viirs_chla_file <- paste0("C:/github/Whale-SDM/data/VIIRSNOAA20/", "NOAA20_VIIRS_chla_20240331_20240731_combined.nc")
viirs_kd490_file <- paste0("C:/github/Whale-SDM/data/VIIRS-kd490/", "VIIRS_Kd490_20240331_20240731_combined.nc")
dynamic_files <- c(L3S_LEO = l3s_file, VIIRS_chla = viirs_chla_file, VIIRS_Kd490 = viirs_kd490_file)

if (any(!file.exists(dynamic_files))) {
  stop(
    paste0(
      "Missing dynamic NetCDF file(s):\n",
      paste(
        dynamic_files[!file.exists(dynamic_files)],
        collapse = "\n"
      )
    )
  )
}

#===================================
#3. Reading the dynamic rasters
#===================================
#L3S contains three subdatasets to extract
l3s_ds <- terra::sds(l3s_file)

get_l3s_variable <- function(dataset, variable_name)
{
  
variable_index <- match(variable_name, names(dataset))
  
if (is.na(variable_index)) {
  stop("L3S variable not found: ", variable_name)
}

dataset[variable_index]
}

#L3S subdata
sst_r <- get_l3s_variable(l3s_ds, "sea_surface_temperature")
sst_gradient_r <- get_l3s_variable(l3s_ds, "sst_gradient_magnitude")
sst_front_r <- get_l3s_variable(l3s_ds, "sst_front_position")

#VIIRS chla and kd490 only have one variable
chla_r <- terra::rast( viirs_chla_file)
kd490_r <- terra::rast( viirs_kd490_file)

#Check data columns
availability_names <- c("l3s_date_available", "chla_date_available", "kd490_date_available")
cat("\nDate-availability columns present:\n")

print(availability_names %in% names(dynamic_model_pts))

cat("\nRows with unavailable product dates:\n")

print(dynamic_model_pts %>% summarise(
      l3s = sum(!l3s_date_available),
      chla = sum(!chla_date_available),
      kd490 = sum(!kd490_date_available))
)

#===================================
#5. Create point geometry 
#===================================
dynamic_points_sv <- terra::vect(dynamic_model_pts, geom = c( "long", "lat"), crs = "EPSG:4326")

#===================================
#6. Extract exact date to layer func
#===================================
#Key
#Missing raster date = NA
#Missing raster cell = NA
#No nearest-day substitution allowed
extract_daily_exact <- function(points_df, points_sv, raster_object, output_name)
{
raster_dates <- validate_daily_raster(raster_object, output_name)
output_values <- rep(NA_real_, nrow(points_df))
# Reproject points only if required
if (terra::same.crs(points_sv, raster_object)
) {

raster_points <- points_sv
} else {
  
raster_points <- terra::project(points_sv, terra::crs(raster_object))
}

matching_point_dates <- sort(unique(points_df$date[points_df$date %in% raster_dates]))

for (i in seq_along(matching_point_dates)) {
  current_date <- matching_point_dates[i]
  point_rows <- which(points_df$date == current_date)
  raster_layer <- match(
  current_date,
  raster_dates)
  
extracted_values <- terra::extract(raster_object[[raster_layer]], raster_points[point_rows, ], method = "simple", ID = FALSE)
output_values[point_rows] <- extracted_values[[1]]
}

points_df[[output_name]] <- output_values

message(output_name, ": matched ", length(matching_point_dates), " point dates; ", length(setdiff(unique(points_df$date), raster_dates)),
" required dates absent from product."
)

return(points_df)
}

#===================================
#7. Extract dynamic variables to track and PA
#===================================
dynamic_model_pts <- extract_daily_exact(points_df = dynamic_model_pts, points_sv = dynamic_points_sv, raster_object = sst_r, output_name = "sst_c")
dynamic_model_pts <- extract_daily_exact(points_df = dynamic_model_pts, points_sv = dynamic_points_sv, raster_object = sst_gradient_r, output_name = "sst_gradient_c_km")
dynamic_model_pts <- extract_daily_exact(points_df = dynamic_model_pts, points_sv = dynamic_points_sv, raster_object = sst_front_r, output_name = "sst_front_position")
dynamic_model_pts <- extract_daily_exact(points_df = dynamic_model_pts, points_sv = dynamic_points_sv, raster_object = chla_r, output_name = "chla_mg_m3")
dynamic_model_pts <- extract_daily_exact(points_df = dynamic_model_pts, points_sv = dynamic_points_sv, raster_object = kd490_r, output_name = "kd490_m_inverse")

#test to confirm existence of necessary columns
dynamic_names <- c("sst_c", "sst_gradient_c_km", "sst_front_position", "chla_mg_m3", "kd490_m_inverse")
cat("\nDynamic predictor columns present:\n")
print(dynamic_names
      %in%
      names(dynamic_model_pts))

cat("\nDynamic predictor summaries:\n")

print(summary(dynamic_model_pts[, dynamic_names]))
summary(dynamic_model_pts$chla_mg_m3)

quantile(dynamic_model_pts$chla_mg_m3, probs = c(
    0,
    0.01,
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    0.99,
    1
),
na.rm = TRUE
)

summary( dynamic_model_pts$sst_c, dynamic_model_pts$sst_gradient_c_km, dynamic_model_pts$sst_front_position, dynamic_model_pts$chla_mg_m3, dynamic_model_pts$kd490_m_inverse)

quantile(
dynamic_model_pts$chla_mg_m3,
probs = c(
    0,
    0.01,
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    0.99,
    1
),
na.rm = TRUE)

#===================================
#8. Calculate dynamic variable missingness
#===================================
dynamic_missingness_by_type <- dynamic_model_pts %>%
  group_by(PresAbs, location_type) %>%
  summarise(n_rows = n(),
    
across(all_of(dynamic_names),
  list(n_missing = ~ sum(is.na(.)),
  prop_missing = ~ mean(is.na(.)))),
  .groups = "drop")

print(dynamic_missingness_by_type, width = Inf)

#explicitly separate absent dates from unavailable cells
dynamic_missingness_source <- dynamic_model_pts %>% group_by(PresAbs) %>%
  summarise(n_rows = n(),
  
#L3S missingness (days)
l3s_missing_date = sum(!l3s_date_available, na.rm = TRUE),

#L3S  missingness on dates that exist
sst_missing_cell = sum(l3s_date_available & is.na(sst_c), na.rm = TRUE),
sst_gradient_missing_cell = sum(l3s_date_available & is.na(sst_gradient_c_km), na.rm = TRUE),
sst_front_missing_cell = sum(l3s_date_available & is.na(sst_front_position), na.rm = TRUE),

#Chla missingness
chla_missing_date = sum(!chla_date_available, na.rm = TRUE),
chla_missing_cell = sum(chla_date_available & is.na(chla_mg_m3), na.rm = TRUE),

#Kd490 missingness
kd490_missing_date = sum(!kd490_date_available, na.rm = TRUE),

kd490_missing_cell = sum(kd490_date_available & is.na(kd490_m_inverse), na.rm = TRUE), .groups = "drop")
print(dynamic_missingness_source,width = Inf)

#save  unfiltered dynamic table, master file
dynamic_out_file <- file.path(out_dir, "observed_crw_four_whales_static_dynamic_unfiltered.csv")
readr::write_csv(dynamic_model_pts, dynamic_out_file, na = "NA")
cat("\nSaved dynamic-extracted master table to:\n", dynamic_out_file, "\n")
