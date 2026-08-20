library(dplyr)
library(readr)
library(lubridate)
library(tidyr)
library(terra)
library(ncdf4)

#=================================
#1. READ COMBINED OBSERVED + CRW LOCATIONS
#=================================
setwd("C:/github/Whale-SDM/startcases/CCC/6H-Gapless/Model-data")
crw_dir <- paste0("C:/github/Whale-SDM/", "startcases/CCC/6H-Gapless/CRW")
pts_file <- file.path(crw_dir, "observed_and_crw_all_four_whales.csv")
out_dir <- paste0("C:/github/Whale-SDM/", "startcases/CCC/6H-Gapless/Model-data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pts_file)) {
  stop("Combined observed + CRW file not found:\n",pts_file)
}

#read dTime as character
model_pts <- readr::read_csv(
  pts_file,
  col_types = readr::cols(
  .default = readr::col_guess(),
  dTime = readr::col_character()
  ),
  show_col_types = FALSE
)


#Repair midnight timestamps, with write.csv() represented midnight timestamps as:
#2024-05-15
#Other times were represented as:
#2024-05-15 06:00:00
model_pts <- model_pts %>% mutate(dTime_raw = dTime,
    date_only_timestamp = grepl(
      "^\\d{4}-\\d{2}-\\d{2}$",
      dTime_raw),

  dTime_character = if_else(date_only_timestamp,
    paste0(
      dTime_raw,
      " 00:00:00"
    ),
    dTime_raw
  ),
  
  dTime = as.POSIXct(
    dTime_character,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )
)


#Confirm parsing

cat("\nNumber of repaired midnight rows:\n")

print(sum(model_pts$date_only_timestamp
  )
)

cat("\nNumber of failed timestamps after repair:\n")

print(sum(is.na(model_pts$dTime))
)

if (anyNA(model_pts$dTime)) {
  
print(
  model_pts %>%
    filter(
      is.na(dTime)
    ) %>%
    distinct(
      dTime_raw
    ) %>%
    head(30)
)

stop(
  "Some timestamp values still could not be parsed."
)
}

# ------------------------------------------------------------
# Standardize identifiers and convenience fields
# ------------------------------------------------------------

model_pts <- model_pts %>%
  mutate(
    long = as.numeric(long),
    lat = as.numeric(lat),
    
  #Continuous track segment
  track_id = as.character(whale_id),
  
  #whale's id
  whale_id = sub(
    "_.*$",
    "",
    track_id
  ),
  
  PresAbs = as.integer(used),
  
  date = as.Date(dTime),
  year = lubridate::year(dTime),
  month = lubridate::month(dTime),
  yday = lubridate::yday(dTime),
  hour = lubridate::hour(dTime)
) %>%
select(
  whale_id,
  track_id,
  choice_id,
  step_id,
  dTime,
  long,
  lat,
  PresAbs,
  sim_id,
  iteration,
  date,
  year,
  month,
  yday,
  hour,
  location_type,
  everything(),
  -used,
  -dTime_raw,
  -dTime_character,
  -date_only_timestamp
) %>%
arrange(
  whale_id,
  track_id,
  step_id,
  desc(PresAbs),
  iteration
)


#summarize/check strcutre
required_columns <- c(
  "whale_id",
  "track_id",
  "choice_id",
  "step_id",
  "dTime",
  "long",
  "lat",
  "PresAbs"
)

structural_na_summary <- model_pts %>%
  summarise(across(all_of(required_columns),
  ~ sum(is.na(.))
  )
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "column",
    values_to = "n_missing"
)
print(structural_na_summary)

#Check whales/track segments
cat("\nBiological whales:\n")
print(model_pts %>% distinct(
      whale_id
    )
)

cat("\nContinuous track segments:\n")

print(
  model_pts %>%
    distinct(whale_id, track_id
  ) %>%
    arrange(whale_id, track_id
  )
)


#Check each choice set
choice_summary <- model_pts %>%
  group_by(whale_id, track_id, choice_id
  ) %>%
  summarise(
    n_used = sum(PresAbs == 1L),
    n_available = sum(PresAbs == 0L),
    n_rows = n(),
    .groups = "drop"
  )
cat("\nChoice-set summary:\n")

print(summary(choice_summary$n_available
  )
)

bad_choices <- choice_summary %>%
  filter(
    n_used != 1L |
      n_available < 1L
  )

if (nrow(bad_choices) > 0) {
  
  warning(
    nrow(bad_choices),
    " choice sets do not contain exactly one used point ",
    "and at least one available point."
  )
  
  print(
    head(
      bad_choices,
      20
    )
  )
}


#count sets/rows

cat("\nRows by biological whale and location type:\n")

print(
  model_pts %>%
    count(whale_id, PresAbs,location_type
    )
)

cat("\nDate range:\n")
print(
  range(model_pts$dTime,
    na.rm = TRUE
  )
)
#=================================
#2. Calculate from GEBCO statis subvariables/derivatives
#=================================

gebco_file <- paste0("C:/github/Whale-SDM/data/GEBCO/", "gebco_2025_n45.0_s25.0_w-130.0_e-110.0.nc")

static_dir <- paste0("C:/github/Whale-SDM/data/GEBCO/", "derived_predictors")

dir.create(static_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(gebco_file)) {
  stop("GEBCO file not found:\n", gebco_file)
}


#Read GEBCO

gebco <- terra::rast(gebco_file)
cat("\nGEBCO raster:\n")
print(gebco)

cat("\nGEBCO layer names:\n")
print(names(gebco))

#Crop to  extent covering all observed and CRW points in the reduced CCC space
#Reduced Central California model domain
model_xlim <- c(-124.0, -121.0)
model_ylim <- c(35.5, 38.0)

# Extra area used ONLY for calculating static derivatives
static_buffer_deg <- 1.0

point_extent <- terra::ext(
  model_xlim[1] - static_buffer_deg,
  model_xlim[2] + static_buffer_deg,
  model_ylim[1] - static_buffer_deg,
  model_ylim[2] + static_buffer_deg
)

gebco_crop <- terra::crop(
  gebco,
  point_extent,
  snap = "out"
)

gebco_crop <- terra::crop(gebco, point_extent,snap = "out")

#Convert negative elevation to positive water depth
#GEBCO:
#ocean = negative
#land  = positive 
#Model predictor:
#water depth = positive meters
#land        = NA
depth_lonlat <- terra::ifel(
  gebco_crop < 0,
  -gebco_crop,
  NA)

names(depth_lonlat) <- "depth_m"

#Project to metric coordinates
# Monterey Bay and most of  domain are in UTM 10N.
#Set  regular 500-m output resolution
depth_utm <- terra::project(
  depth_lonlat,
  "EPSG:32610",
  method = "bilinear",
  res = 500)
names(depth_utm) <- "depth_m"


#Calculate Slope
slope_deg <- terra::terrain(
  depth_utm,
  v = "slope",
  unit = "degrees",
  neighbors = 8
)
names(slope_deg) <- "slope_deg"

#Calculate Local rugosity proxy
#5 cells × 500 m ~ 2.5-km-wide neighborhood.
#This is local depth SD, not a surface-area ratio.
rugosity_sd <- terra::focal(depth_utm, w = matrix(
    1,
    nrow = 5,
    ncol = 5
  ),
  fun = "sd",
  na.rm = TRUE,
  fillvalue = NA)
names(rugosity_sd) <- "rugosity_sd_2p5km_m"


#Create the 200-m isobath shelf proxy
shelf_200m <- terra::as.contour(depth_utm, levels = 200)

if (nrow(shelf_200m) == 0) {
  stop(
    "No 200-m contour was produced. ",
    "Check the GEBCO depth range and study extent."
)
}

#Ensure all derivatives exist at valid depth

#Distance from every cell to the 200-m isobath
dist_200m <- terra::distance(depth_utm, shelf_200m)

#Convert meters to kilometers
dist_200m <- dist_200m / 1000
names(dist_200m) <- "dist_200m_isobath_km"

#Apply same ocean mask to all derivatives so they dont land in land cells
slope_deg <- terra::mask(slope_deg, depth_utm)

rugosity_sd <- terra::mask(rugosity_sd, depth_utm)

dist_200m <- terra::mask(dist_200m, depth_utm)

#restore clear layer names after masking
names(depth_utm) <- "depth_m"
names(slope_deg) <- "slope_deg"
names(rugosity_sd) <- "rugosity_sd_2p5km_m"
names(dist_200m) <- "dist_200m_isobath_km"


#Combine derivative predictors

static_predictors <- c(depth_utm, slope_deg, rugosity_sd, dist_200m)

names(static_predictors) <- c(
  "depth_m",
  "slope_deg",
  "rugosity_sd_2p5km_m",
  "dist_200m_isobath_km"
)
#Save derivatives
terra::writeRaster(
  static_predictors,
  filename = file.path(
    static_dir,
    "GEBCO_static_predictors_500m.tif"
  ),
  overwrite = TRUE
)

terra::writeVector(
  shelf_200m,
  filename = file.path(
    static_dir,
    "GEBCO_200m_isobath.gpkg"
  ),
  overwrite = TRUE
)

cat("\nSaved GEBCO-derived predictors to:\n")
cat(static_dir, "\n")

#=================================
#3. EXTRACT STATIC PREDICTORS
#=================================
pts_sv <- terra::vect(model_pts, geom = c(
  "long",
  "lat"
  ),
 crs = "EPSG:4326"
)

pts_utm <- terra::project(pts_sv, terra::crs(static_predictors))


#mames of static predictor columns
static_names <- names(static_predictors)

#remove old ver  if this section is rerun interactively
model_pts <- model_pts %>%
select(
 -any_of(static_names))


#extract once
static_values <- terra::extract(
  static_predictors,
  pts_utm,
  ID = FALSE
)

model_pts <- bind_cols(
  model_pts,
  as.data.frame(static_values)
)


#Static predictor checks
cat("\nStatic predictor summaries:\n")
print(summary(model_pts[, static_names]))
cat("\nStatic predictor missingness:\n")

print(sapply(model_pts[, static_names], function(x) 
{mean(is.na(x))
  }
)
)

#Generate/save compelte staitic-predictor + PA/track table
#Identify choice sets whose observed location has at least one missing static predictor
bad_used_choice_ids <- model_pts %>%
filter(PresAbs == 1L, if_any(all_of(static_names), is.na)
  ) %>%
  distinct(choice_id) %>%
  pull(choice_id)

cat("\nObserved choice sets removed:\n")
print(length(bad_used_choice_ids))
print(bad_used_choice_ids)

#Removes:
#Entire choice tracks with incomplete observed point AND associated Incomplete CRW alternatives from otherwise valid sets
static_model_pts <- model_pts %>%
  filter(
!choice_id %in% bad_used_choice_ids
) %>%
filter(
  if_all(
    all_of(static_names),
     ~ !is.na(.)
  )
)


#Summarize retained choice sets after static filtering
static_choice_summary <- static_model_pts %>%
  group_by(
    whale_id,
    track_id,
    choice_id
  ) %>%
  summarise(
    n_used = sum(PresAbs == 1L),
    n_available = sum(PresAbs == 0L),
    available_retained_prop = n_available / 200,
    .groups = "drop"
  )

cat("\nNumber of retained observed choice sets:\n")
print(sum(static_model_pts$PresAbs == 1L))
cat("\nRetained CRW alternatives per choice set:\n")
print(summary(static_choice_summary$n_available))
cat("\nProportion of CRW alternatives retained:\n")
print(summary(static_choice_summary$available_retained_prop))

#Confirm every retained set still has one observed point
bad_static_choices <- static_choice_summary %>% filter(
  n_used != 1L |
   n_available < 1L
)
cat("\nInvalid choice sets after static filtering:\n")

print(nrow(bad_static_choices))

if (nrow(bad_static_choices) > 0) {
  print(bad_static_choices, n = 50, width = Inf)
}

#Confirm static predictors are now complete
cat("\nRemaining static predictor NAs:\n")
print(colSums( is.na(static_model_pts[, static_names])))

#Save static-complete table
static_out_file <- file.path(out_dir, "observed_crw_four_whales_static_complete.csv")
readr::write_csv(static_model_pts, static_out_file)
cat("\nSaved static-complete model table to:\n", static_out_file, "\n")

