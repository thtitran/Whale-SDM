library(dplyr)
library(readr)
library(lubridate)

# 1. load files
obs_file <- "C:/github/Whale-SDM/tracks/outputnewgamma/254026-RawArgos_pred_6h.csv"
pa_file  <- "C:/github/Whale-SDM/Output/CRW_test/crw_sim_254026_0.csv"

out_dir  <- "C:/github/Whale-SDM/Output/GAMM_prep"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, "254026_pres_pa_combined_GAMMready.csv")

# 2. Read files
obs <- read.csv(obs_file, stringsAsFactors = FALSE)
pa  <- read.csv(pa_file,  stringsAsFactors = FALSE)

# Quick checks
cat("Observed track dimensions:\n")
print(dim(obs))
cat("\nPseudo-absence dimensions:\n")
print(dim(pa))

cat("\nObserved columns:\n")
print(names(obs))
cat("\nPseudo-absence columns:\n")
print(names(pa))

# Clean and standardize observed data
#
# We only need the core fields for now.
obs2 <- obs %>%
  transmute(
    id        = id,
    dTime     = ymd_hms(date, tz = "UTC"),
    long      = as.numeric(lon),
    lat       = as.numeric(lat),
    iteration = NA_integer_,   # observed points are not simulated iterations
    PresAbs   = 1L
  )

# Clean and standardize pseudo-absence data
# Your PA file currently has:
# iteration, dTime, long, lat, flag
#
#  assign the same whale ID so later to know which whale these PA belong to

pa2 <- pa %>%
  transmute(
    id        = "254026",
    dTime     = ymd_hms(dTime, tz = "UTC"),
    long      = as.numeric(long),
    lat       = as.numeric(lat),
    iteration = as.integer(iteration),
    PresAbs   = 0L
  )


# 5. Optional sanity checks
cat("\nNumber of observed points:\n")
print(nrow(obs2))

cat("\nNumber of pseudo-absence points:\n")
print(nrow(pa2))

cat("\nUnique pseudo-absence iterations:\n")
print(length(unique(pa2$iteration)))

cat("\nObserved date range:\n")
print(range(obs2$dTime, na.rm = TRUE))

cat("\nPseudo-absence date range:\n")
print(range(pa2$dTime, na.rm = TRUE))

# 6 Combine presence + pseudo-absence
model_pts <- bind_rows(obs2, pa2) %>%
  arrange(dTime, PresAbs)

# 7. Add convenience fields
model_pts <- model_pts %>%
  mutate(
    id      = as.factor(id),
    month   = month(dTime),
    year    = year(dTime),
    yday    = yday(dTime)
  )

# 8. Final checks
cat("\nCombined table dimensions:\n")
print(dim(model_pts))

cat("\nCombined PresAbs counts:\n")
print(table(model_pts$PresAbs))

cat("\nMissing values by column:\n")
print(colSums(is.na(model_pts)))

cat("\nFirst few combined rows:\n")
print(head(model_pts))

# 9. Save output
write.csv(model_pts, out_file, row.names = FALSE)

cat("\nSaved combined file to:\n")
cat(out_file, "\n")



#PART 2: environmental extraction
library(terra)
library(dplyr)
library(lubridate)

# 1 File paths
pts_file <- "C:/github/Whale-SDM/Output/GAMM_prep/254026_pres_pa_combined.csv"

# static rasters
bathy_file    <- "C:/github/Whale-SDM/data/static/bathy_static_8km_5070.tif"
distshelf_file<- "C:/github/Whale-SDM/data/static/dist_shelf_km_8km_5070.tif"
rugosity_file <- "C:/github/Whale-SDM/data/static/rugosity_static_8km_5070.tif"

# dynamic netCDFs
phys_file <- "C:/github/Whale-SDM/data/copernicus/surface_phys/cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-uo-vo-so_130.00W-110.00W_25.00N-45.00N_0.49m_2024-03-01-2024-09-30.nc"
zos_file  <- "C:/github/Whale-SDM/data/copernicus/zos/cmems_mod_glo_phy_my_0.083deg_P1D-m_zos_130.00W-110.00W_25.00N-45.00N_2024-03-01-2024-09-30.nc"

out_dir  <- "C:/github/Whale-SDM/Output/GAMM_env"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, "254026_pres_pa_env_GAMMready.csv")

# 2 Read point table
pts <- read.csv(pts_file, stringsAsFactors = FALSE)

# parse datetime
pts$dTime <- ymd_hms(pts$dTime, tz = "UTC")
pts$date  <- as.Date(pts$dTime)

cat("Point table dimensions:\n")
print(dim(pts))
cat("\nColumns:\n")
print(names(pts))
cat("\nDate range:\n")
print(range(pts$date, na.rm = TRUE))

# 3. Convert points to SpatVector
pts_sv <- vect(pts, geom = c("long", "lat"), crs = "EPSG:4326")

# 4. Read static rasters
bathy    <- rast(bathy_file)
distshelf<- rast(distshelf_file)
rugosity <- rast(rugosity_file)

cat("\nStatic raster summaries:\n")
print(bathy)
print(distshelf)
print(rugosity)

# If static rasters are in projected CRS, reproject points to match each raster

# Bathy
pts_bathy <- project(pts_sv, crs(bathy))
bathy_vals <- terra::extract(bathy, pts_bathy)

# Dist to shelf
pts_distshelf <- project(pts_sv, crs(distshelf))
distshelf_vals <- terra::extract(distshelf, pts_distshelf)

# Rugosity
pts_rugosity <- project(pts_sv, crs(rugosity))
rugosity_vals <- terra::extract(rugosity, pts_rugosity)

# Append static values
pts$bathy_km        <- bathy_vals[, 2]
pts$dist_shelf_km   <- distshelf_vals[, 2]
pts$rugosity        <- rugosity_vals[, 2]

# 5 Read dynamic netCDFs
# These may contain multiple variables and many time layers.
# We inspect names first so you can verify the variable names.
phys <- rast(phys_file)
zos  <- rast(zos_file)

cat("\nSurface physics raster:\n")
print(phys)
cat("\nSurface physics layer names:\n")
print(names(phys))

cat("\nZOS raster:\n")
print(zos)
cat("\nZOS layer names:\n")
print(names(zos))

# Try to read time metadata
phys_time <- terra::time(phys)
zos_time  <- terra::time(zos)

cat("\nFirst few physics times:\n")
print(head(phys_time))
cat("\nFirst few zos times:\n")
print(head(zos_time))

# 6 Helper: find layer indices by variable prefix
get_var_layers <- function(r, prefix) {
  idx <- grep(paste0("^", prefix), names(r))
  if (length(idx) == 0) {
    idx <- grep(prefix, names(r))
  }
  return(idx)
}

thetao_idx <- get_var_layers(phys, "thetao")
uo_idx     <- get_var_layers(phys, "uo")
vo_idx     <- get_var_layers(phys, "vo")
so_idx     <- get_var_layers(phys, "so")
zos_idx    <- get_var_layers(zos,  "zos")

cat("\nLayer counts by variable:\n")
cat("thetao:", length(thetao_idx), "\n")
cat("uo    :", length(uo_idx), "\n")
cat("vo    :", length(vo_idx), "\n")
cat("so    :", length(so_idx), "\n")
cat("zos   :", length(zos_idx), "\n")

# 7. helper: extract one variable by matching date
extract_daily_var <- function(points_df, points_sv, raster_obj, layer_idx, var_name) {
  
  if (length(layer_idx) == 0) {
    warning(paste("No layers found for", var_name))
    points_df[[var_name]] <- NA_real_
    return(points_df)
  }
  
  r_sub <- raster_obj[[layer_idx]]
  r_time <- terra::time(r_sub)
  
  if (is.null(r_time)) {
    stop(paste("No time metadata found for", var_name))
  }
  
  r_dates <- as.Date(r_time)
  
  # Preallocate
  out_vals <- rep(NA_real_, nrow(points_df))
  
  # Loop by unique point date
  udates <- sort(unique(points_df$date))
  
  for (d in udates) {
    row_idx <- which(points_df$date == d)
    
    # nearest matching raster date
    j <- which.min(abs(as.numeric(r_dates - d)))
    
    r_day <- r_sub[[j]]
    
    # extract at those rows only
    pts_day <- points_sv[row_idx]
    vals_day <- terra::extract(r_day, pts_day)
    
    out_vals[row_idx] <- vals_day[, 2]
    
    cat("Extracted", var_name, "for point date", as.character(d),
        "using raster date", as.character(r_dates[j]), "\n")
  }
  
  points_df[[var_name]] <- out_vals
  return(points_df)
}

# 8 Extract dynamic variables
# reprojects points to match each dynamic raster CRS if needed

pts_phys <- project(pts_sv, crs(phys))
pts_zos  <- project(pts_sv, crs(zos))

pts <- extract_daily_var(pts, pts_phys, phys, thetao_idx, "thetao")
pts <- extract_daily_var(pts, pts_phys, phys, uo_idx,     "uo")
pts <- extract_daily_var(pts, pts_phys, phys, vo_idx,     "vo")
pts <- extract_daily_var(pts, pts_phys, phys, so_idx,     "so")
pts <- extract_daily_var(pts, pts_zos,  zos,  zos_idx,    "zos")

# 9. Derive current speed
pts$current_speed <- sqrt((pts$uo^2) + (pts$vo^2))

# 10 Final checks
cat("\nOutput columns:\n")
print(names(pts))

cat("\nMissing values by column:\n")
print(colSums(is.na(pts)))

cat("\nPresAbs table:\n")
print(table(pts$PresAbs))

cat("\nHead of extracted data:\n")
print(head(pts))

summary(pts[, c("bathy_km", "dist_shelf_km", "rugosity",
                "thetao", "uo", "vo", "so", "zos", "current_speed")])


# 11 Save output
write.csv(pts, out_file, row.names = FALSE)

cat("\nSaved extracted environmental table to:\n")
cat(out_file, "\n")