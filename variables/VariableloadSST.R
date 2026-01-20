##
# SST download + weekly/monthly rolling means 
# Example: (Spring 2024)
# California Current bbox: -130 to -110 lon, 25 to 45 lat
# Data source: NOAA NCEI ERDDAP OISST (0.25°, daily)
##

library(terra)  
library(rerddap)  
library(zoo)

#Settings
bbox_ll <- c(w = -130, e = -110, s = 25, n = 45)          # lon/lat bbox in -180..180
date_rng <- c(start = "2024-03-01", end = "2024-05-31")   # Spring 2024
erddap_url <- "https://www.ncei.noaa.gov/erddap/"
dataset_id <- "ncdc_oisst_v2_avhrr_by_time_zlev_lat_lon"
out_nc <- "oisst_spring2024_raw.nc"                       # final file name we want

# Convert bbox longitudes to the dataset lon convention (0..360)
# ERDDAP ocean products use 0-360. Track bbox is -130..-110, so add 360.
bbox_0360_lon <- unname(bbox_ll[c("w", "e")] + 360)

# ---- 3) Download SST subset as netCDF ----
# griddap() asks ERDDAP to crop on the SERVER (fast + avoids huge downloads).
# This dataset includes a "depth" dimension even though it's surface-only;
#  fix via having 0 by giving depth = c(0, 0).
# NOTE: store=disk(path=".") saves to a hashed filename; we rename it after.
res <- griddap(
  datasetx  = dataset_id,
  url       = erddap_url,
  time      = unname(date_rng),
  depth     = c(0, 0),
  latitude  = unname(bbox_ll[c("s", "n")]),
  longitude = bbox_0360_lon,
  fields    = "sst",
  fmt       = "nc",
  store     = disk(path = ".")
)

# Rename the hashed file to something stable and memorable.
# (If the destination exists, remove it first to avoid rename failures.)
if (file.exists(out_nc)) file.remove(out_nc)
file.rename(from = res$summary$filename, to = out_nc)

# ---- 4) Read the netCDF as a SpatRaster and shift longitudes back to -180..180 ----
# The downloaded file is in 0..360 lon; shifting by -360 puts it back into
# the common -180..180 system that matches most sf track data and maps.
sst_raw   <- rast(out_nc)
sst_fixed <- shift(sst_raw, dx = -360)

# one test plot
plot(sst_fixed[[1]], main = "Raw daily SST (°C) — Spring 2024 (first day)")

#  Build rolling means over time (7 day and 30 day trailing windows)
# app() applies a function to each grid cell’s full time series.
# rollapply() returns a vector the same length as the time series, where:
# - width=7 means "last 7 days"
# - align='right' means the value at day t uses days (t-6 ... t)
# - fill=NA means early days without enough history become NA
sst_7d <- app(sst_fixed, fun = function(v) {
  zoo::rollapply(v, width = 7, FUN = mean, align = "right", fill = NA, na.rm = TRUE)
})

sst_30d <- app(sst_fixed, fun = function(v) {
  zoo::rollapply(v, width = 30, FUN = mean, align = "right", fill = NA, na.rm = TRUE)
})

# Keep timestamps consistent (terra has its own time() methods)
terra::time(sst_7d)  <- terra::time(sst_fixed)
terra::time(sst_30d) <- terra::time(sst_fixed)

# Friendly layer names (optional, but helps when debugging)
names(sst_7d)  <- paste0("sst_7d_",  seq_len(nlyr(sst_7d)))
names(sst_30d) <- paste0("sst_30d_", seq_len(nlyr(sst_30d)))

# visual check (raw vs 7d vs 30d on the same day) ----
i <- 40  # pick a mid-spring day
plot(sst_fixed[[i]], main = paste("Raw SST (°C)", terra::time(sst_fixed)[i]))
plot(sst_7d[[i]],    main = paste("7-day mean SST (°C)", terra::time(sst_7d)[i]))
plot(sst_30d[[i]],   main = paste("30-day mean SST (°C)", terra::time(sst_30d)[i]))

# Save the rolling-mean rasters so you don’t recompute
# writeRaster(sst_7d,  "sst_7d_spring2024.tif",  overwrite = TRUE)
# writeRaster(sst_30d, "sst_30d_spring2024.tif", overwrite = TRUE)
