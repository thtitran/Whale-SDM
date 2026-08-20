library(httr)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

start_date <- as.Date("2024-03-31")
end_date   <- as.Date("2024-07-31")

lat_min <- 33
lat_max <- 40

lon_min <- -126
lon_max <- -119

output_dir <- "C:/github/Whale-SDM/data/VIIRSNOAA20GAPLESS"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

variable <- c(
  "chlor_a"
)

base_url <- paste0(
  "https://coastwatch.pfeg.noaa.gov/erddap/griddap/",
  "nesdisVHNnoaaSNPPnoaa20chlaGapfilledDaily.nc"
)

# ------------------------------------------------------------
# Create monthly date intervals
# ------------------------------------------------------------

month_starts <- seq(
  from = as.Date(format(start_date, "%Y-%m-01")),
  to   = as.Date(format(end_date, "%Y-%m-01")),
  by   = "month"
)

date_intervals <- data.frame(
  interval_start = pmax(month_starts, start_date),
  interval_end = pmin(
    c(month_starts[-1] - 1, end_date),
    end_date
  )
)

print(date_intervals)

# ------------------------------------------------------------
# Download each monthly interval
# ------------------------------------------------------------

for (i in seq_len(nrow(date_intervals))) {
  
  interval_start <- date_intervals$interval_start[i]
  interval_end   <- date_intervals$interval_end[i]
  
  # Dimension order:
  # time, altitude, latitude, longitude
  #
  # Latitude is stored north-to-south, so use lat_max first.
  dimension_constraint <- paste0(
    "[(", interval_start, "T12:00:00Z):1:(",
    interval_end,   "T12:00:00Z)]",
    "[(0.0)]",
    "[(", lat_max, "):1:(", lat_min, ")]",
    "[(", lon_min, "):1:(", lon_max, ")]"
  )
  
  constraint <- paste0(
    variable,
    dimension_constraint
  )
  
  download_url <- paste0(
    base_url,
    "?",
    constraint
  )
  
  output_file <- file.path(
    output_dir,
    paste0(
      "NOAA20_VIIRS_GL_chla_",
      format(interval_start, "%Y%m%d"),
      "_",
      format(interval_end, "%Y%m%d"),
      ".nc"
    )
  )
  
  message(
    "\nDownloading ",
    interval_start,
    " through ",
    interval_end
  )
  
  cat(download_url, "\n")
  
  if (file.exists(output_file)) {
    unlink(output_file)
  }
  
  response <- RETRY(
    verb = "GET",
    url = download_url,
    write_disk(output_file, overwrite = TRUE),
    progress(),
    timeout(1800),
    times = 5,
    pause_base = 5,
    pause_cap = 60,
    terminate_on = c(400, 401, 403, 404)
  )
  
  status <- status_code(response)
  
  if (status != 200) {
    
    if (file.exists(output_file)) {
      unlink(output_file)
    }
    
    warning(
      "Download failed for ",
      interval_start,
      " through ",
      interval_end,
      "; HTTP status ",
      status
    )
    
    next
  }
  
  file_size <- file.info(output_file)$size
  
  # Guard against small ERDDAP error-response files
  if (is.na(file_size) || file_size < 1000) {
    unlink(output_file)
    
    warning(
      "Downloaded file was unexpectedly small for ",
      interval_start,
      " through ",
      interval_end
    )
    
    next
  }
  
  message(
    "Saved: ",
    output_file,
    "\nSize: ",
    round(file_size / 1024^2, 2),
    " MB"
  )
}

#Combining netCDF files

library(ncdf4)

# ------------------------------------------------------------
# Files
# ------------------------------------------------------------

input_dir <- "C:/github/Whale-SDM/data/VIIRSNOAA20GAPLESS"

input_files <- list.files(
  input_dir,
  pattern = "^NOAA20_VIIRS_GL_chla_[0-9]{8}_[0-9]{8}\\.nc$",
  full.names = TRUE
)

input_files <- sort(input_files)

input_files <- input_files[
  !grepl("_combined\\.nc$", input_files)
]

if (length(input_files) == 0) {
  stop("No VIIRS chlorophyll NetCDF files were found.")
}

output_file <- file.path(
  input_dir,
  "NOAA20_VIIRS_GL_chla_20240331_20240731_combined.nc"
)

print(input_files)

# ------------------------------------------------------------
# Read coordinates and metadata from first file
# ------------------------------------------------------------

nc_first <- nc_open(input_files[1])

lon <- ncvar_get(nc_first, "longitude")
lat <- ncvar_get(nc_first, "latitude")
alt <- ncvar_get(nc_first, "altitude")

lon_units <- ncatt_get(
  nc_first,
  "longitude",
  "units"
)$value

lat_units <- ncatt_get(
  nc_first,
  "latitude",
  "units"
)$value

alt_units <- ncatt_get(
  nc_first,
  "altitude",
  "units"
)$value

time_units <- ncatt_get(
  nc_first,
  "time",
  "units"
)$value

chla_units <- ncatt_get(
  nc_first,
  "chlor_a",
  "units"
)$value

chla_longname <- ncatt_get(
  nc_first,
  "chlor_a",
  "long_name"
)$value

original_dim_order <- sapply(
  nc_first$var$chlor_a$dim,
  function(x) x$name
)

print(original_dim_order)

nc_close(nc_first)

# ------------------------------------------------------------
# Build a manifest of every time record in every file
# ------------------------------------------------------------

time_manifest <- do.call(
  rbind,
  lapply(input_files, function(f) {
    
    nc <- nc_open(f)
    
    time_i <- ncvar_get(nc, "time")
    
    nc_close(nc)
    
    data.frame(
      file = f,
      file_index = seq_along(time_i),
      time = time_i,
      stringsAsFactors = FALSE
    )
  })
)

time_manifest$date <- as.POSIXct(
  time_manifest$time,
  origin = "1970-01-01",
  tz = "UTC"
)

cat("Total records:", nrow(time_manifest), "\n")
cat("Unique records:", length(unique(time_manifest$time)), "\n")

#Display all duplicated dates:

duplicated_times <- time_manifest[
  duplicated(time_manifest$time) |
    duplicated(time_manifest$time, fromLast = TRUE),
]

duplicated_times[
  order(duplicated_times$time),
  c("date", "file", "file_index")
]

# Sort
time_manifest <- time_manifest[
  order(time_manifest$time),
]

# Deduplicate
time_manifest_unique <- time_manifest[
  !duplicated(time_manifest$time),
]

# Remove dates outside March 31–July 31
time_manifest_unique <- time_manifest_unique[
  as.Date(time_manifest_unique$date) >= start_date &
    as.Date(time_manifest_unique$date) <= end_date,
]

# Define the output time axis
all_time <- time_manifest_unique$time

cat("Combined unique records:", length(all_time), "\n")
print(range(time_manifest_unique$date))

# Define dimensions
lon_dim <- ncdim_def(
  name = "longitude",
  units = lon_units,
  vals = lon
)

lat_dim <- ncdim_def(
  name = "latitude",
  units = lat_units,
  vals = lat
)

alt_dim <- ncdim_def(
  name = "altitude",
  units = alt_units,
  vals = alt
)

time_dim <- ncdim_def(
  name = "time",
  units = time_units,
  vals = all_time,
  unlim = TRUE
)
# ------------------------------------------------------------
# Define chlorophyll variable
# ------------------------------------------------------------
chla_var <- ncvar_def(
  name = "chlor_a",
  units = chla_units,
  dim = list(
    lon_dim,
    lat_dim,
    alt_dim,
    time_dim
  ),
  missval = NA_real_,
  longname = chla_longname,
  prec = "float",
  compression = 4
)

# ------------------------------------------------------------
# Create combined NetCDF
# ------------------------------------------------------------

if (file.exists(output_file)) {
  unlink(output_file)
}

nc_out <- nc_create(
  filename = output_file,
  vars = list(chla_var),
  force_v4 = TRUE
)

# ------------------------------------------------------------
# Append each source file
# ------------------------------------------------------------

time_start <- 1

for (f in input_files) {
  
  retained <- time_manifest_unique[
    time_manifest_unique$file == f,
  ]
  
  if (nrow(retained) == 0) {
    message("No retained dates from: ", basename(f))
    next
  }
  
  keep_index <- retained$file_index
  n_time <- length(keep_index)
  
  message(
    "Adding ",
    n_time,
    " retained dates from: ",
    basename(f)
  )
  
  nc_in <- nc_open(f)
  
  chla_i <- ncvar_get(
    nc_in,
    "chlor_a",
    collapse_degen = FALSE
  )
  
  nc_close(nc_in)
  
  chla_i <- chla_i[
    ,
    ,
    ,
    keep_index,
    drop = FALSE
  ]
  
  expected_dim <- c(
    length(lon),
    length(lat),
    length(alt),
    n_time
  )
  
  if (!identical(dim(chla_i), expected_dim)) {
    stop(
      "Unexpected dimensions in ",
      basename(f),
      "\nFound: ",
      paste(dim(chla_i), collapse = " x "),
      "\nExpected: ",
      paste(expected_dim, collapse = " x ")
    )
  }
  
  ncvar_put(
    nc_out,
    varid = "chlor_a",
    vals = chla_i,
    start = c(
      1,
      1,
      1,
      time_start
    ),
    count = c(
      length(lon),
      length(lat),
      length(alt),
      n_time
    )
  )
  
  time_start <- time_start + n_time
}

nc_close(nc_out)


# ------------------------------------------------------------
# Verify combined NetCDF
# ------------------------------------------------------------

nc <- nc_open(output_file)

names(nc$var)
names(nc$dim)

sapply(
  nc$var,
  function(x) x$varsize
)

time_vals <- ncvar_get(
  nc,
  "time"
)

time_units_check <- ncatt_get(
  nc,
  "time",
  "units"
)$value

dates <- as.POSIXct(
  time_vals,
  origin = "1970-01-01",
  tz = "UTC"
)

cat("Time units:", time_units_check, "\n")
cat("Number of dates:", length(dates), "\n")
cat("Unique dates:", length(unique(dates)), "\n")
print(range(dates))

chla_check <- ncvar_get(
  nc,
  "chlor_a"
)

cat(
  "Chlorophyll range:",
  range(chla_check, na.rm = TRUE),
  "\n"
)

cat(
  "Missing-value proportion:",
  mean(is.na(chla_check)),
  "\n"
)

nc_close(nc)



# ------------------------------------------------------------
# Read data while retaining the altitude dimension
# ------------------------------------------------------------

nc <- nc_open(output_file)

chla <- ncvar_get(
  nc,
  "chlor_a",
  collapse_degen = FALSE
)

time_vals <- ncvar_get(nc, "time")
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")

nc_close(nc)

available_dates <- as.Date(
  as.POSIXct(
    time_vals,
    origin = "1970-01-01",
    tz = "UTC"
  )
)

dim(chla)
# Expected: 85 85 1 97


# ------------------------------------------------------------
# 1. Identify dates absent from the NetCDF
# ------------------------------------------------------------

expected_dates <- seq(
  from = start_date,
  to   = end_date,
  by   = "day"
)

absent_dates <- expected_dates[
  !expected_dates %in% available_dates
]

cat("Expected dates:", length(expected_dates), "\n")
cat("Dates stored:", length(available_dates), "\n")
cat("Dates completely absent:", length(absent_dates), "\n")

print(absent_dates)


# ------------------------------------------------------------
# 2. Calculate missingness for each stored date
# ------------------------------------------------------------

chla_dim <- dim(chla)

n_spatial_cells <- prod(chla_dim[1:3])
n_stored_dates  <- chla_dim[4]

# Rows = spatial cells
# Columns = stored dates

na_matrix <- matrix(
  is.na(chla),
  nrow = n_spatial_cells,
  ncol = n_stored_dates
)

date_missing_count <- colSums(na_matrix)
date_missing_prop  <- colMeans(na_matrix)

date_summary <- data.frame(
  date = available_dates,
  missing_cells = date_missing_count,
  total_cells = n_spatial_cells,
  missing_proportion = date_missing_prop,
  entire_layer_missing = date_missing_prop == 1
)

print(date_summary)

cat(
  "Stored dates with an entirely missing layer:",
  sum(date_summary$entire_layer_missing),
  "\n"
)

summary(date_summary$missing_proportion)


# ------------------------------------------------------------
# Plot missingness through time
# ------------------------------------------------------------

plot(
  date_summary$date,
  date_summary$missing_proportion,
  type = "b",
  pch = 16,
  xlim = range(expected_dates),
  ylim = c(0, 1),
  xlab = "Date",
  ylab = "Proportion of grid cells missing",
  main = "Chlorophyll missingness through time"
)

# Completely absent dates shown at 100% missing
points(
  absent_dates,
  rep(1, length(absent_dates)),
  pch = 4,
  col = "red",
  lwd = 2
)

abline(
  h = mean(is.na(chla)),
  col = "blue",
  lty = 2
)

legend(
  "bottomleft",
  legend = c(
    "Stored date",
    "Date absent",
    "Mean within-date missingness"
  ),
  col = c("black", "red", "blue"),
  pch = c(16, 4, NA),
  lty = c(1, NA, 2),
  bty = "n"
)


# ------------------------------------------------------------
# 3. Calculate missingness for each grid cell
# ------------------------------------------------------------

# Exclude stored dates whose complete layer is NA
usable_dates <- !date_summary$entire_layer_missing

cell_missing_prop <- rowMeans(
  na_matrix[, usable_dates, drop = FALSE]
)

cat(
  "Cells never missing:",
  sum(cell_missing_prop == 0),
  "\n"
)

cat(
  "Cells intermittently missing:",
  sum(cell_missing_prop > 0 & cell_missing_prop < 1),
  "\n"
)

cat(
  "Cells always missing:",
  sum(cell_missing_prop == 1),
  "\n"
)


# Convert back to longitude × latitude matrix
# This works because altitude has length 1

cell_missing_map <- matrix(
  cell_missing_prop,
  nrow = length(lon),
  ncol = length(lat)
)


# ------------------------------------------------------------
# Reorder coordinates for image()
# ------------------------------------------------------------

# order() makes both coordinate vectors increase
lon_order <- order(lon)
lat_order <- order(lat)

lon_plot <- lon[lon_order]
lat_plot <- lat[lat_order]

cell_missing_map_plot <- cell_missing_map[
  lon_order,
  lat_order,
  drop = FALSE
]

stopifnot(
  all(diff(lon_plot) > 0),
  all(diff(lat_plot) > 0)
)


# ------------------------------------------------------------
# Plot spatial missingness
# ------------------------------------------------------------

image(
  x = lon_plot,
  y = lat_plot,
  z = cell_missing_map_plot,
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Proportion of available dates missing per cell",
  col = hcl.colors(50, "YlOrRd"),
  zlim = c(0, 1),
  useRaster = FALSE
)


# ------------------------------------------------------------
# 4. Partition missing cell-days over the requested period
# ------------------------------------------------------------

n_absent_dates <- length(absent_dates)

n_all_na_dates <- sum(
  date_summary$entire_layer_missing
)

# Missing because an entire calendar date is unavailable
temporal_missing_cell_days <-
  (n_absent_dates + n_all_na_dates) *
  n_spatial_cells

# Missing cells within otherwise usable dates
within_date_missing_cell_days <- sum(
  date_missing_count[
    !date_summary$entire_layer_missing
  ]
)

total_expected_cell_days <-
  length(expected_dates) *
  n_spatial_cells

missingness_partition <- data.frame(
  missingness_type = c(
    "Missing dates / entirely missing layers",
    "Missing cells within usable dates",
    "All missingness combined"
  ),
  missing_cell_days = c(
    temporal_missing_cell_days,
    within_date_missing_cell_days,
    temporal_missing_cell_days +
      within_date_missing_cell_days
  ),
  proportion_of_expected_dataset = c(
    temporal_missing_cell_days,
    within_date_missing_cell_days,
    temporal_missing_cell_days +
      within_date_missing_cell_days
  ) / total_expected_cell_days
)

print(missingness_partition)