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

output_dir <- "C:/github/Whale-SDM/data/VIIRS-kd490GAPLESS"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

variable <- "kd_490"

base_url <- paste0(
  "https://coastwatch.pfeg.noaa.gov/erddap/griddap/",
  "nesdisNPPN20S3AkdSCIDINEOFDaily.nc"
)

# ------------------------------------------------------------
# Create monthly intervals
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
      "VIIRS_Kd490_GL_",
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

input_dir <- "C:/github/Whale-SDM/data/VIIRS-kd490GAPLESS"

input_files <- list.files(
  input_dir,
  pattern = "^VIIRS_Kd490_GL_[0-9]{8}_[0-9]{8}\\.nc$",
  full.names = TRUE
)

input_files <- sort(input_files)

input_files <- input_files[
  !grepl("_combined\\.nc$", input_files)
]

if (length(input_files) == 0) {
  stop("No VIIRS Kd490 NetCDF files were found.")
}

output_file <- file.path(
  input_dir,
  "VIIRS_Kd490_GL_20240331_20240731_combined.nc"
)

print(input_files)

# ------------------------------------------------------------
# Read coordinates and metadata
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

kd_units <- ncatt_get(
  nc_first,
  "kd_490",
  "units"
)$value

kd_longname <- ncatt_get(
  nc_first,
  "kd_490",
  "long_name"
)$value

original_dim_order <- sapply(
  nc_first$var$kd_490$dim,
  function(x) x$name
)

print(original_dim_order)

nc_close(nc_first)

# ------------------------------------------------------------
# Build time manifest
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

# Keep requested dates only
time_manifest <- time_manifest[
  as.Date(time_manifest$date) >= start_date &
    as.Date(time_manifest$date) <= end_date,
]

# Sort and retain one copy of each timestamp
time_manifest <- time_manifest[
  order(time_manifest$time),
]

time_manifest_unique <- time_manifest[
  !duplicated(time_manifest$time),
]

all_time <- time_manifest_unique$time

cat("Total source records:", nrow(time_manifest), "\n")
cat("Unique retained records:", length(all_time), "\n")

# ------------------------------------------------------------
# Validate grids and time units
# ------------------------------------------------------------

for (f in input_files) {
  
  nc <- nc_open(f)
  
  lon_i <- ncvar_get(nc, "longitude")
  lat_i <- ncvar_get(nc, "latitude")
  alt_i <- ncvar_get(nc, "altitude")
  
  time_units_i <- ncatt_get(
    nc,
    "time",
    "units"
  )$value
  
  nc_close(nc)
  
  if (!isTRUE(all.equal(lon_i, lon))) {
    stop("Longitude grid differs in: ", basename(f))
  }
  
  if (!isTRUE(all.equal(lat_i, lat))) {
    stop("Latitude grid differs in: ", basename(f))
  }
  
  if (!isTRUE(all.equal(alt_i, alt))) {
    stop("Altitude differs in: ", basename(f))
  }
  
  if (!identical(time_units_i, time_units)) {
    stop("Time units differ in: ", basename(f))
  }
}

# ------------------------------------------------------------
# Define output NetCDF
# ------------------------------------------------------------

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

kd490_var <- ncvar_def(
  name = "kd_490",
  units = kd_units,
  dim = list(
    lon_dim,
    lat_dim,
    alt_dim,
    time_dim
  ),
  missval = NA_real_,
  longname = kd_longname,
  prec = "float",
  compression = 4
)

if (file.exists(output_file)) {
  unlink(output_file)
}

nc_out <- nc_create(
  filename = output_file,
  vars = list(kd490_var),
  force_v4 = TRUE
)

# ------------------------------------------------------------
# Append files
# ------------------------------------------------------------

time_start <- 1

for (f in input_files) {
  
  retained <- time_manifest_unique[
    time_manifest_unique$file == f,
  ]
  
  if (nrow(retained) == 0) {
    message("No unique records retained from: ", basename(f))
    next
  }
  
  keep_index <- retained$file_index
  n_time <- length(keep_index)
  
  message(
    "Adding ",
    n_time,
    " records from: ",
    basename(f)
  )
  
  nc_in <- nc_open(f)
  
  kd_i <- ncvar_get(
    nc_in,
    "kd_490",
    collapse_degen = FALSE
  )
  
  nc_close(nc_in)
  
  kd_i <- kd_i[
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
  
  if (!identical(dim(kd_i), expected_dim)) {
    stop(
      "Unexpected dimensions in ",
      basename(f),
      "\nFound: ",
      paste(dim(kd_i), collapse = " x "),
      "\nExpected: ",
      paste(expected_dim, collapse = " x ")
    )
  }
  
  ncvar_put(
    nc_out,
    varid = "kd_490",
    vals = kd_i,
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
# Verify
# ------------------------------------------------------------

nc <- nc_open(output_file)

names(nc$var)
names(nc$dim)

sapply(
  nc$var,
  function(x) x$varsize
)

time_vals <- ncvar_get(nc, "time")

dates <- as.POSIXct(
  time_vals,
  origin = "1970-01-01",
  tz = "UTC"
)

kd_check <- ncvar_get(
  nc,
  "kd_490"
)

cat("Number of dates:", length(dates), "\n")
cat("Unique dates:", length(unique(dates)), "\n")
print(range(dates))

cat(
  "Kd490 range:",
  range(kd_check, na.rm = TRUE),
  "\n"
)

cat(
  "Missing proportion:",
  mean(is.na(kd_check)),
  "\n"
)

nc_close(nc)