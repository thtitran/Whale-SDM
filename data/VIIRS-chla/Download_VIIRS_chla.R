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

output_dir <- "C:/github/Whale-SDM/data/VIIRS-chla"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

variable <- c(
  "chla"
)

base_url <- paste0(
  "https://coastwatch.pfeg.noaa.gov/erddap/griddap/",
  "erdVHNchla1day.nc"
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
      "VIIRS_chla_",
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

input_dir <- "C:/github/Whale-SDM/data/VIIRS-chla"

input_files <- list.files(
  input_dir,
  pattern = "^VIIRS_chla_[0-9]{8}_[0-9]{8}\\.nc$",
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
  "VIIRS_chla_20240331_20240731_combined.nc"
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
  "chla",
  "units"
)$value

chla_longname <- ncatt_get(
  nc_first,
  "chla",
  "long_name"
)$value

# Inspect original variable dimension order
original_dim_order <- sapply(
  nc_first$var$chla$dim,
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

#Keep date copy
time_manifest <- time_manifest[
  order(time_manifest$time),
]

time_manifest_unique <- time_manifest[
  !duplicated(time_manifest$time),
]

all_time <- time_manifest_unique$time

cat("Combined unique records:", length(all_time), "\n")

time_dim <- ncdim_def(
  name = "time",
  units = time_units,
  vals = all_time,
  unlim = TRUE
)

# ------------------------------------------------------------
# Define dimensions
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

# ------------------------------------------------------------
# Define chlorophyll variable
# ------------------------------------------------------------

chla_var <- ncvar_def(
  name = "chla",
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
  
  # Identify the nonduplicated time slices retained from this file
  retained <- time_manifest_unique[
    time_manifest_unique$file == f,
  ]
  
  if (nrow(retained) == 0) {
    message("No unique dates retained from: ", basename(f))
    next
  }
  
  keep_index <- retained$file_index
  n_time <- length(keep_index)
  
  message(
    "Adding ",
    n_time,
    " unique dates from: ",
    basename(f)
  )
  
  nc_in <- nc_open(f)
  
  chla_i <- ncvar_get(
    nc_in,
    "chla",
    collapse_degen = FALSE
  )
  
  nc_close(nc_in)
  
  # Retain only time slices not already supplied by another file
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
    varid = "chla",
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

# Confirm chlorophyll contains valid values
chla_check <- ncvar_get(
  nc,
  "chla"
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