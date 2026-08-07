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

output_dir <- "C:/github/Whale-SDM/data/L3S-LEO"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

variables <- c(
  "sea_surface_temperature",
  "sst_gradient_magnitude",
  "sst_front_position"
)

base_url <- paste0(
  "https://coastwatch.noaa.gov/erddap/griddap/",
  "noaacwLEOACSPOSSTL3SCDaily.nc"
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
  
  dimension_constraint <- paste0(
    "[(", interval_start, "T12:00:00Z):1:(",
    interval_end,   "T12:00:00Z)]",
    "[(", lat_min, "):1:(", lat_max, ")]",
    "[(", lon_min, "):1:(", lon_max, ")]"
  )
  
  constraint <- paste0(
    variables,
    dimension_constraint,
    collapse = ","
  )
  
  download_url <- paste0(
    base_url,
    "?",
    constraint
  )
  
  output_file <- file.path(
    output_dir,
    paste0(
      "L3S_LEO_",
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
  
  # Remove an earlier incomplete/error file
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
  
  file_size_mb <- file.info(output_file)$size / 1024^2
  
  message(
    "Saved: ",
    output_file,
    "\nSize: ",
    round(file_size_mb, 2),
    " MB"
  )
}



#Combining netCDF files

library(ncdf4)

# ------------------------------------------------------------
# Files
# ------------------------------------------------------------

input_dir <- "C:/github/Whale-SDM/data/L3S-LEO"

input_files <- list.files(
  input_dir,
  pattern = "^L3S_LEO_[0-9]{8}_[0-9]{8}\\.nc$",
  full.names = TRUE
)

input_files <- sort(input_files)

# Make sure the previous incorrect combined file is excluded
input_files <- input_files[
  !grepl("_combined\\.nc$", input_files)
]

output_file <- file.path(
  input_dir,
  "L3S_LEO_20240331_20240731_combined.nc"
)

input_files

# Read longitude and latitude from the first file
nc_first <- nc_open(input_files[1])

lon <- ncvar_get(nc_first, "longitude")
lat <- ncvar_get(nc_first, "latitude")

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

nc_close(nc_first)

# Read all time values
all_time <- numeric()

for (f in input_files) {
  
  nc <- nc_open(f)
  
  time_i <- ncvar_get(nc, "time")
  
  all_time <- c(
    all_time,
    time_i
  )
  
  nc_close(nc)
}

length(all_time)

time_units_all <- character(length(input_files))

for (i in seq_along(input_files)) {
  
  nc <- nc_open(input_files[i])
  
  time_units_all[i] <- ncatt_get(
    nc,
    "time",
    "units"
  )$value
  
  nc_close(nc)
}

unique(time_units_all)

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

time_dim <- ncdim_def(
  name = "time",
  units = time_units,
  vals = all_time,
  unlim = TRUE
)

nc_test <- nc_open(input_files[1])

nc_test$var$sea_surface_temperature$dimids
sapply(
  nc_test$var$sea_surface_temperature$dim,
  function(x) x$name
)

nc_close(nc_test)

sst_var <- ncvar_def(
  name = "sea_surface_temperature",
  units = "degree_C",
  dim = list(
    lon_dim,
    lat_dim,
    time_dim
  ),
  missval = NA_real_,
  longname = "Sea surface temperature",
  prec = "float",
  compression = 4
)

gradient_var <- ncvar_def(
  name = "sst_gradient_magnitude",
  units = "degree_C km-1",
  dim = list(
    lon_dim,
    lat_dim,
    time_dim
  ),
  missval = NA_real_,
  longname = "Sea surface temperature gradient magnitude",
  prec = "float",
  compression = 4
)

front_var <- ncvar_def(
  name = "sst_front_position",
  units = "1",
  dim = list(
    lon_dim,
    lat_dim,
    time_dim
  ),
  missval = NA_real_,
  longname = "Sea surface temperature front position",
  prec = "float",
  compression = 4
)

#Inspect original units
nc_test <- nc_open(input_files[1])

ncatt_get(
  nc_test,
  "sea_surface_temperature",
  "units"
)$value

ncatt_get(
  nc_test,
  "sst_gradient_magnitude",
  "units"
)$value

ncatt_get(
  nc_test,
  "sst_front_position",
  "units"
)$value

nc_close(nc_test)


if (file.exists(output_file)) {
  unlink(output_file)
}

nc_out <- nc_create(
  filename = output_file,
  vars = list(
    sst_var,
    gradient_var,
    front_var
  ),
  force_v4 = TRUE
)

time_start <- 1

for (f in input_files) {
  
  message(
    "Adding: ",
    basename(f)
  )
  
  nc_in <- nc_open(f)
  
  time_i <- ncvar_get(
    nc_in,
    "time"
  )
  
  n_time <- length(time_i)
  
  sst_i <- ncvar_get(
    nc_in,
    "sea_surface_temperature"
  )
  
  gradient_i <- ncvar_get(
    nc_in,
    "sst_gradient_magnitude"
  )
  
  front_i <- ncvar_get(
    nc_in,
    "sst_front_position"
  )
  
  # Ensure one-day files retain an explicit time dimension
  if (length(dim(sst_i)) == 2) {
    dim(sst_i) <- c(
      length(lon),
      length(lat),
      1
    )
  }
  
  if (length(dim(gradient_i)) == 2) {
    dim(gradient_i) <- c(
      length(lon),
      length(lat),
      1
    )
  }
  
  if (length(dim(front_i)) == 2) {
    dim(front_i) <- c(
      length(lon),
      length(lat),
      1
    )
  }
  
  ncvar_put(
    nc_out,
    "sea_surface_temperature",
    sst_i,
    start = c(
      1,
      1,
      time_start
    ),
    count = c(
      length(lon),
      length(lat),
      n_time
    )
  )
  
  ncvar_put(
    nc_out,
    "sst_gradient_magnitude",
    gradient_i,
    start = c(
      1,
      1,
      time_start
    ),
    count = c(
      length(lon),
      length(lat),
      n_time
    )
  )
  
  ncvar_put(
    nc_out,
    "sst_front_position",
    front_i,
    start = c(
      1,
      1,
      time_start
    ),
    count = c(
      length(lon),
      length(lat),
      n_time
    )
  )
  
  time_start <- time_start + n_time
  
  nc_close(nc_in)
}

nc_close(nc_out)


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

time_units <- ncatt_get(
  nc,
  "time",
  "units"
)$value

dates <- as.POSIXct(
  time_vals,
  origin = "1970-01-01",
  tz = "UTC"
)

length(dates)
length(unique(dates))
range(dates)

nc_close(nc)