suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(sf)
  library(geosphere)
  library(readr)
})

mpm_path <- "C:/github/Whale-SDM/tracks/outputnewgamma/combined_mpm.csv"
out_dir  <- "C:/github/Whale-SDM/tracks/pseudo_absences"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#PROJECTION FOR CALIFORNIA
sim_crs <- 3310  # NAD83 / California Albers

#shapefile
land_shp <-
  st_read("C:/path/to/your/California_or_global_coastline.shp") %>%
  st_transform(4326)

#reads gamma classified tracks
mpm_all <- read_csv(mpm_path) %>%
  mutate(
    date = as.POSIXct(date, tz = "UTC"),
    id   = as.character(id)
  )

#  sanity check
table(mpm_all$behav)
#levels should include "low_persist" and "high_persist"

#simulation func

sim_fun <- function(track, land_shp, projCRS = sim_crs,
                    mult = 1, check_land = TRUE, quant = 1) {
  
  # land in projected CRS
  land_utm <- st_transform(land_shp, projCRS)
  
  # compute movement metrics on observed track
  track_sf <- track %>%
    mutate(
      turn_angle = bearing(cbind(lon, lat),
                           cbind(lead(lon), lead(lat))),       # degrees
      step_dist_m = distVincentySphere(cbind(lag(lon), lag(lat)),
                                       cbind(lon, lat))        # metres
    ) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    st_transform(projCRS)
  
  coords <- st_coordinates(track_sf)
  
  track_sum <- track_sf %>%
    mutate(
      x = coords[, 1],
      y = coords[, 2],
      on_land = as.integer(st_intersects(., land_utm))
    ) %>%
    st_drop_geometry() %>%
    filter(is.na(on_land))   # keep only water points
  
  if (nrow(track_sum) < 2) return(NULL)
  
  turning_angle <- track_sum$turn_angle
  step_distance <- track_sum$step_dist_m
  n_steps       <- nrow(track_sum)
  
  # trim extremes (optional)
  turning_angle <- subset(
    turning_angle,
    turning_angle <= quantile(turning_angle, 1 - (quant/2), na.rm = TRUE) &
      turning_angle >= quantile(turning_angle, quant/2,     na.rm = TRUE)
  )
  step_distance <- subset(
    step_distance,
    step_distance <= quantile(step_distance, 1 - quant, na.rm = TRUE)
  )
  
  if (length(turning_angle) == 0 || length(step_distance) == 0) return(NULL)
  
  pb <- txtProgressBar(max = n_steps, style = 3)
  
  # start at first observed point (projected)
  sim_track <- track_sum %>%
    arrange(date) %>%
    slice(1) %>%
    transmute(step = 0L, x, y)
  
  for (s in 1:n_steps) {
    repeat {
      heading_deg <- sample(na.omit(turning_angle), 1)
      step_m      <- sample(na.omit(step_distance), 1)
      
      theta <- heading_deg * pi/180  # degrees → radians
      
      new_point <- tibble(
        step = s,
        x    = (cos(theta) * step_m * mult) + dplyr::last(sim_track$x),
        y    = (sin(theta) * step_m * mult) + dplyr::last(sim_track$y)
      ) %>%
        st_as_sf(coords = c("x", "y"), crs = projCRS, remove = FALSE) %>%
        mutate(on_land = as.integer(st_intersects(., land_utm))) %>%
        st_drop_geometry()
      
      if (!check_land) break
      if (is.na(new_point$on_land)) break  # water = accept, break repeat{}
    }
    
    sim_track <- bind_rows(sim_track, select(new_point, -on_land))
    setTxtProgressBar(pb, s)
  }
  close(pb)
  
  # align with original id/date/behav
  sim_out <- track_sum %>%
    transmute(id, date, behav) %>%
    bind_cols(sim_track %>% slice(1:nrow(track_sum)))
  
  sim_out
}

rep_fun <- function(track, land_shp, n_rep = 5, ...) {
  map_dfr(seq_len(n_rep), function(i) {
    message("   rep ", i)
    sim_fun(track, land_shp, ...) %>% mutate(rep = i)
  })
}

#simbulation by behavior and id

## low-persistence (foraging-ish) pseudo-absences
low_dat <- mpm_all %>%
  filter(behav == "low_persist") %>%
  arrange(id, date)

low_by_id <- split(low_dat, low_dat$id)

sim_low <- map_dfr(low_by_id, ~ {
  message("\nSimulating LOW for ID = ", .x$id[1])
  rep_fun(.x, land_shp, n_rep = 5, mult = 1, quant = 0.2, check_land = TRUE)
})

write_csv(sim_low, file.path(out_dir, "pseudo_low_persist_tracks.csv"))

## high-persistence (migration-ish) pseudo-absences
high_dat <- mpm_all %>%
  filter(behav == "high_persist") %>%
  arrange(id, date)

high_by_id <- split(high_dat, high_dat$id)

sim_high <- map_dfr(high_by_id, ~ {
  message("\nSimulating HIGH for ID = ", .x$id[1])
  rep_fun(.x, land_shp, n_rep = 5, mult = 1, quant = 0.2, check_land = TRUE)
})

write_csv(sim_high, file.path(out_dir, "pseudo_high_persist_tracks.csv"))
