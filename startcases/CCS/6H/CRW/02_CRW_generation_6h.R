library(adehabitatLT)
library(maps)
library(mapdata)
library(sf)
library(dplyr)

setwd("C:/github/Whale-SDM/startcases/CCS/6H/CRW")

track_dir <- "C:/github/Whale-SDM/processed tracks"

track_files <- c(
  file.path(track_dir, "254025-RawArgos_pred_6h.csv"),
  file.path(track_dir, "254026-RawArgos_pred_6h.csv"),
  file.path(track_dir, "254027-RawArgos_pred_6h.csv"),
  file.path(track_dir, "254030-RawArgos_pred_6h.csv")
)


out_dir <- paste0(
  "C:/github/Whale-SDM/startcases/CCS/6H/CRW")
  
#tracks per whale
n_crw <- 200

#limits used for land/boundary rejection
xlim_land <- c(-150, -90)
ylim_land <- c(0, 55)


#func to read and standardize one observed track
read_whale_track <- function(file) {
  
  dat <- read.csv(
    file,
    stringsAsFactors = FALSE
  )
  
  required_columns <- c(
    "id",
    "lon",
    "lat",
    "date"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(dat)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "File ",
        basename(file),
        " is missing these columns: ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  dat %>%
    transmute(
      tags = as.character(id),
      long = as.numeric(lon),
      lat = as.numeric(lat),
      dTime = as.POSIXct(
        date,
        format = "%Y-%m-%dT%H:%M:%SZ",
        tz = "UTC"
      ),
      source_file = basename(file)
    )
}

#read + combine all observed tracks
tags <- bind_rows(
  lapply(
    track_files,
    read_whale_track
  )
) %>%
  filter(
    !is.na(tags),
    !is.na(long),
    !is.na(lat),
    !is.na(dTime)
  ) %>%
  distinct(
    tags,
    dTime,
    .keep_all = TRUE
  ) %>%
  arrange(
    tags,
    dTime
  ) %>%
  filter(
    !is.na(tags),
    !is.na(long),
    !is.na(lat),
    !is.na(dTime))%>%
  
  distinct(
    tags,
    dTime,
    .keep_all = TRUE) %>%
  arrange(
    tags,
    dTime)


#see how many locations each whale has
track_summary <- tags %>%
  group_by(tags) %>%
  summarise(
    n_locations = n(),
    
    start_time = min(dTime),
    end_time = max(dTime),
    
    median_interval_hours = if (n() > 1) {
      median(
        as.numeric(
          diff(dTime),
          units = "hours"
        ),
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    
    .groups = "drop"
  )

print(track_summary)

#land polygon for all whales to use
make_land_mask <- function(
    xlim_land,
    ylim_land
) {
  
  # maps/mapdata polygons can contain self-intersections that
  # S2 rejects. Use planar GEOS geometry while repairing them.
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  
  on.exit(
    sf::sf_use_s2(old_s2),
    add = TRUE
  )
  
  land_map <- maps::map(
    "worldHires",
    fill = TRUE,
    plot = FALSE,
    xlim = xlim_land,
    ylim = ylim_land
  )
  
  land_sf <- sf::st_as_sf(
    land_map
  )
  
  # The maps package coordinates are already longitude/latitude.
  # Assign EPSG:4326 only if no CRS has been recorded.
  if (is.na(sf::st_crs(land_sf))) {
    
    sf::st_crs(land_sf) <- 4326
    
  } else {
    
    land_sf <- sf::st_transform(
      land_sf,
      4326
    )
  }
  
  # Repair invalid rings and self-intersections
  land_sf <- sf::st_make_valid(
    land_sf
  )
  
  # GEOS zero-width buffering repairs some remaining ring issues
  land_sf <- sf::st_buffer(
    land_sf,
    dist = 0
  )
  
  # Dissolve the repaired polygons into one land geometry
  land_sf <- sf::st_union(
    land_sf
  )
  
  # Final repair after dissolving
  land_sf <- sf::st_make_valid(
    land_sf
  )
  
  return(land_sf)
}
  
land_sf <- make_land_mask(
  xlim_land = xlim_land,
  ylim_land = ylim_land
)

#CRW generation function
create_crw <- function(
    tags,
    tagid,
    land_sf,
    out_dir,
    n_sim = 200,
    xlim_land = c(-150, -90),
    ylim_land = c(0, 55),
    max_tries = 5000,
    plot_n = 20,
    seed = NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Temporarily disable spherical geometry for these
  # intersection operations
  old_s2 <- sf::sf_use_s2()
  
  sf::sf_use_s2(FALSE)
  
  on.exit(
    sf::sf_use_s2(old_s2),
    add = TRUE
  )
  
  message(
    "\nGenerating CRWs for whale: ",
    tagid
  )
  
#Extract whale 
  tag <- tags %>%
    filter(.data$tags == tagid) %>%
    arrange(dTime) %>%
    select(long, lat, dTime)
  
  tag <- tag %>%
    filter(
      !is.na(long), !is.na(lat), !is.na(dTime)) %>%
    distinct(
      dTime,
      .keep_all = TRUE)
  
n_tag <- nrow(tag)
  
  if (n_tag < 3) {
    stop(
      paste(
        "Whale",
        tagid,
        "has fewer than three valid locations."
      )
    )
  }

#Projects track into UTM coordinates 
tag_sf <- sf::st_as_sf(
  tag,
  coords = c(
    "long",
    "lat"
  ),
  crs = 4326,
  remove = FALSE
)

mean_lon <- mean(
  tag$long,
  na.rm = TRUE
)

mean_lat <- mean(
  tag$lat,
  na.rm = TRUE
)

utm_zone <- floor(
  (mean_lon + 180) / 6
) + 1

# Northern versus southern hemisphere UTM EPSG
if (mean_lat >= 0) {
  utm_epsg <- 32600 + utm_zone
} else {
  utm_epsg <- 32700 + utm_zone
}

tag_utm <- sf::st_transform(
  tag_sf,
  crs = utm_epsg
)

tag_xy <- sf::st_coordinates(
  tag_utm
)

tag$X <- tag_xy[, 1]
tag$Y <- tag_xy[, 2]

#movement calculation
observed_traj <- adehabitatLT::as.ltraj(
  xy = cbind(
    tag$X,
    tag$Y),
  date = tag$dTime,
  id = as.character(tagid))

tr <- observed_traj[[1]]

valid_steps <- which(
  !is.na(tr$dist) &
    is.finite(tr$dist) &
    tr$dist > 0 &
    !is.na(tr$rel.angle) &
    is.finite(tr$rel.angle))

if (length(valid_steps) < 2) {
  stop(
    paste(
      "Whale",
      tagid,
      "does not have enough usable movement steps."
    )
  )
}


# Initial heading of the observed whale
valid_headings <- which(
  !is.na(tr$abs.angle) &
    is.finite(tr$abs.angle)
)

if (length(valid_headings) == 0) {
  stop(
    paste(
      "Whale",
      tagid,
      "has no valid initial heading."
    )
  )
}

angle0 <- tr$abs.angle[
  valid_headings[1]
]


#prep simulation list
simulation_list <- vector(
  mode = "list",
  length = n_sim)

#generate CRW
for (k in seq_len(n_sim)) {
  
  if (
    k == 1 ||
    k %% 10 == 0 ||
    k == n_sim
  ) {
    message(
      "  Simulation ",
      k,
      " of ",
      n_sim
    )
  }
  
#Each CRW gets attached corresponding timestamp
  sim <- data.frame(
    whale_id = rep(
      as.character(tagid),
      n_tag),
    
    iteration = rep(
      k,
      n_tag),
    
    sim_id = rep(
      sprintf(
        "%s_crw_%03d",
        tagid,
        k),
      n_tag),
    
    step_id = seq_len(
      n_tag),
    
    choice_id = paste(
      tagid,
      format(
        tag$dTime,
        "%Y%m%dT%H%M%S"),
      sep = "_"),
    
    x = numeric(
      n_tag),
    
    y = numeric(
      n_tag),
    
    # Exact observed timestamps
    dTime = tag$dTime,
    
    stringsAsFactors = FALSE)
  
  
  # All simulations begin at the first observed location
  sim$x[1] <- tag$X[1]
  sim$y[1] <- tag$Y[1]
  angle <- angle0
  
#Generate each step loop
  for (j in 2:n_tag) {
    accepted <- FALSE
    tries <- 0
    while (!accepted) {
      tries <- tries + 1
      if (tries > max_tries) {
        stop(
          paste(
            "The land-rejection loop exceeded",
            max_tries,
            "attempts for whale",
            tagid,
            "simulation",
            k,
            "step",
            j)
        )
      }
#Random movement step from track
    sampled_step <- sample(
      valid_steps,
      size = 1)
    step_distance <- tr$dist[
      sampled_step
    ]
    
    relative_turn <- tr$rel.angle[
      sampled_step]
    
    
#Candidate heading where  current heading is not permanently changed unless the proposed location  accepted
    candidate_angle <- angle + relative_turn
    
#Proposed next location in UTM meters
    candidate_x <- sim$x[j - 1] +
      step_distance *
      cos(candidate_angle)
    
    candidate_y <- sim$y[j - 1] +
      step_distance *
      sin(candidate_angle)
    
    
#converts candidate point into lat/lon
    candidate_utm <- sf::st_sfc(
      sf::st_point(
        c(
          candidate_x,
          candidate_y)),
      crs = utm_epsg)
    
    candidate_ll <- sf::st_transform(
      candidate_utm,
      crs = 4326)
    
    candidate_coords <- sf::st_coordinates(
      candidate_ll)
    candidate_lon <- candidate_coords[1, 1]
    candidate_lat <- candidate_coords[1, 2]

#Reject locations outside  boundary
    outside_bounds <-
      candidate_lon < xlim_land[1] ||
      candidate_lon > xlim_land[2] ||
      candidate_lat < ylim_land[1] ||
      candidate_lat > ylim_land[2]
    
    if (outside_bounds) {
      next
    }
    
    
#Reject land points
    on_land <- suppressMessages(
      suppressWarnings(
        any(
          sf::st_intersects(
            candidate_ll,
            land_sf,
            sparse = FALSE
        )
     )
  )
)
    
    if (on_land) {
      next
    }
    
#Accept the proposed point
    sim$x[j] <- candidate_x
    sim$y[j] <- candidate_y
#only update  track heading after acceptance
    angle <- candidate_angle
    accepted <- TRUE
  }
}

#convert complete CRW  back to lon/lat
sim_sf <- sf::st_as_sf(
    sim,
    coords = c(
      "x",
      "y"),
    crs = utm_epsg,
    remove = FALSE)
  
  sim_ll <- sf::st_transform(
    sim_sf,
    crs = 4326)
  
  simulated_coordinates <- sf::st_coordinates(
    sim_ll)
  
  sim$long <- simulated_coordinates[, 1]
  sim$lat <- simulated_coordinates[, 2]
  sim$used <- 0L
  sim$location_type <- "available_crw"
  
simulation_list[[k]] <- sim %>%
  select(
    whale_id,
    iteration,
    sim_id,
    step_id,
    choice_id,
    dTime,
    long,
    lat,
    used,
    location_type
  )
}

#Combine simu for this whale
sim_all <- bind_rows(
  simulation_list)
#write csv
whale_csv <- file.path(
  out_dir,
  paste0(
    "crw_sim_",
    tagid,
    ".csv"))

write.csv(
  sim_all,
  whale_csv,
  row.names = FALSE)

#diagnostic plot
whale_png <- file.path(
  out_dir,
  paste0(
    "crw_sim_tracks_",
    tagid,
    ".png"))

plot_iterations <- sort(
  unique(
    sim_all$iteration))

set.seed(
  9000 + sum(
    utf8ToInt(
      as.character(tagid)
    )
  )
)

plot_iterations <- sample(
  plot_iterations,
  size = min(
    plot_n,
    length(plot_iterations)
  ),
  replace = FALSE
)

plot_buffer <- 0.35

plot_xlim <- range(
  c(
    tag$long,
    sim_all$long),
  na.rm = TRUE
) + c(
  -plot_buffer,
  plot_buffer)

plot_ylim <- range(
  c(
    tag$lat,
    sim_all$lat
  ),
  na.rm = TRUE
) + c(
  -plot_buffer,
  plot_buffer
)

png(
  whale_png,
  width = 1800,
  height = 1800,
  res = 300,
  bg = "white")

par(
  mar = c(
    4,
    4,
    3,
    1
  ),
  cex.axis = 0.9,
  cex.lab = 1)

maps::map(
  "worldHires",
  xlim = plot_xlim,
  ylim = plot_ylim,
  fill = TRUE,
  col = "grey90",
  border = "grey50")

maps::map.axes()

#Draw simulated tracks
for (current_iteration in plot_iterations) {
  current_sim <- sim_all %>%
    filter(
      iteration == current_iteration) %>%
    arrange(
      step_id)
  lines(
    current_sim$long,
    current_sim$lat,
    col = grDevices::adjustcolor(
      "deepskyblue3",
      alpha.f = 0.35),
    lwd = 0.7
  )
}

#Draw observed tracks
lines(
  tag$long,
  tag$lat,
  col = "black",
  lwd = 3.5)

points(
  tag$long,
  tag$lat,
  pch = 16,
  cex = 0.35,
  col = "black")

# Observed start
points(
  tag$long[1],
  tag$lat[1],
  pch = 16,
  cex = 1.3,
  col = "chartreuse3")

# Observed end
points(
  tag$long[n_tag],
  tag$lat[n_tag],
  pch = 16,
  cex = 1.3,
  col = "red3")

title(
  main = paste(
    "Observed track and CRW simulations:",
    tagid),
  cex.main = 0.9)

legend(
  "topright",
  legend = c(
    "CRW simulated tracks",
    "Observed track",
    "Observed start",
    "Observed end"),
  col = c(
    "deepskyblue3",
    "black",
    "chartreuse3",
    "red3"),
  lwd = c(
    1, 3.5, NA, NA),
  pch = c(
    NA, NA, 16, 16),
  pt.cex = c(NA, NA, 1.1, 1.1),
  bg = "white",
  box.lwd = 0.8,
  cex = 0.8)

dev.off()


message(
  "Finished whale ",
  tagid,
  "\n  CSV: ",
  whale_csv,
  "\n  Plot: ",
  whale_png
)

return(sim_all)
}


#CRW generation for all whales
tagids <- sort(
  unique(
    tags$tags
  ))

message(
  "\nWhale or track IDs detected:\n",
  paste(
    tagids,
    collapse = "\n"
  ))

all_crw_list <- vector(
  mode = "list",
  length = length(tagids))

names(all_crw_list) <- tagids


for (i in seq_along(tagids)) {
  
current_tagid <- tagids[i]
  
all_crw_list[[i]] <- create_crw(
  tags = tags,
  tagid = current_tagid,
  land_sf = land_sf,
  out_dir = out_dir,
  n_sim = n_crw,
  xlim_land = xlim_land,
  ylim_land = ylim_land,
  max_tries = 5000,
  plot_n = 20,

#reproducible random seed per whale
  seed = 202600 + i
)
}


#Combine allwhales CRWs

all_crw <- bind_rows(
  all_crw_list)

combined_crw_file <- file.path(
  out_dir,
  "crw_sim_all_four_whales.csv")

write.csv(
  all_crw,
  combined_crw_file,
  row.names = FALSE)

#Prepare  observed locations in matching format
observed_locations <- tags %>%
  group_by(tags) %>%
  arrange(
    dTime,
    .by_group = TRUE
  ) %>%
  mutate(
    whale_id = as.character(tags),
    
    iteration = 0L,
    
    sim_id = paste0(
      whale_id,
      "_observed"),
    step_id = row_number(),
    choice_id = paste(
      whale_id,
      format(
        dTime,
        "%Y%m%dT%H%M%S"),
      sep = "_"),
    
    used = 1L,
    
    location_type = "observed"
) %>%
ungroup() %>%
select(
  whale_id,
  iteration,
  sim_id,
  step_id,
  choice_id,
  dTime,
  long,
  lat,
  used,
  location_type
)


# 11. Combine observed + available locations
model_locations <- bind_rows(
  observed_locations,
  all_crw
) %>%
  filter(
    step_id > 1
) %>%
  arrange(
    whale_id,
    step_id,
    used,
    iteration
)

model_locations_file <- file.path(
  out_dir,
  "observed_and_crw_all_four_whales.csv")

write.csv(
  model_locations,
  model_locations_file,
  row.names = FALSE)

#  Quality summaries

cat(
  "\n============================================\n")

cat(
  "CRW generation complete\n")

cat(
  "============================================\n")

cat(
  "\nObserved locations by whale:\n")

print(
  observed_locations %>%
    count(
      whale_id,
      name = "observed_locations"))

cat(
  "\nSimulated CRW rows by whale:\n")

print(
  all_crw %>%
    count(
      whale_id,
      name = "simulated_rows"))

cat(
  "\nNumber of complete CRW tracks by whale:\n")

print(
  all_crw %>%
    distinct(
      whale_id,
      sim_id
    ) %>%
    count(
      whale_id,
      name = "complete_crw_tracks")
)

cat(
  "\nRows in combined model-location table:\n")

print(
  model_locations %>%
    count(
      whale_id,
      used,
      location_type)
)

cat(
  "\nCombined CRW file:\n",
  combined_crw_file,
  "\n")

cat(
  "\nObserved + available model table:\n",
  model_locations_file,
  "\n")

