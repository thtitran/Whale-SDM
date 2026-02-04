require(adehabitatLT) #turns a sequence of points into a trajectory object and computes step lengths + turning angles.
require(maps) 
require(mapdata)
library(sf)
library(dplyr) #sorting/filtering

#Example csv of humpback
hbtag <- read.csv("C:/github/Whale-SDM/Output/Track processing output/254025-RawArgos_pred_6h.csv", stringsAsFactors = FALSE)


#Humpback csv requires conversion to UTC POSIXct format:
hbtag$dTime <- as.POSIXct(hbtag$date, tz = "UTC")
#Changing column names to suite code
hbtag$tags <- hbtag$id
hbtag$long <- hbtag$lon
names(hbtag)

#New tag data frame
tags <- hbtag[, c("tags", "long", "lat", "dTime")]

out.dir <- "C:/github/Whale-SDM/Output/CRW_test"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)


#takes one observed track to generate pseudo-absences
#tags = dataframe containing all animals, tagid=from blue_whale_dataset.r wrapper
createCRW <- function(tags, tagid, n.sim = 200, reverse = FALSE) {
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  print(tagid)
  tag <- tags[which(tags$tags == tagid), c('long', 'lat', 'dTime')] #extracts an individual’s track, tags$tags == tagid creates a logical vector
  #which(...) converts that logical vector into integer row indices (e.g., c(1, 2, 5, 9, ...))
  #Dataframe indexing: tags[ rows , cols ]
  #Rows: which(tags$tags == tagid)
  #Cols: c('long','lat','dTime') (a character vector of column names)
  tag <- arrange(tag, dTime) #Sorts the rows of tag by time (dTime) ascending to ensure correct chronological order
  #tag: From  big dataframe tags, selects only  rows where the animal ID = tagid, and keeps only three columns: long, lat, dTime. Stores  result in  new object: 'tag'
  out.alltags.csv <- sprintf('%s/crw_sim_%s.csv', out.dir, tagid) #Constructs an output filename string for where the simulated points will be written
  #sprintf() formats strings like C’s printf
  if (reverse) {
    out.alltags.csv <- sprintf('%s/crw_sim_reverse_%s.csv', out.dir, tagid)
    tag <- tag[seq(dim(tag)[1], 1), ] #Reorders tag so its rows go from last to first. reverse-walk simulations (sometimes used as a check or for symmetric availability, depending on method).
  }
  #Starts a conditional block that only runs if reverse is TRUE.
  #If reversing, change the output filename so the file name reflects that it’s a reverse simulation
  
  # Remove non-unique and NA values
  tag <- tag[!duplicated(tag$dTime), ] #Keeps only the first row for each unique timestamp dTime and drops dupes
  tag <- tag[!is.na(tag$lat) & !is.na(tag$long) & !is.na(tag$dTime), ] #Drops any rows where latitude, longitude, or datetime is missing.
  
  #If after cleaning you have fewer than 3 points, it stops and returns NULL because turning angle needs 3 points min
  if (nrow(tag) < 3) {
    warning(paste("Not enough points to simulate tag:", tagid))
    return(NULL)
  }
  
  
  # Create trajectory
  #Creates  trajectory object (ltraj) that includes step-by-step movement calcs then extracts  first (and usually only) trajectory element
  tr <- as.ltraj(cbind(tag$long, tag$lat), date = tag$dTime, id = tagid) #makes a 2-column matrix of coordinates, from adehabitatLT constructs an ltraj object
  #tr is actually a list of trajectories (one per id).
  tr1 <- tr[[1]]
  #tr1 (the actual trajectory data for this individual)
  #tr1 contains derived columns such as:
  #dist = step length between successive points
  #dt = time difference between points (seconds)
  #rel.angle = turning angle between successive steps (radians)
  #abs.angle = heading direction (radians)
  
  #Creates an index of rows in tr1 that have all the movement quantities needed for simulation.
  #tr1[['dist']] accesses the dist column, Using [['...']] instead of $dist is just another access style
  #!is.na(...) keeps only rows where the value exists.
  #which(...) converts TRUE/FALSE into row numbers.
  i.tr1.nona <- which(!is.na(tr1[['dist']]) & !is.na(tr1[['rel.angle']]) & !is.na(tr1[['dt']]))
  #New object: i.tr1.nona (integer vector of row indices)
  #Purpose: the simulation will randomly sample from these indices to get:
  #a step length (dist)
  #a turning angle (rel.angle)
  #a time step (dt)
  if (length(i.tr1.nona) < 2) {
    warning(paste("Not enough usable steps for tag:", tagid))
    return(NULL)
  }
  
  
  CC.map <- maps::map(
    "worldHires",
    fill = TRUE,
    plot = FALSE,
    xlim = c(-150, -90),
    ylim = c(0, 55)
  )
  
  CC.sf <- sf::st_as_sf(CC.map)
  CC.sf <- sf::st_set_crs(CC.sf, 4326)
  CC.sf <- sf::st_make_valid(CC.sf)   # optional but recommended
  CC.sf <- sf::st_union(CC.sf)        # dissolve to one geometry (faster)
  #Loads coastline polygons from the maps package (high-res world), but doesn’t plot them.
  #plot = FALSE returns the map object instead of drawing it.
  #ylim = c(0, 50) restricts latitudes included (saves work; focuses region).
  #New object: CC.map (a map object containing polygon coordinates and names)
  
  sim.alltags <- NULL #Initializes an empty container & will repeatedly rbind() simulation results onto this
  
  for (k in 1:n.sim) {
    print(sprintf("  k'th simulation: %d", k))
    #Repeats the entire simulation process n.sim times, k is the loop counter (thus not stored after the loop finishes).
    
    n.tag <- nrow(tag)
    
    sim <- data.frame(x = numeric(n.tag),y = numeric(n.tag),t = as.POSIXct(rep(NA, n.tag), tz = "UTC"),flag = NA_real_,iteration = rep(k, n.tag))
    #Creates a dataframe sim with the same number of rows as the observed track, intended to hold the simulated track.
    sim[1, "x"] <- tag$long[1]
    sim[1, "y"] <- tag$lat[1]
    sim[1, "t"] <- tr1$date[1]
    angle <- tr1[2, "abs.angle"]
    #Initial conditions:
    #The simulated track starts at the observed first location.
    #The simulated start time is the observed first time.
    #The initial heading angle is pulled from the observed trajectory (abs.angle at step 2).
    #New objects
    #sim now has row 1 initialized.
    #New: angle (current heading, in radians)
    
    for (j in 2:n.tag) {
      # We’re building ONE simulated track, point by point.
      # j indexes position along this simulated track:
      #   j = 1 is the anchored start point (already set)
      #   j = 2..n.tag are simulated steps
      
      on.land <- TRUE
      tries <- 0
      # Rejection sampler:
      # Keep proposing a candidate point until it is NOT on land.
      while (on.land) {
        
        # Safety cap: prevents infinite loops if something goes wrong
        # (e.g., step sizes too large, land mask mismatch, etc.)
        tries <- tries + 1
        if (tries > 5000) {
          stop(paste(
            "Stuck in land-rejection loop at step", j,
            "simulation", k,
            "- check step lengths or land mask"
          ))
        }
        
        # Sample ONE observed movement “step” from the empirical pool:
        # i.tr1.nona contains indices where dist, dt, rel.angle are all defined.
        i <- sample(i.tr1.nona, 1)
        
        # Step length (dist) and turning angle (rel.angle) come from the observed track.
        dist  <- tr1[i, "dist"]
        angle <- angle + tr1[i, "rel.angle"]
        # angle is the *current heading* (radians) we carry forward through the simulation.
        # We update it by adding a sampled turning angle.
        
        # Propose a new location using trig (a correlated random walk step):
        x <- sim[j - 1, "x"] + dist * cos(angle)
        y <- sim[j - 1, "y"] + dist * sin(angle)
        # NOTE: because x/y here are lon/lat, dist is in degrees (not meters).
        
        # Land test using sf (modern spatial workflow):
        # 1) make a point geometry at the proposed (x, y)
        # 2) check whether it intersects any land polygon in CC.sf
        if (x < -150 || x > -90 || y < 0 || y > 55) next #Quick bounds reject (saves tons of land tests)
        pt <- sf::st_sfc(sf::st_point(c(x, y)), crs = 4326)
        on.land <- suppressMessages(suppressWarnings(any(sf::st_intersects(pt, CC.sf, sparse = FALSE))))
        # on.land = TRUE  → reject and try again
        # on.land = FALSE → accept below
        
        if (!on.land) {
          # Accept the proposed point: store it in the simulated track
          sim[j, "x"] <- x
          sim[j, "y"] <- y
          
          # Advance time by a sampled time step from the observed track
          dt_i <- tr1[i, "dt"]
          
          # Guardrail: if dt is missing or absurdly huge, force it to 1 day
          if (is.na(dt_i) || dt_i > (60 * 60 * 24 * 30)) {
            dt_i <- 60 * 60 * 24
          }
          
          sim[j, "t"] <- sim[j - 1, "t"] + dt_i
          
          # We do NOT need to manually set on.land <- FALSE here.
          # It is already FALSE (that’s why we entered this if-block),
          # so the while(on.land) will stop naturally on the next check.
        }
      }
    }
    
    n.tag <- nrow(tag) # repeat in case n.tag from loop is not picked up
    
    # Calculate flag
    # Create minimal 3-point tracks (as.ltraj needs >= 3 points to compute angles)
    tagstart <- rbind(tag[1, ], tag[2, ], tag[n.tag, ])
    tagsim   <- rbind(sim[1, ], sim[2, ], sim[n.tag, ])
    
    # Small jitter to avoid degenerate geometry (e.g., identical points -> undefined angles)
    # Nudge the second point slightly south so angle calculations don't blow up
    tagstart$lat[2] <- tagstart$lat[1] - 0.01
    tagsim$y[2]     <- tagsim$y[1] - 0.01
    
    
    trstart <- as.ltraj(cbind(tagstart$long, tagstart$lat), date = tagstart$dTime, id = tagid)
    trsim <- as.ltraj(cbind(tagsim$x, tagsim$y), date = tagsim$t, id = tagid)
    #Creates two trajectory objects: one for observed mini-track and one for simulated mini-track.
    
    distdiff <- abs(trstart[[1]]$dist[2] - trsim[[1]]$dist[2])
    angdiff <- abs(trstart[[1]]$rel.angle[2] * 180 / pi - trsim[[1]]$rel.angle[2] * 180 / pi)
    #Compares the first step of observed vs simulated:
    #dist[2] corresponds to the step from point 1 → 2 (because the first row often has NA dist).
    
    #rel.angle[2] is a turning angle in radians; they convert to degrees by * 180 / pi.
    #So:
    #distdiff = absolute difference in step length
    #angdiff = absolute difference in turning angle
    
    sim$flag <- distdiff / trstart[[1]]$dist[2] * 3 + angdiff / 45
    #Computes a single score (“flag”) and assigns it to every row in this simulated track.
    #penalize distance mismatch relative to observed distance, scaled by 3
    #plus penalize angle mismatch scaled by 45 degrees
    #So lower flag ≈ “the first move looks like the observed first move”.
    
    sim.alltags <- rbind(sim.alltags, sim)
    #Appends this simulation’s dataframe sim onto the big collector sim.alltags and ends the outer loop iteration k.
  }
  
  
  message("Finished sims for tag ", tagid, " — writing outputs now.")
  write.csv(sim.alltags, out.alltags.csv, row.names = FALSE)
  #Writes the big dataframe sim.alltags to a CSV file on disk.
  #sim.alltags = all simulated points from all n.sim iterations concatenated together.
  #out.alltags.csv = the filename string created earlier (based on out.dir and tagid).
  #row.names = FALSE prevents R from adding a leftmost “1,2,3…” column to the CSV.
  
  # Output map to PNG
  out.png <- sprintf('%s/crw_sim_%s.png', out.dir, tagid)
  res <- 72 #Sets a resolution value (dots per inch) used for the PNG device.
  png(filename = out.png, width = 8*res, height = 8*res, res = res)
  
  # Choose map extent based on observed + simulated points (with a small buffer)
  xlim <- range(c(tag$long, sim.alltags$x), na.rm = TRUE) + c(-1, 1)
  ylim <- range(c(tag$lat,  sim.alltags$y), na.rm = TRUE) + c(-1, 1)
  maps::map("worldHires", xlim=xlim, ylim=ylim, fill=TRUE, col="grey90")
  # Redraw the map and the original track
  map.axes()
  
  # Observed track (thicker grey)
  lines(tag$long, tag$lat, col = "grey40", lwd = 2)
  
  # Sample up to 10 simulated iterations and plot them as thin black lines
  set.seed(1)  # so the same 10 sims get drawn each run (optional)
  iters <- sort(unique(sim.alltags$iteration))
  iters_to_plot <- sample(iters, size = min(10, length(iters)), replace = FALSE)
  
  for (kk in iters_to_plot) {
    ssub <- sim.alltags[sim.alltags$iteration == kk, ]
    lines(ssub$x, ssub$y, col = "black", lwd = 1)
  }
  
  # Start point marker (observed)
  points(tag$long[1], tag$lat[1], pch = 16, cex = 1.2)
  
  title(main = sprintf("Observed + %d CRW sims (tagid %s)", length(iters_to_plot), tagid))
  
  dev.off()
  return(sim.alltags)
}

# ---- TEST RUN (outside the function) ----
tagid <- unique(tags$tags)[1]
sim_test <- createCRW(tags, tagid, n.sim = 1)

#full CRW set for one whale
tagid <- unique(tags$tags)[1]
sim_all <- createCRW(tags, tagid, n.sim = 200)

#run all tags
tagids <- unique(tags$tags)
all_sims <- lapply(tagids, function(id) {
  createCRW(tags, id, n.sim = 200)
})
# Optional: combine everything into one big data frame
all_sims_df <- dplyr::bind_rows(all_sims)
