require(adehabitatLT) #turns a sequence of points into a trajectory object and computes step lengths + turning angles.
require(maps) 
require(mapdata)
require(maptools) #builds a land polygon and tests if a simulated point is on land.
require(sp)
require(raster)
library(dplyr) #sorting/filtering

#converting tag to useable format
hbtag <- "C:/github/Whale-SDM/Output/Track processing output/254025-RawArgos_pred_6h.csv"
dat <- read.csv(hbtag, stringsAsFactors = FALSE)
#Time column to POSIXct format
dat

#takes one observed track to generate pseudo-absences
#tags = dataframe containing all animals, tagid=from blue_whale_dataset.r wrapper
createCRW <- function(tags, tagid, n.sim = 200, reverse = FALSE) {
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

  
  CC.map <- maps::map('worldHires', fill = TRUE, col = 'transparent', plot = FALSE, ylim = c(0, 50))
    #Loads coastline polygons from the maps package (high-res world), but doesn’t plot them.
        #plot = FALSE returns the map object instead of drawing it.
        #ylim = c(0, 50) restricts latitudes included (saves work; focuses region).
    #New object: CC.map (a map object containing polygon coordinates and names)
  CC.IDs <- sapply(strsplit(CC.map$names, ":"), function(x) x[1]) #Creates polygon IDs from the map names.
  CC.sp <- map2SpatialPolygons(CC.map, IDs = CC.IDs, proj4string = CRS("+proj=longlat +datum=WGS84"))
    #Converts the map polygons into an sp::SpatialPolygons object with a WGS84 lon/lat CRS.
  
  sim.alltags <- NULL #Initializes an empty container & will repeatedly rbind() simulation results onto this
  
  for (k in 1:n.sim) {
    print(sprintf("  k'th simulation: %d", k))
        #Repeats the entire simulation process n.sim times, k is the loop counter (thus not stored after the loop finishes).
    
    n.tag <- nrow(tag)
    sim <- data.frame(x = numeric(n.tag), y = numeric(n.tag), t = as.POSIXct(rep(NA, n.tag)), flag = NA, iteration = rep(k, n.tag))
      #Creates a dataframe sim with the same number of rows as the observed track, intended to hold the simulated track.
    sim[1, 'x'] <- tag$long[1]
    sim[1, 'y'] <- tag$lat[1]
    sim[1, 't'] <- tr1$date[1]
    angle <- tr1[2, 'abs.angle']
      #Initial conditions:
          #The simulated track starts at the observed first location.
          #The simulated start time is the observed first time.
          #The initial heading angle is pulled from the observed trajectory (abs.angle at step 2).
      #New objects
          #sim now has row 1 initialized.
          #New: angle (current heading, in radians)
    
    for (j in 2:n.tag) {
      on.land <- TRUE
      while (on.land) {
       #For each position index j (2..n.tag), it keeps trying to propose a step until that proposal lands in the ocean.
          #on.land <- TRUE starts a “rejection sampling” loop.
          #while(on.land) repeats until you set on.land <- FALSE.
         
          i <- sample(i.tr1.nona, 1)
        dist <- tr1[i, 'dist']
        angle <- angle + tr1[i, 'rel.angle']
        #Randomly picks one observed step index i, and uses:
          #its step length (dist)
          #its turning angle (rel.angle)
          #to update the heading.
        
        x <- sim[j - 1, 'x'] + dist * cos(angle)
        y <- sim[j - 1, 'y'] + dist * sin(angle)
        #Moves from the previous simulated point to a new proposed point using basic trig, creates new y and x
        
        pt <- SpatialPoints(matrix(c(x, y), nrow = 1), proj4string = CRS("+proj=longlat +datum=WGS84"))
        place <- over(pt, CC.sp)
        #Checks whether the proposed point is on land.
          #Creates a spatial point object at (x,y).
          #Uses over() to see if it lies inside any land polygon.
        #New objects
          #pt (SpatialPoints)
          #place (result of spatial overlay)
        
        if (is.na(place)) {
          sim[j, 'x'] <- x
          sim[j, 'y'] <- y
          #Accept the proposed point and store it in the simulated track.
          
          dt_i <- tr1[i, 'dt']
          if (is.na(dt_i) || dt_i > (60 * 60 * 24 * 30)) {
            dt_i <- 60 * 60 * 24  # default to 1 day
          }
          sim[j, 't'] <- sim[j - 1, 't'] + dt_i
          on.land <- FALSE
          #Advances simulated time by the time gap (dt) from the sampled observed step.
              #If dt is missing or absurdly huge (> 30 days), force dt to 1 day.
              #Set the simulated time for row j.
              #Set on.land <- FALSE to exit the while-loop (we accepted a point).
        }
      }
    }
    
    # Calculate flag
    tagstart <- rbind(tag[1, ], tag[2, ], tag[n.tag, ])
    tagsim <- rbind(sim[1, ], sim[2, ], sim[n.tag, ])
      #Creates two tiny 3-point “tracks”:
        #tagstart: first, second, last observed point
        #tagsim: first, second, last simulated point
        #Why? Because as.ltraj() needs at least 3 points to compute a turning angle/step structure.
    
    tagstart$lat[2] <- tagstart$lat[1] - 0.01
    tagsim$y[2] <- tagsim$y[1] - 0.01
    #This is a hack to prevent degenerate geometry issues when computing angles.
      #It nudges the second point slightly south to ensure something about the angle computation isn’t undefined.
      #It’s… hacky, but it’s trying to avoid an edge-case where angle calculations fail (e.g., if points coincide).
    
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
  
  
  
  write.csv(sim.alltags, out.alltags.csv, row.names = FALSE)
  #Writes the big dataframe sim.alltags to a CSV file on disk.
    #sim.alltags = all simulated points from all n.sim iterations concatenated together.
    #out.alltags.csv = the filename string created earlier (based on out.dir and tagid).
  
  # Output map to PNG
  out.png <- sprintf('%s/crw_sim_%s.png', out.dir, tagid)
  res <- 72
  png(filename = out.png, width = 8*res, height = 8*res, res = res)
  
  # Redraw the map and the original track
  maps::map('worldHires', xlim=c(-140, -100), ylim=c(15, 50))
  map.axes()
  lines(tag$lon, tag$lat, col='grey')
  points(sim[1, 'x'], sim[1, 'y'], col='blue', pch=2, cex=2)
  title(main = sprintf("CRW Simulations for tagid %s", tagid))
  lines(sim$x, sim$y, col='black')
  
  dev.off()
  return(sim.alltags)
}
