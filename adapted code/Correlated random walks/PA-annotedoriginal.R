require(adehabitatLT) #turns a sequence of points into a trajectory object and computes step lengths + turning angles.
require(maps) 
require(mapdata)
library(sf)
library(dplyr) #sorting/filtering

#Example csv of humpback
hbtag <- read.csv("C:/github/Whale-SDM/Output/Track processing output/254026-RawArgos_pred_6h.csv", stringsAsFactors = FALSE)


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
createCRW <- function(tags, tagid, n.sim = 200, reverse = FALSE,
                      xlim_land = c(-150, -90), ylim_land = c(0, 55),
                      max_tries = 5000,
                      png_px = 2400, png_bg = "white") {
  
  # --- sf s2 off for speed/robustness in some intersection ops ---
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  
  message("Running tag: ", tagid)
  
  # ---- Extract + clean observed track ----
  tag <- tags[which(tags$tags == tagid), c("long", "lat", "dTime")]
  tag <- dplyr::arrange(tag, dTime)
  
  if (reverse) {
    tag <- tag[seq(nrow(tag), 1), ]
  }
  
  # Remove non-unique and NA values
  tag <- tag[!duplicated(tag$dTime), ]
  tag <- tag[!is.na(tag$lat) & !is.na(tag$long) & !is.na(tag$dTime), ]
  
  if (nrow(tag) < 3) {
    warning(paste("Not enough points to simulate tag:", tagid))
    return(NULL)
  }
  
  # ---- Land polygon (California Current-ish window) ----
  CC.map <- maps::map(
    "worldHires",
    fill = TRUE,
    plot = FALSE,
    xlim = xlim_land,
    ylim = ylim_land
  )
  CC.sf <- sf::st_as_sf(CC.map)
  CC.sf <- sf::st_set_crs(CC.sf, 4326)
  CC.sf <- sf::st_make_valid(CC.sf)
  CC.sf <- sf::st_union(CC.sf)  # dissolve for faster intersects
  
  # ---- Project observed track to a metric CRS (UTM) ----
  tag_sf <- sf::st_as_sf(tag, coords = c("long", "lat"), crs = 4326, remove = FALSE)
  
  mean_lon <- mean(tag$long, na.rm = TRUE)
  utm_zone <- floor((mean_lon + 180) / 6) + 1
  utm_crs  <- paste0("+proj=utm +zone=", utm_zone, " +datum=WGS84 +units=m +no_defs")
  
  tag_utm <- sf::st_transform(tag_sf, crs = utm_crs)
  xy <- sf::st_coordinates(tag_utm)
  tag$X <- xy[, 1]
  tag$Y <- xy[, 2]
  
  # ---- Build trajectory in meters ----
  tr <- adehabitatLT::as.ltraj(cbind(tag$X, tag$Y), date = tag$dTime, id = tagid)
  tr1 <- tr[[1]]
  
  i.tr1.nona <- which(!is.na(tr1[["dist"]]) & !is.na(tr1[["rel.angle"]]) & !is.na(tr1[["dt"]]))
  if (length(i.tr1.nona) < 2) {
    warning(paste("Not enough usable steps for tag:", tagid))
    return(NULL)
  }
  
  # ---- Validate the observed first-step distance (fixes the stray `next`) ----
  # We'll compute a minimal 3-point "start" trajectory in UTM to get dist[2].
  n.tag <- nrow(tag)
  tagstart <- rbind(tag[1, ], tag[2, ], tag[n.tag, ])
  
  # tiny jitter to avoid degenerate geometry
  tagstart$Y[2] <- tagstart$Y[1] - 10  # 10 m south
  
  trstart <- adehabitatLT::as.ltraj(cbind(tagstart$X, tagstart$Y), date = tagstart$dTime, id = tagid)
  d0 <- trstart[[1]]$dist[2]
  if (is.na(d0) || d0 <= 0) {
    warning(paste("Bad initial step distance for tag:", tagid, "— skipping tag"))
    return(NULL)
  }
  
  # ---- Output filenames ----
  out.alltags.csv <- sprintf("%s/crw_sim_%s.csv", out.dir, tagid)
  if (reverse) out.alltags.csv <- sprintf("%s/crw_sim_reverse_%s.csv", out.dir, tagid)
  out.png <- sprintf("%s/crw_sim_%s.png", out.dir, tagid)
  if (reverse) out.png <- sprintf("%s/crw_sim_reverse_%s.png", out.dir, tagid)
  
  # ---- Run simulations ----
  sim.alltags <- NULL
  
  # initial heading (radians) from observed trajectory (abs.angle at step 2)
  angle_candidates <- tr1[["abs.angle"]]
  angle0 <- angle_candidates[which(!is.na(angle_candidates))[1]]
  
  if (is.na(angle0)) {
    warning(paste("No valid initial heading for tag:", tagid))
    return(NULL)
  }  
  for (k in 1:n.sim) {
    message("  simulation k = ", k)
    
    n.tag <- nrow(tag)
    
    sim <- data.frame(
      x = numeric(n.tag),   # UTM meters
      y = numeric(n.tag),   # UTM meters
      t = as.POSIXct(rep(NA, n.tag), tz = "UTC"),
      flag = NA_real_,
      iteration = rep(k, n.tag)
    )
    
    # anchor start to observed first location/time
    sim[1, "x"] <- tag$X[1]
    sim[1, "y"] <- tag$Y[1]
    sim[1, "t"] <- tr1$date[1]
    angle <- angle0
    
    for (j in 2:n.tag) {
      on.land <- TRUE
      tries <- 0
      
      while (on.land) {
        tries <- tries + 1
        if (tries > max_tries) {
          stop(paste(
            "Stuck in land-rejection loop at step", j,
            "simulation", k,
            "- check step lengths or land mask"
          ))
        }
        
        # sample an empirical movement step
        i <- sample(i.tr1.nona, 1)
        dist  <- tr1[i, "dist"]                 # meters now
        angle <- angle + tr1[i, "rel.angle"]    # radians
        
        # propose new point in UTM meters
        x <- sim[j - 1, "x"] + dist * cos(angle)
        y <- sim[j - 1, "y"] + dist * sin(angle)
        
        # convert candidate point to lon/lat for bounds + land intersection
        pt_utm <- sf::st_sfc(sf::st_point(c(x, y)), crs = utm_crs)
        pt_ll  <- sf::st_transform(pt_utm, 4326)
        ll <- sf::st_coordinates(pt_ll)
        x_ll <- ll[1, 1]
        y_ll <- ll[1, 2]
        
        # quick geographic bounds reject
        if (x_ll < xlim_land[1] || x_ll > xlim_land[2] || y_ll < ylim_land[1] || y_ll > ylim_land[2]) {
          next
        }
        
        on.land <- suppressMessages(suppressWarnings(
          any(sf::st_intersects(pt_ll, CC.sf, sparse = FALSE))
        ))
        
        if (!on.land) {
          sim[j, "x"] <- x
          sim[j, "y"] <- y
          
          dt_i <- tr1[i, "dt"]
          if (is.na(dt_i) || dt_i > (60 * 60 * 24 * 30)) dt_i <- 60 * 60 * 24
          sim[j, "t"] <- sim[j - 1, "t"] + dt_i
        }
      }
    }
    
    # ---- Flag calculation (compare first move observed vs simulated) ----
    tagsim <- rbind(sim[1, ], sim[2, ], sim[n.tag, ])
    
    # jitter simulated mini-track too
    tagsim$y[2] <- tagsim$y[1] - 10  # 10 m south
    
    trsim <- adehabitatLT::as.ltraj(cbind(tagsim$x, tagsim$y), date = tagsim$t, id = tagid)
    
    distdiff <- abs(trstart[[1]]$dist[2] - trsim[[1]]$dist[2])
    angdiff  <- abs(trstart[[1]]$rel.angle[2] * 180 / pi - trsim[[1]]$rel.angle[2] * 180 / pi)
    
    sim$flag <- distdiff / trstart[[1]]$dist[2] * 3 + angdiff / 45
    
    # ---- Convert simulated track to lon/lat columns for output/plot ----
    sim_pts_utm <- sf::st_as_sf(sim, coords = c("x", "y"), crs = utm_crs)
    sim_pts_ll  <- sf::st_transform(sim_pts_utm, 4326)
    xy_ll <- sf::st_coordinates(sim_pts_ll)
    sim$long <- xy_ll[, 1]
    sim$lat  <- xy_ll[, 2]
    
    sim.alltags <- rbind(sim.alltags, sim)
  }
  
  message("Finished sims for tag ", tagid, " — writing outputs now.")
  sim_out <- sim.alltags[, c("iteration", "t", "long", "lat", "flag")]
  names(sim_out)[names(sim_out) == "t"] <- "dTime"
  utils::write.csv(sim_out, out.alltags.csv, row.names = FALSE)
  
  # ---- Plot PNG with sane size (fixes the huge image issue) ----
  grDevices::png(out.png, width = png_px, height = png_px, bg = png_bg)
  
  # extent based on observed + simulated lon/lat (small buffer in degrees)
  xlim <- range(c(tag$long, sim.alltags$long), na.rm = TRUE) + c(-1, 1)
  ylim <- range(c(tag$lat,  sim.alltags$lat),  na.rm = TRUE) + c(-1, 1)
  
  maps::map("worldHires", xlim = xlim, ylim = ylim, fill = TRUE, col = "grey90")
  maps::map.axes()
  
  # observed track in red
  lines(tag$long, tag$lat, col = "red", lwd = 2)
  
  # plot a subset of simulated tracks (blue)
  set.seed(1)
  iters <- sort(unique(sim.alltags$iteration))
  iters_to_plot <- sample(iters, size = min(10, length(iters)), replace = FALSE)
  
  for (kk in iters_to_plot) {
    ssub <- sim.alltags[sim.alltags$iteration == kk, ]
    lines(ssub$long, ssub$lat, col = "blue", lwd = 1)
  }
  
  points(tag$long[1], tag$lat[1], pch = 16, cex = 1.2)
  title(main = sprintf("Observed + %d CRW sims (tagid %s)", length(iters_to_plot), tagid))
  
  legend(
    "topright",
    legend = c("Observed track", "Simulated tracks (subset)", "Observed start point"),
    col = c("red", "blue", "black"),
    lwd = c(2, 1, NA),
    pch = c(NA, NA, 16),
    pt.cex = c(NA, NA, 1),
    bg = "white",
    box.lwd = 0.8
  )
  
  grDevices::dev.off()
  
  if (interactive()) {
    maps::map("worldHires", xlim = xlim, ylim = ylim, fill = TRUE, col = "grey90")
    maps::map.axes()
    lines(tag$long, tag$lat, col = "red", lwd = 2)
    for (kk in iters_to_plot) {
      ssub <- sim.alltags[sim.alltags$iteration == kk, ]
      lines(ssub$long, ssub$lat, col = "blue", lwd = 1)
    }
    points(tag$long[1], tag$lat[1], pch = 16, cex = 1.2)
    title(main = sprintf("Observed + %d CRW sims (tagid %s)", length(iters_to_plot), tagid))
  }
  
  return(sim.alltags)
}



# ---- TEST RUN (outside the function) ----
tagid <- unique(tags$tags)[1]
sim_test <- createCRW(tags, tagid, n.sim = 10)

#full CRW set for one whale
tagid <- unique(tags$tags)[1]
sim_all <- createCRW(tags, tagid, n.sim = 20)

#run all tags
tagids <- unique(tags$tags)
all_sims <- lapply(tagids, function(id) {
  createCRW(tags, id, n.sim = 200)
})
# Optional: combine everything into one big data frame
all_sims_df <- dplyr::bind_rows(all_sims)


