Sys.setenv(COPERNICUSMARINE_SERVICE_USERNAME = "ttran5")
Sys.setenv(COPERNICUSMARINE_SERVICE_PASSWORD = "Deusomnissiah1!")

toolbox <- "C:/Users/Theo/Downloads/copernicusmarine.exe"

out_dir_zos  <- "C:/github/Whale-SDM/data/copernicus/zos"
out_dir_phys <- "C:/github/Whale-SDM/data/copernicus/surface_phys"

dir.create(out_dir_zos,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_phys, recursive = TRUE, showWarnings = FALSE)

#ZOS
cmd_zos <- paste(
  shQuote(toolbox),
  "subset",
  "--dataset-id cmems_mod_glo_phy_my_0.083deg_P1D-m",
  "--variable zos",
  "--start-datetime 2024-03-01T00:00:00",
  "--end-datetime 2024-09-30T00:00:00",
  "--minimum-longitude -130",
  "--maximum-longitude -110",
  "--minimum-latitude 25",
  "--maximum-latitude 45",
  "-o", shQuote(out_dir_zos)
)

system(cmd_zos)

# Surface thetao/uo/vo/so
cmd_phys <- paste(
  shQuote(toolbox),
  "subset",
  "--dataset-id cmems_mod_glo_phy_my_0.083deg_P1D-m",
  "--variable thetao",
  "--variable uo",
  "--variable vo",
  "--variable so",
  "--minimum-depth 0.49",
  "--maximum-depth 0.50",
  "--start-datetime 2024-03-01T00:00:00",
  "--end-datetime 2024-09-30T00:00:00",
  "--minimum-longitude -130",
  "--maximum-longitude -110",
  "--minimum-latitude 25",
  "--maximum-latitude 45",
  "-o", shQuote(out_dir_phys)
)

system(cmd_phys)

library(terra)

f <- list.files("C:/github/Whale-SDM/data/copernicus",
                pattern = "\\.nc$",
                full.names = TRUE)

r <- rast(f[1])
r
names(r)
time(r)


#NOAA 200 m contour
library(sf)
library(terra)


url_200m <- paste0(
  "https://gis.ngdc.noaa.gov/arcgis/rest/services/",
  "nccos/USWestCoast_DeepSeaCoralsSponges/MapServer/2/query?",
  "where=1%3D1&",
  "outFields=*&",
  "returnGeometry=true&",
  "f=geojson"
)

shelf_200 <- st_read(url_200m, quiet = FALSE)
shelf_200
st_crs(shelf_200)

if (is.na(st_crs(shelf_200))) {
  st_crs(shelf_200) <- 4326
}

# Study area in lon/lat
bbox_ll <- st_as_sfc(
  st_bbox(c(
    xmin = -130,
    xmax = -110,
    ymin = 25,
    ymax = 45
  ), crs = 4326)
)

# shelf_200 = your NOAA 200 m contour sf object in lon/lat
# if needed:
# shelf_200 <- st_read(url_200m, quiet = FALSE)

# crop line to study area
shelf_200_ca <- st_intersection(shelf_200, bbox_ll)

#Project to a metric CRS
# NAD83/Conus Albers
bbox_5070      <- st_transform(bbox_ll, 5070)
shelf_200_5070 <- st_transform(shelf_200_ca, 5070)

#Build raster template in meters
#  GLORYS-like ~8 km
res_m <- 8000

# roughly 0.25-ish cells  as a coarse approximation.
# res_m <- 25000

bb <- st_bbox(bbox_5070)

r_template <- rast(
  xmin = bb["xmin"],
  xmax = bb["xmax"],
  ymin = bb["ymin"],
  ymax = bb["ymax"],
  resolution = res_m,
  crs = "EPSG:5070"
)

#rasterize the shelf line
shelf_vect <- vect(shelf_200_5070)

shelf_r <- rasterize(
  shelf_vect,
  r_template,
  field = 1,
  background = NA
)

#distance raster in meters
dist_shelf_m <- distance(shelf_r)

plot(dist_shelf_m, main = "Distance to 200 m contour (m)")
lines(shelf_vect, col = "red", lwd = 2)

#optional: convert to km
dist_shelf_km <- dist_shelf_m / 1000
names(dist_shelf_km) <- "dist_shelf_km"

plot(dist_shelf_km, main = "Distance to 200 m contour (km)")
lines(vect(shelf_200_5070), col = "red", lwd = 2)



#GEBCO
gebco_file <- "C:/github/Whale-SDM/data/GEBCO/gebco_2025_n45.0_s25.0_w-130.0_e-110.0.nc"

gebco <- rast(gebco_file)
gebco
names(gebco)
plot(gebco, main = "GEBCO raw bathymetry/elevation")

# ocean depth in meters, positive downward
bathy <- ifel(gebco < 0, -gebco, NA)
names(bathy) <- "bathy_m"

plot(bathy, main = "Bathymetry (m)")

# neighborhood rugosity via GEBCO
rugosity <- focal(
  bathy,
  w = matrix(1, 3, 3),
  fun = sd,
  na.rm = TRUE
)

names(rugosity) <- "rugosity_sd"

plot(rugosity, main = "Rugosity proxy (SD of bathymetry)")


#Regularization into metric
# bathymetry to metric CRS
bathy_5070 <- project(bathy, "EPSG:5070", method = "bilinear")
crs(bathy_5070)

# rugosity FROM the projected bathymetry
rugosity_5070 <- focal(
  bathy_5070,
  w = matrix(1, 3, 3),
  fun = sd,
  na.rm = TRUE
)
names(rugosity_5070) <- "rugosity_sd"
crs(rugosity_5070)

#Project the shelf contour vector 
shelf_200_5070 <- st_transform(shelf_200_ca, 5070)
st_crs(shelf_200_5070)

#resample bathy and rugosity to the distance raster since dist_shelf_km is metric template
bathy_static <- resample(bathy_5070, dist_shelf_km, method = "bilinear")
rugosity_static <- resample(rugosity_5070, dist_shelf_km, method = "bilinear")

compareGeom(bathy_static, rugosity_static, stopOnError = FALSE)
compareGeom(bathy_static, dist_shelf_km, stopOnError = FALSE)

#check extent
crs(bathy_static)
crs(rugosity_static)
crs(dist_shelf_km)

compareGeom(bathy_static, rugosity_static, stopOnError = FALSE)
compareGeom(bathy_static, dist_shelf_km, stopOnError = FALSE)
compareGeom(rugosity_static, dist_shelf_km, stopOnError = FALSE)

res(bathy_static)
res(rugosity_static)
res(dist_shelf_km)

ext(bathy_static)
ext(rugosity_static)
ext(dist_shelf_km)

#visual overlay check
par(mfrow = c(1, 3))
plot(bathy_static, main = "Bathymetry")
plot(rugosity_static, main = "Rugosity")
plot(dist_shelf_km, main = "Distance to shelf (km)")
nlyr(dist_shelf_km)


#write raster files
dir.create("C:/github/Whale-SDM/data/static", recursive = TRUE, showWarnings = FALSE)

writeRaster(bathy_static,
            "C:/github/Whale-SDM/data/static/bathy_static_8km_5070.tif",
            overwrite = TRUE)

writeRaster(rugosity_static,
            "C:/github/Whale-SDM/data/static/rugosity_static_8km_5070.tif",
            overwrite = TRUE)

writeRaster(dist_shelf_km,
            "C:/github/Whale-SDM/data/static/dist_shelf_km_8km_5070.tif",
            overwrite = TRUE)