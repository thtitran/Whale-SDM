library(mgcv) #GAMMs
library(lubridate) #datetimes

# 1. Load and prep data

# Response variable:
# PresAbs = 1 means a used/presence point
# PresAbs = 0 means an available/pseudo-absence point
tracks <- read.csv("C:/github/Whale-SDM/Output/GAMM_env/254026_pres_pa_env_GAMMready.csv")

# Show  first 6 rows to visually confirm file loaded correctly
head(tracks)

# Show the structure of the dataframe:
str(tracks)

# Convert the whale ID column to a factor for categorical variables
tracks$id <- as.factor(tracks$id)

# Convert dTime into a proper date-time object (POSIXct)
tracks$dTime <- as.POSIXct(tracks$dTime, tz = "UTC")

# Convert the date column into a proper Date object
tracks$date <- as.Date(tracks$date)

# Create a new month column from the date column
# 1 = January, 2 = February, ..., 12 = December
tracks$month <- month(tracks$date)

# Remove rows with any missing values for convienece
tracks <- na.omit(tracks)

# 2.Fit candidate GAMMs ####
# Individuals (variable name 'ptt') are nested as a random effect

# 2. Fit candidate GAMMs ----

# Model 1:
# temperature, sea surface height, eastward current, northward current,
# salinity, bathymetry, rugosity, distance to shelf
gam.mod1 <- mgcv::gam(
  PresAbs ~
    s(thetao, bs = "ts") +
    s(zos, bs = "ts") +
    s(uo, bs = "ts") +
    s(vo, bs = "ts") +
    s(so, bs = "ts") +
    s(bathy_km, bs = "ts") +
    s(rugosity, bs = "ts") +
    s(dist_shelf_km, bs = "ts"),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 2:
# Model 1 + spatial surface
gam.mod2 <- mgcv::gam(
  PresAbs ~
    s(thetao, bs = "ts") +
    s(zos, bs = "ts") +
    s(uo, bs = "ts") +
    s(vo, bs = "ts") +
    s(so, bs = "ts") +
    s(bathy_km, bs = "ts") +
    s(rugosity, bs = "ts") +
    s(dist_shelf_km, bs = "ts") +
    te(long, lat, bs = c("ts", "ts")),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 3:
# static habitat variables + dynamic variables, with temperature-by-space interaction
gam.mod3 <- mgcv::gam(
  PresAbs ~
    s(zos, bs = "ts") +
    s(uo, bs = "ts") +
    s(vo, bs = "ts") +
    s(so, bs = "ts") +
    s(bathy_km, bs = "ts") +
    s(rugosity, bs = "ts") +
    s(dist_shelf_km, bs = "ts") +
    te(thetao, lat, bs = c("ts", "ts")),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 4:
# simple baseline: temperature + bathymetry + rugosity + shelf distance
gam.mod4 <- mgcv::gam(
  PresAbs ~
    s(thetao, bs = "ts") +
    s(bathy_km, bs = "ts") +
    s(rugosity, bs = "ts") +
    s(dist_shelf_km, bs = "ts"),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 5:
# current speed + temperature + sea surface height + salinity + static habitat
gam.mod5 <- mgcv::gam(
  PresAbs ~
    s(current_speed, bs = "ts") +
    s(thetao, bs = "ts") +
    s(zos, bs = "ts") +
    s(so, bs = "ts") +
    s(bathy_km, bs = "ts") +
    s(rugosity, bs = "ts") +
    s(dist_shelf_km, bs = "ts"),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 6:
# temperature + sea surface height + bathymetry + rugosity + shelf distance + eastward current
gam.mod6 <- mgcv::gam(
  PresAbs ~
    s(thetao, bs = "ts") +
    s(zos, bs = "ts") +
    s(bathy_km, bs = "ts") +
    s(rugosity, bs = "ts") +
    s(dist_shelf_km, bs = "ts") +
    s(uo, bs = "ts"),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 7:
# temperature + sea surface height + bathymetry + rugosity + shelf distance + northward current
gam.mod7 <- mgcv::gam(
  PresAbs ~
    s(thetao, bs = "ts") +
    s(zos, bs = "ts") +
    s(bathy_km, bs = "ts") +
    s(rugosity, bs = "ts") +
    s(dist_shelf_km, bs = "ts") +
    s(vo, bs = "ts"),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 8:
# temperature + sea surface height + bathymetry-by-rugosity interaction + shelf distance
gam.mod8 <- mgcv::gam(
  PresAbs ~
    s(thetao, bs = "ts") +
    s(zos, bs = "ts") +
    te(bathy_km, rugosity, bs = c("ts", "ts")) +
    s(dist_shelf_km, bs = "ts"),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 9:
# salinity + sea surface height + bathymetry-by-rugosity interaction + shelf distance + current speed
gam.mod9 <- mgcv::gam(
  PresAbs ~
    s(so, bs = "ts") +
    s(zos, bs = "ts") +
    te(bathy_km, rugosity, bs = c("ts", "ts")) +
    s(dist_shelf_km, bs = "ts") +
    s(current_speed, bs = "ts"),
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

# Model 10:
# temperature + sea surface height + bathymetry + rugosity + shelf distance + month
#treat month as a simple factor because of error trying to fit a curve with more flexibility than the variable actually contains
tracks$month_f <- as.factor(tracks$month)

gam.mod10 <- mgcv::gam(
  PresAbs ~
    s(thetao, bs = "ts") +
    s(zos, bs = "ts") +
    s(bathy_km, bs = "ts") +
    s(rugosity, bs = "ts") +
    s(dist_shelf_km, bs = "ts") +
    month_f,
  family = binomial,
  data = tracks,
  method = "REML",
  select = TRUE
)

plot(gam.mod4, pages = 1, shade = TRUE, residuals = TRUE)

