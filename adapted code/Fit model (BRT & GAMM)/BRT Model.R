# Load libraries ----
library(dismo)
library(gbm)
library(lubridate)

# 1. Load and prep data ----
# Presence + CRW pseudo-absence dataset for whale 254026
# Response variable: PresAbs (1 = used, 0 = available)

tracks <- read.csv("C:/github/Whale-SDM/Output/GAMM_env/254026_pres_pa_env_GAMMready.csv")

head(tracks)
str(tracks)

#Convert ID to factor
tracks$id <- as.factor(tracks$id)

#Parse date columns
tracks$dTime <- as.POSIXct(tracks$dTime, tz = "UTC")
tracks$date  <- as.Date(tracks$date)
tracks$month <- month(tracks$date)

#Remove rows with missing values for a first-pass model
tracks_brtvars <- tracks[, c("PresAbs", "thetao", "zos", "so", "bathy_km", "rugosity", "dist_shelf_km")]
tracks_brtvars <- na.omit(tracks_brtvars)

table(tracks_brtvars$PresAbs)


#Temporarily thin PAs for BRT testing
#split the full dataset into:
# - pres = actual whale tracks
# - pa   = pseudo-absence points
pres <- subset(tracks, PresAbs == 1)
pa   <- subset(tracks, PresAbs == 0)

#check how many of each
table(tracks$PresAbs)

# Set a random seed so you get the same sampled pseudo-absences each time (reproducible)
set.seed(123)

# Chooses a pseudo-absence:presence ratio for this test
# Here, 5 means: keep 5 PA for every 1 presence
pa_ratio <- 5

# Randomly sample subset of PA
#size = number of presences * chosen ratio
pa_sub <- pa[sample(nrow(pa), size = pa_ratio * nrow(pres)), ]

# Recombine the full set of presences with the sampled pseudo-absences
tracks_brt <- rbind(pres, pa_sub)

# Check that the class balance now look reasonable
table(tracks_brt$PresAbs)
prop.table(table(tracks_brt$PresAbs))

table(tracks_brt$PresAbs)


# Fit first BRT
brt.mod1 <- dismo::gbm.step(
  data = tracks_brt,
  gbm.x = c("thetao", "zos", "so", "bathy_km", "rugosity", "dist_shelf_km"),
  gbm.y = "PresAbs",
  family = "bernoulli",
  tree.complexity = 2,
  learning.rate = 0.005,
  bag.fraction = 0.5
)

brt.mod2 <- dismo::gbm.step(
  data = tracks_brt,
  gbm.x = c("thetao", "zos", "so", "bathy_km", "rugosity", "dist_shelf_km"),
  gbm.y = "PresAbs",
  family = "bernoulli",
  tree.complexity = 2,
  learning.rate = 0.01,
  bag.fraction = 0.5
)


# Inspect variable importance
summary(brt.mod1)
summary(brt.mod2)


# Function to calculate percent deviance explained
dev_eval <- function(model_object) {
  null <- model_object$self.statistics$mean.null
  res  <- model_object$self.statistics$mean.resid
  dev  <- ((null - res) / null) * 100
  return(dev)
}

# Deviance explained
dev_eval(brt.mod1)

# Save variable contribution table
summary_res <- summary(brt.mod1)
write.csv(summary_res,
          file = "C:/github/Whale-SDM/Output/BRT/summary_results_brt_mod1.csv",
          row.names = FALSE)

# Save deviance explained
summary_dev <- dev_eval(brt.mod1)
write.csv(data.frame(deviance_explained = summary_dev),
          file = "C:/github/Whale-SDM/Output/BRT/summary_dev_brt_mod1.csv",
          row.names = FALSE)

# Plot partial dependence curves
png("C:/github/Whale-SDM/Output/BRT/partial_curves_brt_mod1.png",
    width = 12, height = 8, units = "in", res = 300, bg = "white")
gbm.plot(brt.mod1, smooth = TRUE, write.title = TRUE)
dev.off()

# Save model object
saveRDS(brt.mod1, "C:/github/Whale-SDM/Output/BRT/brt_mod1.rds")