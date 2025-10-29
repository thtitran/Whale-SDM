library(dplyr)
library(lubridate)
# ==CLEANING/PROCESSING/SIMPLIFYING DATA TO MATCH ORIGINAL FORMAT==
#1 Load  raw Argos data
raw_df <- read.csv("C:/github/Whale-SDM/tracks/254025/254025-RawArgos.csv")

#2 Remove rows with missing PassDate/PassTime
clean <- raw_df %>%
  filter(!is.na(PassDate), !is.na(PassTime))

#3 Combine and parse date/time
clean <- clean %>%
  mutate(
    combined = paste(PassDate, PassTime),
    date = suppressWarnings(dmy_hms(combined, tz = "UTC"))
  )

#4 Drop any rows that still failed to parse
clean <- clean %>%
  filter(!is.na(date)) %>%
  transmute(
    id  = PTT,
    date,
    lc  = toupper(trimws(Class)),
    lon = Longitude,
    lat = Latitude
  ) %>%
  arrange(id, date)

#5 Clean out NA's
clean <- clean %>%
  filter(!is.na(lon), !is.na(lat), !is.na(date))

#6 Save the cleaned track
write.csv(clean, "C:/github/Whale-SDM/tracks/254025/254025 track for SSM.csv", row.names = FALSE)

cat("✅ Cleaned file saved as '2540225track for SSM.csv' \n")



# ==MODEL PROCESSING==
library(aniMotum)
# Read the author-format CSV (UTC datetimes) ----
trk <- read.csv("C:/github/Whale-SDM/tracks/254025/254025 track for SSM.csv") |>
  mutate(date = ymd_hms(date, tz = "UTC")) |>
  arrange(id, date)

# 1 but optional: split long gaps to avoid weird interpolation across outages ----
gap_days <- 7  # good baseline gap day but use 5–10 depending
trk <- trk |>
  group_by(id) |>
  arrange(date, .by_group = TRUE) |>
  mutate(
    g = as.numeric(difftime(date, dplyr::lag(date), units = "days")),
    seg = cumsum(if_else(is.na(g) | g <= gap_days, 0L, 1L)),
    id = paste0(id, "_", seg)
  ) |>
  ungroup() |>
  select(id, date, lc, lon, lat)

# 2 Drop tiny segments (<5 pts) that break the smoother
trk <- trk |>
  group_by(id) |>
  filter(n() >= 5) |>
  ungroup()

# 3 Fit SSM
# time.step ~ median sampling interval; Argos-only often 6–24h. Start with 6h.
fit <- fit_ssm(
  trk,                #  cleaned 5-col data (optionally segmented)
  model     = "crw",  # correlated random walk for mobile animals
  time.step = 6,      # hours; regularization interval
  vmax      = 10,     # m/s; permissive cap for large whales
  spdf      = FALSE   # return data.frame/tibble
)

# 4 Regularized positions for analysis/mapping
#Produces error filtered + time-regularized locations for analysis.
pred <- grab(fit, what = "predicted")
write.csv(pred, "tracks_pred_6h.csv", row.names = FALSE)

# 5 Move-persistence (γ) + behavior class
#fit_mpm() computes gamma on  SSM regularized points
#g < 0.8 is the paper’s threshold (low persistence =~ area-restricted search vs. high persistence ≈ transit). Should be adjusted?
mp <- fit_mpm(fit, what = "predicted")
mpred <- grab(mp, what = "fitted")
# Handle differences in aniMotum versions (g vs gamma, list columns, etc.)
if (!"gamma" %in% names(mpred) && "g" %in% names(mpred)) {
  mpred <- dplyr::rename(mpred, gamma = g)
}
mpred <- mpred %>%
  dplyr::mutate(
    gamma = {
      x <- gamma
      if (is.list(x)) x <- unlist(x, use.names = FALSE)
      if (is.matrix(x)) x <- as.numeric(x)
      as.numeric(x)
    },
    behav = dplyr::if_else(gamma < 0.8, "low_persist", "high_persist")
  )
# Write output for downstream analysis
write.csv(mpred, "tracks_pred_6h_mpm.csv", row.names = FALSE)


