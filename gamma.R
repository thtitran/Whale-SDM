suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(aniMotum)
  library(purrr)
})

# ==== PATHS ====
input_dir <- "C:/github/Whale-SDM/combinedtracks"        # folder containing your CSVs
out_dir   <- "C:/github/Whale-SDM/tracks/outputnewgamma" # results will be written here
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==== PARAMETERS ====
gap_days        <- 7     # split gaps > this many days; set NULL to disable splitting
min_points      <- 5     # drop segments with fewer points (SSM needs >= 5)
time_step_hours <- 6     # SSM regularization interval (hours); tune to median fix rate
vmax_ms         <- 10    # plausible max speed (m/s); tighten if many points flagged
gamma_thresh    <- 0.34   # behavior threshold (gamma < thresh => low_persist)

# ==== HELPERS ====
parse_date_utc <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  d <- suppressWarnings(ymd_hms(x, tz = "UTC"))
  d[is.na(d)] <- suppressWarnings(dmy_hms(x[is.na(d)], tz = "UTC"))
  d
}

segment_and_filter <- function(df, gap_split_days, min_n) {
  if (is.null(gap_split_days)) {
    out <- df
  } else {
    out <- df %>%
      group_by(id) %>%
      arrange(date, .by_group = TRUE) %>%
      mutate(
        gap = as.numeric(difftime(date, dplyr::lag(date), units = "days")),
        seg = cumsum(if_else(is.na(gap) | gap <= gap_split_days, 0L, 1L)),
        id  = paste0(id, "_", seg)
      ) %>%
      ungroup() %>%
      select(id, date, lc, lon, lat)
  }
  out %>% group_by(id) %>% filter(n() >= min_n) %>% ungroup()
}

process_one <- function(path,
                        time_step = time_step_hours,
                        vmax = vmax_ms,
                        gap_split_days = gap_days,
                        min_n = min_points,
                        g_thresh = gamma_thresh) {
  
  message("\n---- Processing: ", basename(path), " ----")
  base <- tools::file_path_sans_ext(basename(path))
  
  # Skip our own outputs if they live in input_dir
  if (grepl("_pred_\\d+h(.*)\\.csv$|_mpm\\.csv$|combined_(pred|mpm)\\.csv$", base)) {
    message("Skipping derived/output file: ", basename(path))
    return(NULL)
  }
  
  # 1) Read & coerce to author format if needed
  df <- read.csv(path, stringsAsFactors = FALSE)
  
  if (all(c("PassDate","PassTime") %in% names(df))) {
    # Raw Argos -> author format
    df <- df %>%
      mutate(date = parse_date_utc(paste(PassDate, PassTime))) %>%
      transmute(
        id  = as.character(PTT),
        date,
        lc  = toupper(trimws(Class)),
        lon = as.numeric(Longitude),
        lat = as.numeric(Latitude)
      )
  } else {
    # Author format assumed
    stopifnot(all(c("id","date","lc","lon","lat") %in% names(df)))
    df <- df %>%
      mutate(
        id   = as.character(id),
        date = parse_date_utc(date),
        lc   = toupper(trimws(lc)),
        lon  = as.numeric(lon),
        lat  = as.numeric(lat)
      )
  }
  
  # 2) Basic clean
  df <- df %>%
    filter(!is.na(date), !is.na(lon), !is.na(lat)) %>%
    arrange(id, date) %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon))
  
  if (nrow(df) < min_n) {
    warning("Skipping (too few total points): ", basename(path))
    return(NULL)
  }
  
  # 3) Segment & drop tiny — attempt A: with gap splitting
  dfA <- segment_and_filter(df, gap_split_days, min_n)
  
  # Attempt B: WITHOUT gap splitting (fallback)
  dfB <- df
  if (nrow(dfA) < min_n) {
    dfB <- segment_and_filter(df, gap_split_days = NULL, min_n)
  }
  
  df_use <- if (nrow(dfA) >= min_n) dfA else dfB
  if (nrow(df_use) < min_n) {
    warning("Skipping (no segment has ≥ ", min_n, " points): ", basename(path))
    return(NULL)
  }
  
  kept_ids <- df_use %>% count(id, name = "n_pts") %>% arrange(desc(n_pts))
  message("Fitting on ", nrow(kept_ids), " segment(s): ",
          paste0(kept_ids$id, "(", kept_ids$n_pts, ")", collapse = ", "))
  
  # 4) Fit SSM
  fit <- fit_ssm(
    df_use,
    model     = "crw",
    time.step = time_step,
    vmax      = vmax,
    spdf      = FALSE
  )
  
  # 5) Predicted (regularized) positions
  pred <- grab(fit, what = "predicted") %>%
    mutate(source_file = basename(path))
  
  pred_path <- file.path(out_dir, sprintf("%s_pred_%dh.csv", base, time_step))
  write.csv(pred, pred_path, row.names = FALSE)
  
  # 6) Move persistence model
  mp    <- fit_mpm(fit, what = "predicted")
  mpred <- grab(mp, what = "fitted") %>% mutate(source_file = basename(path))
  
  # Normalize gamma column (g vs gamma; list/matrix)
  if (!"gamma" %in% names(mpred) && "g" %in% names(mpred)) {
    mpred <- dplyr::rename(mpred, gamma = g)
  }
  mpred <- mpred %>%
    mutate(
      gamma = {
        x <- gamma
        if (is.list(x))   x <- unlist(x, use.names = FALSE)
        if (is.matrix(x)) x <- as.numeric(x)
        as.numeric(x)
      },
      behav = if_else(gamma < g_thresh, "low_persist", "high_persist")
    )
  
  mpm_path <- file.path(out_dir, sprintf("%s_pred_%dh_mpm.csv", base, time_step))
  write.csv(mpred, mpm_path, row.names = FALSE)
  
  # 7) Diagnostics PDF (non-interactive)
  pdf(file.path(out_dir, sprintf("%s_diagnostics.pdf", base)), width = 8, height = 6)
  op <- par(ask = FALSE); on.exit(par(op), add = TRUE)
  plot(fit)
  plot(mp)
  dev.off()
  
  list(pred = pred, mpm = mpred)
}

# ==== RUN OVER THE FOLDER ====
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

results <- lapply(files, function(f) {
  tryCatch(process_one(f),
           error = function(e) { warning("Failed on ", basename(f), ": ", e$message); NULL })
})

# keep only successful fits
ok <- purrr::compact(results)

if (length(ok) > 0) {
  pred_all <- dplyr::bind_rows(purrr::map(ok, "pred"))
  mpm_all  <- dplyr::bind_rows(purrr::map(ok, "mpm"))
  write.csv(pred_all, file.path(out_dir, "combined_pred.csv"), row.names = FALSE)
  write.csv(mpm_all,  file.path(out_dir, "combined_mpm.csv"),  row.names = FALSE)
  message("Combined outputs written to: ", out_dir)
} else {
  message("No tracks produced valid segments (all skipped). ",
          "Try gap_days = NULL or larger time_step_hours (e.g., 12–24).")
}
