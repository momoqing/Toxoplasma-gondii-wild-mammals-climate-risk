
packages <- c(
  "readxl", "dplyr", "tidyr", "tibble", "terra", "sf",
  "ggplot2", "RColorBrewer", "caret", "randomForest",
  "stringr", "scales"
)

need_install <- packages[!packages %in% rownames(installed.packages())]
if (length(need_install) > 0) {
  install.packages(need_install)
}

invisible(lapply(packages, library, character.only = TRUE))

set.seed(123)


ml_file <- "ml.xlsx"

path_current <- "D:/R/R/climate_data/wc2.1_10m_bio"

path_bcc <- "D:/R/R/climate_data/bcc"
path_ec  <- "D:/R/R/climate_data/EC"

out_dir <- "D:/R/R/Tgondii_future_projection_outputs_corrected"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rf_model_path   <- "Best_RF_Model.rds"
dummies_path    <- "dummies_rf.rds"
nzv_path        <- "nzv_rf.rds"
preProc_path    <- "preProc_rf.rds"
train_mode_path <- "train_mode_rf.rds"

n_points <- 20000

main_period <- "2041-2060"


rf_model   <- readRDS(rf_model_path)
dummies    <- readRDS(dummies_path)
nzv        <- readRDS(nzv_path)
preProc    <- readRDS(preProc_path)
train_mode <- readRDS(train_mode_path)

data_raw <- readxl::read_excel(ml_file, sheet = "Sheet1")

model_data <- data_raw %>%
  dplyr::rename(
  ) %>%
  dplyr::select(
    Prevalence,
    Year,
    Rainfall,
    Temperature,
    Longitude,
    Latitude,
    Total,
    Order,
    Family,
    Diet,
    Habitat,
    Sample_Type,
    Methods,
    Continent,
    Country
  ) %>%
  dplyr::filter(
    !is.na(Prevalence),
    !is.na(Longitude),
    !is.na(Latitude),
    Prevalence >= 0,
    Prevalence <= 1
  )

feature_vars <- colnames(model_data)[colnames(model_data) != "Prevalence"]

cat("\nTraining rainfall summary, expected unit = mm/day:\n")
print(summary(model_data$Rainfall))

cat("\nTraining temperature summary, expected unit = degree Celsius:\n")
print(summary(model_data$Temperature))

parse_ssp <- function(file) {
  token <- stringr::str_extract(tolower(basename(file)), "ssp(126|245|370|585)")
  dplyr::recode(
    token,
    "ssp126" = "SSP1-2.6",
    "ssp245" = "SSP2-4.5",
    "ssp370" = "SSP3-7.0",
    "ssp585" = "SSP5-8.5",
    .default = NA_character_
  )
}

parse_period <- function(file) {
  stringr::str_extract(basename(file), "20[0-9]{2}-20[0-9]{2}")
}

year_from_period <- function(period) {
  dplyr::case_when(
    period == "2021-2040" ~ 2030,
    period == "2041-2060" ~ 2050,
    TRUE ~ NA_real_
  )
}

get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

# Works whether train_mode is a named list, data frame, or named vector.
get_train_mode_value <- function(var_name) {
  if (is.list(train_mode) && var_name %in% names(train_mode)) {
    return(as.character(train_mode[[var_name]][1]))
  }
  if (is.data.frame(train_mode) && var_name %in% colnames(train_mode)) {
    return(as.character(train_mode[[var_name]][1]))
  }
  if (var_name %in% colnames(model_data)) {
    return(as.character(get_mode(model_data[[var_name]])))
  }
  NA_character_
}

fix_temperature_scale <- function(x) {
  # Some WorldClim temperature rasters are stored as degree Celsius × 10.
  # If values look too large, divide by 10.
  q95 <- suppressWarnings(stats::quantile(abs(x), 0.95, na.rm = TRUE))
  if (is.finite(q95) && q95 > 80) {
    message("Temperature appears to be scaled by 10. Dividing by 10.")
    return(x / 10)
  }
  x
}

fix_rainfall_scale <- function(x) {
  # WorldClim bio12 is annual precipitation in mm/year.
  # The machine-learning dataset uses rainfall in mm/day.
  # If values are clearly annual precipitation, divide by 365.
  med <- suppressWarnings(stats::median(x, na.rm = TRUE))
  if (is.finite(med) && med > 50) {
    message("Rainfall appears to be annual precipitation in mm/year. Dividing by 365.")
    return(x / 365)
  }
  x
}

extract_climate_at_points <- function(r, points_df, layer_index) {
  # points_df must contain Longitude and Latitude.
  pts <- terra::vect(
    points_df,
    geom = c("Longitude", "Latitude"),
    crs = "EPSG:4326"
  )

  vals <- terra::extract(r[[layer_index]], pts, ID = FALSE)
  out <- dplyr::bind_cols(points_df, as.data.frame(vals))

  colnames(out) <- c("Longitude", "Latitude", "Temperature", "Rainfall")

  out <- out %>%
    dplyr::mutate(
      Temperature = fix_temperature_scale(Temperature),
      Rainfall = fix_rainfall_scale(Rainfall)
    ) %>%
    dplyr::filter(!is.na(Temperature), !is.na(Rainfall))

  return(out)
}

make_prediction_data <- function(clim_df, year_value) {
  pred_df <- clim_df %>%
    dplyr::mutate(
      Year = year_value,
      Total = stats::median(model_data$Total, na.rm = TRUE),
      Order = get_train_mode_value("Order"),
      Family = get_train_mode_value("Family"),
      Diet = get_train_mode_value("Diet"),
      Habitat = get_train_mode_value("Habitat"),
      Sample_Type = get_train_mode_value("Sample_Type"),
      Methods = get_train_mode_value("Methods"),
      Continent = get_train_mode_value("Continent"),
      Country = get_train_mode_value("Country")
    ) %>%
    dplyr::select(dplyr::all_of(feature_vars))

  # Match factor levels stored by caret::dummyVars where possible.
  if (!is.null(dummies$lvls)) {
    for (v in names(dummies$lvls)) {
      if (v %in% colnames(pred_df)) {
        pred_df[[v]] <- factor(as.character(pred_df[[v]]), levels = dummies$lvls[[v]])
      }
    }
  } else {
    pred_df <- pred_df %>%
      dplyr::mutate(dplyr::across(where(is.character), as.factor))
  }

  pred_df
}

get_required_predictors <- function(model) {
  # caret::train object
  if (inherits(model, "train")) {
    if (!is.null(model$trainingData)) {
      return(setdiff(colnames(model$trainingData), ".outcome"))
    }
    if (!is.null(model$finalModel$importance)) {
      return(rownames(model$finalModel$importance))
    }
  }

  # randomForest object
  if (!is.null(model$importance)) {
    return(rownames(model$importance))
  }

  stop("Cannot determine predictor names from rf_model.")
}

get_preproc_predictors <- function(preProc) {
  cols <- unique(c(
    names(preProc$mean),
    names(preProc$std),
    names(preProc$ranges),
    names(preProc$median),
    names(preProc$bc),
    names(preProc$yj)
  ))

  cols <- cols[!is.na(cols)]

  if (length(cols) == 0) {
    return(NULL)
  }

  cols
}

align_columns <- function(x, required_cols, fill_value = 0) {
  x <- as.data.frame(x)

  missing_cols <- setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0) {
    for (nm in missing_cols) {
      x[[nm]] <- fill_value
    }
    message("Added missing columns with 0: ", paste(missing_cols, collapse = ", "))
  }

  extra_cols <- setdiff(colnames(x), required_cols)
  if (length(extra_cols) > 0) {
    x <- x[, setdiff(colnames(x), extra_cols), drop = FALSE]
  }

  x <- x[, required_cols, drop = FALSE]
  x
}

required_predictors <- get_required_predictors(rf_model)
preproc_predictors <- get_preproc_predictors(preProc)

predict_rf_prevalence <- function(pred_df) {
  # Some dummyVars objects were trained with a data frame containing Prevalence.
  # Add a temporary outcome column and remove it after dummy coding.
  pred_df$Prevalence <- 0

  dummy_mat <- predict(dummies, newdata = pred_df)
  dummy_df <- as.data.frame(dummy_mat)

  # Remove outcome if generated by dummyVars.
  dummy_df <- dummy_df[, !colnames(dummy_df) %in% "Prevalence", drop = FALSE]

  # Remove zero-variance or near-zero-variance variables as in training.
  if (length(nzv) > 0) {
    if (is.numeric(nzv)) {
      idx <- nzv[nzv <= ncol(dummy_df)]
      if (length(idx) > 0) {
        dummy_df <- dummy_df[, -idx, drop = FALSE]
      }
    } else {
      dummy_df <- dummy_df[, setdiff(colnames(dummy_df), nzv), drop = FALSE]
    }
  }

  # Align BEFORE preProc to avoid missing columns during preprocessing.
  if (!is.null(preproc_predictors)) {
    dummy_df <- align_columns(dummy_df, preproc_predictors, fill_value = 0)
  }

  x <- as.data.frame(predict(preProc, dummy_df))

  # Align again BEFORE randomForest prediction.
  x <- align_columns(x, required_predictors, fill_value = 0)

  pred <- as.numeric(predict(rf_model, newdata = x))

  # Keep predictions within prevalence range.
  pred <- pmin(pmax(pred, 0), 1)

  return(pred)
}


bio1_file <- file.path(path_current, "wc2.1_10m_bio_1.tif")
bio12_file <- file.path(path_current, "wc2.1_10m_bio_12.tif")

if (!file.exists(bio1_file)) stop("Cannot find baseline bio1 file: ", bio1_file)
if (!file.exists(bio12_file)) stop("Cannot find baseline bio12 file: ", bio12_file)

bio1_current <- terra::rast(bio1_file)
bio12_current <- terra::rast(bio12_file)

current_stack <- c(bio1_current, bio12_current)
names(current_stack) <- c("Temperature", "Rainfall")

spatial_points <- terra::spatSample(
  current_stack,
  size = n_points,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE,
  xy = TRUE
) %>%
  dplyr::transmute(Longitude = x, Latitude = y)

write.csv(
  spatial_points,
  file.path(out_dir, paste0("固定空间点_", n_points, ".csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


current_clim <- extract_climate_at_points(
  r = current_stack,
  points_df = spatial_points,
  layer_index = c(1, 2)
)

cat("\nCurrent baseline rainfall summary after correction, expected unit = mm/day:\n")
print(summary(current_clim$Rainfall))

cat("\nCurrent baseline temperature summary after correction:\n")
print(summary(current_clim$Temperature))

current_pred_df <- make_prediction_data(
  clim_df = current_clim,
  year_value = 2020
)

prev_current <- predict_rf_prevalence(current_pred_df)

current_result <- current_clim %>%
  dplyr::transmute(
    Longitude,
    Latitude,
    Current_Temperature = Temperature,
    Current_Rainfall = Rainfall,
    Current_Risk = prev_current
  )

write.csv(
  current_result,
  file.path(out_dir, paste0("弓形虫_当前风险基准_", n_points, "points_降水已转换.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

future_files <- dplyr::bind_rows(
  tibble::tibble(
    file = list.files(path_bcc, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE),
    Climate_Model = "BCC-CSM2-MR"
  ),
  tibble::tibble(
    file = list.files(path_ec, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE),
    Climate_Model = "EC-Earth3-Veg"
  )
) %>%
  dplyr::mutate(
    SSP = vapply(file, parse_ssp, character(1)),
    Time_Period = vapply(file, parse_period, character(1))
  ) %>%
  dplyr::filter(!is.na(SSP), !is.na(Time_Period)) %>%
  dplyr::mutate(
    SSP = factor(SSP, levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")),
    Time_Period = factor(Time_Period, levels = c("2021-2040", "2041-2060")),
    Climate_Model = factor(Climate_Model, levels = c("BCC-CSM2-MR", "EC-Earth3-Veg"))
  ) %>%
  dplyr::arrange(Climate_Model, SSP, Time_Period)

if (nrow(future_files) == 0) {
  stop("No future climate files were detected. Check path_bcc, path_ec, and file names.")
}

message("\nDetected future climate files:")
print(future_files %>% dplyr::select(file, Climate_Model, SSP, Time_Period))

if (nrow(future_files) != 16) {
  warning(
    "Expected 16 future climate files, but detected ",
    nrow(future_files),
    ". Please check whether all 2 GCM × 4 SSP × 2 period files are present."
  )
}

process_future_file <- function(file, climate_model, ssp, period) {
  message("Processing: ", basename(file), " | ", climate_model, " | ", ssp, " | ", period)

  r <- terra::rast(file)

  if (terra::nlyr(r) < 12) {
    stop(
      "Future climate file has fewer than 12 layers: ", file,
      "\nExpected a multilayer WorldClim bioclimatic TIFF with bio1 and bio12."
    )
  }

  future_clim <- extract_climate_at_points(
    r = r,
    points_df = spatial_points,
    layer_index = c(1, 12)
  )

  future_pred_df <- make_prediction_data(
    clim_df = future_clim,
    year_value = year_from_period(as.character(period))
  )

  prev_future <- predict_rf_prevalence(future_pred_df)

  future_clim %>%
    dplyr::transmute(
      Longitude,
      Latitude,
      Temperature,
      Rainfall,
      Prevalence = prev_future,
      Climate_Model = as.character(climate_model),
      SSP = as.character(ssp),
      Time_Period = as.character(period),
      File_Name = basename(file)
    )
}

pred_result <- dplyr::bind_rows(lapply(seq_len(nrow(future_files)), function(i) {
  process_future_file(
    file = future_files$file[i],
    climate_model = future_files$Climate_Model[i],
    ssp = future_files$SSP[i],
    period = future_files$Time_Period[i]
  )
})) %>%
  dplyr::mutate(
    SSP = factor(SSP, levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")),
    Time_Period = factor(Time_Period, levels = c("2021-2040", "2041-2060")),
    Climate_Model = factor(Climate_Model, levels = c("BCC-CSM2-MR", "EC-Earth3-Veg"))
  )

cat("\nFuture rainfall summary after correction, expected unit = mm/day:\n")
print(summary(pred_result$Rainfall))

cat("\nFuture temperature summary after correction:\n")
print(summary(pred_result$Temperature))

write.csv(
  pred_result,
  file.path(out_dir, "弓形虫_未来风险预测_16组合_含气候变量_降水已转换.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


pred_ensemble <- pred_result %>%
  dplyr::group_by(Longitude, Latitude, SSP, Time_Period) %>%
  dplyr::summarise(
    Temperature = mean(Temperature, na.rm = TRUE),
    Rainfall = mean(Rainfall, na.rm = TRUE),
    Prevalence = mean(Prevalence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SSP = factor(SSP, levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")),
    Time_Period = factor(Time_Period, levels = c("2021-2040", "2041-2060"))
  )

write.csv(
  pred_ensemble,
  file.path(out_dir, "弓形虫_未来风险预测_GCM平均_8组合_含气候变量_降水已转换.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

print(dim(pred_ensemble))
print(table(pred_ensemble$SSP, pred_ensemble$Time_Period))


eps <- 1e-6

risk_change_ensemble <- pred_ensemble %>%
  dplyr::left_join(current_result, by = c("Longitude", "Latitude")) %>%
  dplyr::mutate(
    Absolute_Change = Prevalence - Current_Risk,
    Relative_Change_Pct = (Absolute_Change / pmax(Current_Risk, eps)) * 100,
    Temperature_Change = Temperature - Current_Temperature,
    Rainfall_Change = Rainfall - Current_Rainfall
  )

write.csv(
  risk_change_ensemble,
  file.path(out_dir, "弓形虫_风险变化_GCM平均_8组合_降水已转换.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

risk_change_16 <- pred_result %>%
  dplyr::left_join(current_result, by = c("Longitude", "Latitude")) %>%
  dplyr::mutate(
    Absolute_Change = Prevalence - Current_Risk,
    Relative_Change_Pct = (Absolute_Change / pmax(Current_Risk, eps)) * 100,
    Temperature_Change = Temperature - Current_Temperature,
    Rainfall_Change = Rainfall - Current_Rainfall
  )

write.csv(
  risk_change_16,
  file.path(out_dir, "弓形虫_风险变化_16组合_含气候变量_降水已转换.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

risk_summary <- pred_ensemble %>%
  dplyr::group_by(SSP, Time_Period) %>%
  dplyr::summarise(
    mean_risk = mean(Prevalence, na.rm = TRUE),
    sd_risk = sd(Prevalence, na.rm = TRUE),
    q05 = stats::quantile(Prevalence, 0.05, na.rm = TRUE),
    q50 = stats::quantile(Prevalence, 0.50, na.rm = TRUE),
    q95 = stats::quantile(Prevalence, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

climate_summary_16 <- pred_result %>%
  dplyr::group_by(Climate_Model, SSP, Time_Period) %>%
  dplyr::summarise(
    mean_temp = mean(Temperature, na.rm = TRUE),
    sd_temp = sd(Temperature, na.rm = TRUE),
    mean_rain = mean(Rainfall, na.rm = TRUE),
    sd_rain = sd(Rainfall, na.rm = TRUE),
    mean_risk = mean(Prevalence, na.rm = TRUE),
    sd_risk = sd(Prevalence, na.rm = TRUE),
    .groups = "drop"
  )

climate_summary_ensemble <- pred_ensemble %>%
  dplyr::group_by(SSP, Time_Period) %>%
  dplyr::summarise(
    mean_temp = mean(Temperature, na.rm = TRUE),
    sd_temp = sd(Temperature, na.rm = TRUE),
    mean_rain = mean(Rainfall, na.rm = TRUE),
    sd_rain = sd(Rainfall, na.rm = TRUE),
    mean_risk = mean(Prevalence, na.rm = TRUE),
    sd_risk = sd(Prevalence, na.rm = TRUE),
    .groups = "drop"
  )

risk_change_summary <- risk_change_ensemble %>%
  dplyr::group_by(SSP, Time_Period) %>%
  dplyr::summarise(
    mean_abs_change = mean(Absolute_Change, na.rm = TRUE),
    sd_abs_change = sd(Absolute_Change, na.rm = TRUE),
    q05_abs_change = stats::quantile(Absolute_Change, 0.05, na.rm = TRUE),
    q50_abs_change = stats::quantile(Absolute_Change, 0.50, na.rm = TRUE),
    q95_abs_change = stats::quantile(Absolute_Change, 0.95, na.rm = TRUE),
    mean_rel_change_pct = mean(Relative_Change_Pct, na.rm = TRUE),
    sd_rel_change_pct = sd(Relative_Change_Pct, na.rm = TRUE),
    q05_rel_change_pct = stats::quantile(Relative_Change_Pct, 0.05, na.rm = TRUE),
    q50_rel_change_pct = stats::quantile(Relative_Change_Pct, 0.50, na.rm = TRUE),
    q95_rel_change_pct = stats::quantile(Relative_Change_Pct, 0.95, na.rm = TRUE),
    positive_fraction = mean(Relative_Change_Pct > 0, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(risk_summary, file.path(out_dir, "诊断_预测风险_summary_降水已转换.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(climate_summary_16, file.path(out_dir, "诊断_气候输入_16组合_summary_降水已转换.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(climate_summary_ensemble, file.path(out_dir, "诊断_气候输入_GCM平均_summary_降水已转换.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(risk_change_summary, file.path(out_dir, "诊断_风险变化_summary_降水已转换.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


plot_prediction_map <- function(df, facet_formula, title_text, out_file,
                                width = 14, height = 8,
                                color_limits = c(0.24, 0.38),
                                color_breaks = seq(0.24, 0.38, by = 0.02)) {
  
  df_sf <- sf::st_as_sf(
    df,
    coords = c("Longitude", "Latitude"),
    crs = 4326
  )
  
  p <- ggplot2::ggplot(df_sf)
  
  p <- p +
    ggplot2::geom_sf(
      ggplot2::aes(color = Prevalence),
      size = 0.55,
      alpha = 0.9
    )
  
  p <- p +
    ggplot2::scale_color_gradientn(
      colors = c("#FFF7BC", "#FEC44F", "#FE9929", "#D95F0E", "#993404"),
      name = "Predicted\nrisk",
      limits = color_limits,
      breaks = color_breaks,
      labels = scales::number_format(accuracy = 0.01),
      oob = scales::oob_squish
    )
  
  p <- p + facet_formula
  
  p <- p +
    ggplot2::labs(title = title_text) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F0F0F0"),
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(
    filename = out_file,
    plot = p,
    width = width,
    height = height,
    dpi = 600,
    device = "pdf",
    useDingbats = FALSE
  )
  
  return(p)
}
plot_relative_change_map <- function(df, facet_formula, title_text, out_file,
                                     width = 14, height = 8,
                                     color_limits = NULL) {
  if (is.null(color_limits)) {
    lim <- stats::quantile(abs(df$Relative_Change_Pct), 0.98, na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) lim <- 10
    color_limits <- c(-lim, lim)
  }

  df_sf <- sf::st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)

  p <- ggplot2::ggplot(df_sf) +
    ggplot2::geom_sf(ggplot2::aes(color = Relative_Change_Pct), size = 0.55, alpha = 0.9) +
    ggplot2::scale_color_gradient2(
      low = "#2C7BB6",
      mid = "#FFFFBF",
      high = "#D7191C",
      midpoint = 0,
      limits = color_limits,
      oob = scales::oob_squish,
      name = "Relative\nchange (%)"
    ) +
    facet_formula +
    ggplot2::labs(title = title_text) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F0F0F0"),
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

  ggplot2::ggsave(out_file, p, width = width, height = height, dpi = 600,
                  device = "pdf", useDingbats = FALSE)
  return(p)
}

main_df <- pred_ensemble %>%
  dplyr::filter(Time_Period == main_period)

plot_prediction_map(
  df = main_df,
  facet_formula = ggplot2::facet_wrap(~ SSP, ncol = 2),
  title_text = paste0("Projected T. gondii infection risk, GCM ensemble mean, ", main_period),
  out_file = file.path(out_dir, paste0("Figure4_Absolute_risk_GCMensemble_", main_period, "_4SSP_降水已转换.pdf")),
  width = 14,
  height = 9,
  color_limits = c(0.20, 0.40),
  color_breaks = seq(0.20, 0.40, by = 0.05)
)

plot_prediction_map(
  df = pred_ensemble,
  facet_formula = ggplot2::facet_grid(SSP ~ Time_Period),
  title_text = "Projected T. gondii infection risk, GCM ensemble mean",
  out_file = file.path(out_dir, "FigureS5_Absolute_risk_GCMensemble_4SSP_2periods_降水已转换.pdf"),
  width = 16,
  height = 14,
  color_limits = c(0.20, 0.40),
  color_breaks = seq(0.20, 0.40, by = 0.05)
)

rel_lim <- stats::quantile(abs(risk_change_ensemble$Relative_Change_Pct), 0.98, na.rm = TRUE)
if (!is.finite(rel_lim) || rel_lim == 0) rel_lim <- 10
rel_limits <- c(-rel_lim, rel_lim)

main_change_df <- risk_change_ensemble %>%
  dplyr::filter(Time_Period == main_period)

plot_relative_change_map(
  df = main_change_df,
  facet_formula = ggplot2::facet_wrap(~ SSP, ncol = 2),
  title_text = paste0("Projected relative change in T. gondii infection risk, ", main_period),
  out_file = file.path(out_dir, paste0("Figure5_Relative_change_GCMensemble_", main_period, "_4SSP_降水已转换.pdf")),
  width = 14,
  height = 9,
  color_limits = rel_limits
)

plot_relative_change_map(
  df = risk_change_ensemble,
  facet_formula = ggplot2::facet_grid(SSP ~ Time_Period),
  title_text = "Projected relative change in T. gondii infection risk, GCM ensemble mean",
  out_file = file.path(out_dir, "FigureS6_Relative_change_GCMensemble_4SSP_2periods_降水已转换.pdf"),
  width = 16,
  height = 14,
  color_limits = rel_limits
)
