
packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "sf",
  "RANN", "scales", "tibble"
)

need_install <- packages[!packages %in% rownames(installed.packages())]
if (length(need_install) > 0) {
  install.packages(need_install)
}

invisible(lapply(packages, library, character.only = TRUE))

set.seed(123)


ml_file <- "ml.xlsx"

out_dir <- "D:/R/R/Tgondii_future_projection_outputs_corrected"

current_file <- file.path(
  out_dir,
  "弓形虫_当前风险基准_20000points_降水已转换.csv"
)

future_file <- file.path(
  out_dir,
  "弓形虫_未来风险预测_GCM平均_8组合_含气候变量_降水已转换.csv"
)

change_file <- file.path(
  out_dir,
  "弓形虫_风险变化_GCM平均_8组合_降水已转换.csv"
)

novelty_out_dir <- file.path(out_dir, "Climate_novelty_analysis_corrected")
dir.create(novelty_out_dir, showWarnings = FALSE, recursive = TRUE)

core_lower <- 0.025
core_upper <- 0.975
mv_threshold_prob <- 0.95
mv_threshold_prob_strict <- 0.975

ml_raw <- readxl::read_excel(ml_file, sheet = "Sheet1")

training_climate <- ml_raw %>%
  dplyr::rename(

  ) %>%
  dplyr::select(Temperature, Rainfall) %>%
  dplyr::filter(
    !is.na(Temperature),
    !is.na(Rainfall)
  )

cat("\nTraining climate summary:\n")
print(summary(training_climate))

if (median(training_climate$Rainfall, na.rm = TRUE) > 50) {
  stop(
    "Training rainfall appears to be annual precipitation, not mm/day. ",
    "Please check ml.xlsx."
  )
}

train_range <- training_climate %>%
  dplyr::summarise(
    temp_min = min(Temperature, na.rm = TRUE),
    temp_q025 = stats::quantile(Temperature, core_lower, na.rm = TRUE),
    temp_q975 = stats::quantile(Temperature, core_upper, na.rm = TRUE),
    temp_max = max(Temperature, na.rm = TRUE),
    rain_min = min(Rainfall, na.rm = TRUE),
    rain_q025 = stats::quantile(Rainfall, core_lower, na.rm = TRUE),
    rain_q975 = stats::quantile(Rainfall, core_upper, na.rm = TRUE),
    rain_max = max(Rainfall, na.rm = TRUE)
  )

cat("\nTraining climate domain:\n")
print(train_range)

write.csv(
  train_range,
  file.path(novelty_out_dir, "Climate_training_domain_summary_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


if (!file.exists(current_file)) stop("Cannot find current_file: ", current_file)
if (!file.exists(future_file))  stop("Cannot find future_file: ", future_file)
if (!file.exists(change_file))  stop("Cannot find change_file: ", change_file)

current_result <- read.csv(current_file, fileEncoding = "UTF-8")
future_result  <- read.csv(future_file,  fileEncoding = "UTF-8")
risk_change    <- read.csv(change_file,  fileEncoding = "UTF-8")

required_current_cols <- c("Longitude", "Latitude", "Current_Temperature",
                           "Current_Rainfall", "Current_Risk")
required_future_cols <- c("Longitude", "Latitude", "SSP", "Time_Period",
                          "Temperature", "Rainfall", "Prevalence")
required_change_cols <- c("Longitude", "Latitude", "SSP", "Time_Period",
                          "Relative_Change_Pct", "Absolute_Change")

missing_current <- setdiff(required_current_cols, colnames(current_result))
missing_future  <- setdiff(required_future_cols, colnames(future_result))
missing_change  <- setdiff(required_change_cols, colnames(risk_change))

if (length(missing_current) > 0) stop("current_result missing columns: ", paste(missing_current, collapse = ", "))
if (length(missing_future) > 0)  stop("future_result missing columns: ", paste(missing_future, collapse = ", "))
if (length(missing_change) > 0)  stop("risk_change missing columns: ", paste(missing_change, collapse = ", "))

cat("\nCurrent projection rainfall summary, expected mm/day:\n")
print(summary(current_result$Current_Rainfall))

cat("\nFuture projection rainfall summary, expected mm/day:\n")
print(summary(future_result$Rainfall))

if (median(current_result$Current_Rainfall, na.rm = TRUE) > 50 ||
    median(future_result$Rainfall, na.rm = TRUE) > 50) {
  stop(
    "Projection rainfall still appears to be annual precipitation in mm/year. ",
    "Please rerun projection code with bio12 divided by 365."
  )
}

current_climate <- current_result %>%
  dplyr::transmute(
    Longitude,
    Latitude,
    SSP = "Current",
    Time_Period = "Current",
    Projection_Type = "Current baseline",
    Temperature = Current_Temperature,
    Rainfall = Current_Rainfall
  )

future_climate <- future_result %>%
  dplyr::transmute(
    Longitude,
    Latitude,
    SSP = as.character(SSP),
    Time_Period = as.character(Time_Period),
    Projection_Type = "Future ensemble",
    Temperature,
    Rainfall
  )

all_projection_climate <- dplyr::bind_rows(
  current_climate,
  future_climate
) %>%
  dplyr::filter(
    !is.na(Temperature),
    !is.na(Rainfall)
  )


temp_min  <- train_range$temp_min[1]
temp_max  <- train_range$temp_max[1]
rain_min  <- train_range$rain_min[1]
rain_max  <- train_range$rain_max[1]

temp_q025 <- train_range$temp_q025[1]
temp_q975 <- train_range$temp_q975[1]
rain_q025 <- train_range$rain_q025[1]
rain_q975 <- train_range$rain_q975[1]

train_center <- colMeans(training_climate[, c("Temperature", "Rainfall")], na.rm = TRUE)
train_scale  <- apply(training_climate[, c("Temperature", "Rainfall")], 2, sd, na.rm = TRUE)

if (any(!is.finite(train_scale)) || any(train_scale == 0)) {
  stop("Training climate scale contains zero or non-finite values.")
}

training_scaled <- scale(
  training_climate[, c("Temperature", "Rainfall")],
  center = train_center,
  scale = train_scale
) %>%
  as.data.frame()

train_nn <- RANN::nn2(
  data = training_scaled,
  query = training_scaled,
  k = 2
)

training_nn_distance <- train_nn$nn.dists[, 2]

mv_threshold_95 <- stats::quantile(
  training_nn_distance,
  mv_threshold_prob,
  na.rm = TRUE
)

mv_threshold_975 <- stats::quantile(
  training_nn_distance,
  mv_threshold_prob_strict,
  na.rm = TRUE
)

cat("\nMultivariate novelty threshold:\n")
cat("95% threshold:", mv_threshold_95, "\n")
cat("97.5% threshold:", mv_threshold_975, "\n")

projection_scaled <- scale(
  all_projection_climate[, c("Temperature", "Rainfall")],
  center = train_center,
  scale = train_scale
) %>%
  as.data.frame()

projection_nn <- RANN::nn2(
  data = training_scaled,
  query = projection_scaled,
  k = 1
)

all_projection_novelty <- all_projection_climate %>%
  dplyr::mutate(
    Climate_Distance = projection_nn$nn.dists[, 1],

    Temp_Outside_Strict = Temperature < temp_min | Temperature > temp_max,
    Rain_Outside_Strict = Rainfall < rain_min | Rainfall > rain_max,
    Outside_Strict_Range = Temp_Outside_Strict | Rain_Outside_Strict,

    Temp_Outside_Core = Temperature < temp_q025 | Temperature > temp_q975,
    Rain_Outside_Core = Rainfall < rain_q025 | Rainfall > rain_q975,
    Outside_Core_Range = Temp_Outside_Core | Rain_Outside_Core,

    Multivariate_Novel_95 = Climate_Distance > mv_threshold_95,
    Multivariate_Novel_975 = Climate_Distance > mv_threshold_975,

    Extrapolation_Class = dplyr::case_when(
      Outside_Strict_Range & Multivariate_Novel_95 ~ "Outside strict range + multivariate novel",
      Outside_Strict_Range ~ "Outside strict range only",
      Multivariate_Novel_95 ~ "Multivariate novel within univariate range",
      Outside_Core_Range ~ "Outside core range only",
      TRUE ~ "Within observed climate space"
    )
  )

write.csv(
  all_projection_novelty,
  file.path(novelty_out_dir, "Climate_novelty_all_projection_points_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


novelty_summary <- all_projection_novelty %>%
  dplyr::group_by(Projection_Type, SSP, Time_Period) %>%
  dplyr::summarise(
    n_points = dplyr::n(),
    mean_temperature = mean(Temperature, na.rm = TRUE),
    sd_temperature = stats::sd(Temperature, na.rm = TRUE),
    mean_rainfall = mean(Rainfall, na.rm = TRUE),
    sd_rainfall = stats::sd(Rainfall, na.rm = TRUE),

    strict_range_fraction = mean(Outside_Strict_Range, na.rm = TRUE),
    temp_strict_fraction = mean(Temp_Outside_Strict, na.rm = TRUE),
    rain_strict_fraction = mean(Rain_Outside_Strict, na.rm = TRUE),

    core_range_fraction = mean(Outside_Core_Range, na.rm = TRUE),
    multivariate_novel_95_fraction = mean(Multivariate_Novel_95, na.rm = TRUE),
    multivariate_novel_975_fraction = mean(Multivariate_Novel_975, na.rm = TRUE),

    mean_climate_distance = mean(Climate_Distance, na.rm = TRUE),
    q95_climate_distance = stats::quantile(Climate_Distance, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nNovelty summary by SSP and period:\n")
print(novelty_summary)

write.csv(
  novelty_summary,
  file.path(novelty_out_dir, "Climate_novelty_summary_by_SSP_period_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

novelty_latitude_summary <- all_projection_novelty %>%
  dplyr::mutate(
    abs_lat = abs(Latitude),
    latitude_band = dplyr::case_when(
      abs_lat < 23.5 ~ "Tropical",
      abs_lat >= 23.5 & abs_lat < 45 ~ "Subtropical",
      abs_lat >= 45 & abs_lat < 66.5 ~ "Mid-high latitude",
      abs_lat >= 66.5 ~ "Polar",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::group_by(Projection_Type, SSP, Time_Period, latitude_band) %>%
  dplyr::summarise(
    n_points = dplyr::n(),
    strict_range_fraction = mean(Outside_Strict_Range, na.rm = TRUE),
    temp_strict_fraction = mean(Temp_Outside_Strict, na.rm = TRUE),
    rain_strict_fraction = mean(Rain_Outside_Strict, na.rm = TRUE),
    core_range_fraction = mean(Outside_Core_Range, na.rm = TRUE),
    multivariate_novel_95_fraction = mean(Multivariate_Novel_95, na.rm = TRUE),
    multivariate_novel_975_fraction = mean(Multivariate_Novel_975, na.rm = TRUE),
    mean_climate_distance = mean(Climate_Distance, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nNovelty summary by latitude band:\n")
print(novelty_latitude_summary)

write.csv(
  novelty_latitude_summary,
  file.path(novelty_out_dir, "Climate_novelty_summary_by_latitude_band_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


novelty_summary_no_antarctica <- all_projection_novelty %>%
  dplyr::filter(Latitude > -60) %>%
  dplyr::group_by(Projection_Type, SSP, Time_Period) %>%
  dplyr::summarise(
    n_points = dplyr::n(),
    strict_range_fraction = mean(Outside_Strict_Range, na.rm = TRUE),
    core_range_fraction = mean(Outside_Core_Range, na.rm = TRUE),
    multivariate_novel_95_fraction = mean(Multivariate_Novel_95, na.rm = TRUE),
    multivariate_novel_975_fraction = mean(Multivariate_Novel_975, na.rm = TRUE),
    mean_climate_distance = mean(Climate_Distance, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  novelty_summary_no_antarctica,
  file.path(novelty_out_dir, "Climate_novelty_summary_excluding_Antarctica_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


novelty_future <- all_projection_novelty %>%
  dplyr::filter(Projection_Type == "Future ensemble") %>%
  dplyr::select(
    Longitude,
    Latitude,
    SSP,
    Time_Period,
    Climate_Distance,
    Outside_Strict_Range,
    Outside_Core_Range,
    Multivariate_Novel_95,
    Multivariate_Novel_975,
    Extrapolation_Class
  ) %>%
  dplyr::mutate(
    SSP = as.character(SSP),
    Time_Period = as.character(Time_Period)
  )

risk_change_novelty <- risk_change %>%
  dplyr::mutate(
    SSP = as.character(SSP),
    Time_Period = as.character(Time_Period)
  ) %>%
  dplyr::left_join(
    novelty_future,
    by = c("Longitude", "Latitude", "SSP", "Time_Period")
  )

risk_change_by_novelty <- risk_change_novelty %>%
  dplyr::group_by(SSP, Time_Period, Extrapolation_Class) %>%
  dplyr::summarise(
    n_points = dplyr::n(),
    mean_relative_change = mean(Relative_Change_Pct, na.rm = TRUE),
    median_relative_change = stats::median(Relative_Change_Pct, na.rm = TRUE),
    positive_fraction = mean(Relative_Change_Pct > 0, na.rm = TRUE),
    q05 = stats::quantile(Relative_Change_Pct, 0.05, na.rm = TRUE),
    q95 = stats::quantile(Relative_Change_Pct, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  risk_change_novelty,
  file.path(novelty_out_dir, "Risk_change_with_climate_novelty_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  risk_change_by_novelty,
  file.path(novelty_out_dir, "Risk_change_summary_by_climate_novelty_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nRisk-change summary by climate novelty class:\n")
print(risk_change_by_novelty)

risk_change_novelty_weighted <- risk_change_novelty %>%
  dplyr::mutate(
    lat_weight = pmax(cos(Latitude * pi / 180), 0),
    positive_change = Relative_Change_Pct > 0
  )

weighted_summary <- risk_change_novelty_weighted %>%
  dplyr::group_by(SSP, Time_Period) %>%
  dplyr::summarise(
    weighted_strict_novel_fraction =
      stats::weighted.mean(Outside_Strict_Range, lat_weight, na.rm = TRUE),
    weighted_core_novel_fraction =
      stats::weighted.mean(Outside_Core_Range, lat_weight, na.rm = TRUE),
    weighted_multivariate_novel_95_fraction =
      stats::weighted.mean(Multivariate_Novel_95, lat_weight, na.rm = TRUE),
    weighted_positive_change_fraction =
      stats::weighted.mean(positive_change, lat_weight, na.rm = TRUE),
    weighted_mean_relative_change =
      stats::weighted.mean(Relative_Change_Pct, lat_weight, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  weighted_summary,
  file.path(novelty_out_dir, "Weighted_climate_novelty_and_risk_change_summary_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nWeighted novelty and risk-change summary:\n")
print(weighted_summary)


novelty_palette <- c(
  "Within observed climate space" = "#4DAF4A",
  "Outside core range only" = "#FEE08B",
  "Multivariate novel within univariate range" = "#FDAE61",
  "Outside strict range only" = "#D73027",
  "Outside strict range + multivariate novel" = "#7F0000"
)

# Ensure class order
all_projection_novelty <- all_projection_novelty %>%
  dplyr::mutate(
    Extrapolation_Class = factor(
      Extrapolation_Class,
      levels = names(novelty_palette)
    )
  )

current_novelty_sf <- all_projection_novelty %>%
  dplyr::filter(Projection_Type == "Current baseline") %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

p_current_novelty <- ggplot2::ggplot(current_novelty_sf) +
  ggplot2::geom_sf(ggplot2::aes(color = Extrapolation_Class), size = 0.45, alpha = 0.9) +
  ggplot2::scale_color_manual(
    values = novelty_palette,
    name = "Climate novelty",
    drop = FALSE
  ) +
  ggplot2::coord_sf(expand = FALSE) +
  ggplot2::labs(title = "Climate novelty of contemporary baseline prediction points") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    axis.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(novelty_out_dir, "Climate_novelty_current_baseline_corrected.pdf"),
  p_current_novelty,
  width = 12,
  height = 6,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)

future_novelty_sf <- all_projection_novelty %>%
  dplyr::filter(Projection_Type == "Future ensemble") %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

p_future_novelty <- ggplot2::ggplot(future_novelty_sf) +
  ggplot2::geom_sf(ggplot2::aes(color = Extrapolation_Class), size = 0.35, alpha = 0.9) +
  ggplot2::scale_color_manual(
    values = novelty_palette,
    name = "Climate novelty",
    drop = FALSE
  ) +
  ggplot2::facet_grid(SSP ~ Time_Period) +
  ggplot2::coord_sf(expand = FALSE) +
  ggplot2::labs(title = "Climate novelty under future climate scenarios") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    axis.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(novelty_out_dir, "Climate_novelty_future_GCMensemble_4SSP_2periods_corrected.pdf"),
  p_future_novelty,
  width = 16,
  height = 14,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)

future_distance_sf <- all_projection_novelty %>%
  dplyr::filter(Projection_Type == "Future ensemble") %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

p_future_distance <- ggplot2::ggplot(future_distance_sf) +
  ggplot2::geom_sf(ggplot2::aes(color = Climate_Distance), size = 0.35, alpha = 0.9) +
  ggplot2::scale_color_gradientn(
    colors = c("#F7FCF0", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8"),
    name = "Climate\ndistance",
    oob = scales::oob_squish
  ) +
  ggplot2::facet_grid(SSP ~ Time_Period) +
  ggplot2::coord_sf(expand = FALSE) +
  ggplot2::labs(title = "Distance to observed training climate space") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    axis.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(novelty_out_dir, "Climate_distance_future_GCMensemble_4SSP_2periods_corrected.pdf"),
  p_future_distance,
  width = 16,
  height = 14,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)

