
packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "caret",
  "randomForest", "RColorBrewer", "scales", "tibble",
  "viridisLite"
)

need_install <- packages[!packages %in% rownames(installed.packages())]
if (length(need_install) > 0) {
  install.packages(need_install)
}

invisible(lapply(packages, library, character.only = TRUE))

set.seed(123)



ml_file <- "ml.xlsx"

rf_model_path   <- "Best_RF_Model.rds"
dummies_path    <- "dummies_rf.rds"
nzv_path        <- "nzv_rf.rds"
preProc_path    <- "preProc_rf.rds"
train_mode_path <- "train_mode_rf.rds"

out_dir <- "D:/R/R/Tgondii_response_surface_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_temp <- 120
n_rain <- 120

core_lower <- 0.025
core_upper <- 0.975

surface_year <- 2020

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
    !is.na(Temperature),
    !is.na(Rainfall),
    Prevalence >= 0,
    Prevalence <= 1
  )

if (median(model_data$Rainfall, na.rm = TRUE) > 50) {
  stop("Rainfall in ml.xlsx appears to be annual precipitation, not mm/day. Please check units.")
}

feature_vars <- colnames(model_data)[colnames(model_data) != "Prevalence"]

cat("\nTraining climate summary:\n")
print(summary(model_data[, c("Temperature", "Rainfall")]))

climate_domain <- model_data %>%
  dplyr::summarise(
    temp_min = min(Temperature, na.rm = TRUE),
    temp_q025 = quantile(Temperature, core_lower, na.rm = TRUE),
    temp_q975 = quantile(Temperature, core_upper, na.rm = TRUE),
    temp_max = max(Temperature, na.rm = TRUE),
    rain_min = min(Rainfall, na.rm = TRUE),
    rain_q025 = quantile(Rainfall, core_lower, na.rm = TRUE),
    rain_q975 = quantile(Rainfall, core_upper, na.rm = TRUE),
    rain_max = max(Rainfall, na.rm = TRUE)
  )

write.csv(
  climate_domain,
  file.path(out_dir, "Response_surface_training_climate_domain.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

print(climate_domain)


get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

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

get_required_predictors <- function(model) {
  if (inherits(model, "train")) {
    if (!is.null(model$trainingData)) {
      return(setdiff(colnames(model$trainingData), ".outcome"))
    }
    if (!is.null(model$finalModel$importance)) {
      return(rownames(model$finalModel$importance))
    }
  }

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

  if (length(cols) == 0) return(NULL)
  cols
}

align_columns <- function(x, required_cols, fill_value = 0) {
  x <- as.data.frame(x)

  missing_cols <- setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0) {
    for (nm in missing_cols) x[[nm]] <- fill_value
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

make_prediction_data <- function(clim_df,
                                 year_value = 2020,
                                 diet_value = NULL,
                                 habitat_value = NULL,
                                 label_value = "Overall") {

  pred_df <- clim_df %>%
    dplyr::mutate(
      Year = year_value,
      Total = median(model_data$Total, na.rm = TRUE),
      Longitude = median(model_data$Longitude, na.rm = TRUE),
      Latitude = median(model_data$Latitude, na.rm = TRUE),
      Order = get_train_mode_value("Order"),
      Family = get_train_mode_value("Family"),
      Diet = ifelse(is.null(diet_value), get_train_mode_value("Diet"), diet_value),
      Habitat = ifelse(is.null(habitat_value), get_train_mode_value("Habitat"), habitat_value),
      Sample_Type = get_train_mode_value("Sample_Type"),
      Methods = get_train_mode_value("Methods"),
      Continent = get_train_mode_value("Continent"),
      Country = get_train_mode_value("Country"),
      Scenario_Label = label_value
    ) %>%
    dplyr::select(dplyr::all_of(feature_vars), Scenario_Label)

  pred_features <- pred_df %>%
    dplyr::select(dplyr::all_of(feature_vars))

  if (!is.null(dummies$lvls)) {
    for (v in names(dummies$lvls)) {
      if (v %in% colnames(pred_features)) {
        pred_features[[v]] <- factor(
          as.character(pred_features[[v]]),
          levels = dummies$lvls[[v]]
        )
      }
    }
  } else {
    pred_features <- pred_features %>%
      dplyr::mutate(dplyr::across(where(is.character), as.factor))
  }

  list(
    pred_features = pred_features,
    scenario_label = pred_df$Scenario_Label
  )
}

predict_rf_prevalence <- function(pred_df) {
  pred_df$Prevalence <- 0

  dummy_mat <- predict(dummies, newdata = pred_df)
  dummy_df <- as.data.frame(dummy_mat)

  dummy_df <- dummy_df[, !colnames(dummy_df) %in% "Prevalence", drop = FALSE]

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

  if (!is.null(preproc_predictors)) {
    dummy_df <- align_columns(dummy_df, preproc_predictors, fill_value = 0)
  }

  x <- as.data.frame(predict(preProc, dummy_df))
  x <- align_columns(x, required_predictors, fill_value = 0)

  pred <- as.numeric(predict(rf_model, newdata = x))
  pred <- pmin(pmax(pred, 0), 1)

  pred
}

predict_surface <- function(clim_grid,
                            diet_value = NULL,
                            habitat_value = NULL,
                            label_value = "Overall") {

  tmp <- make_prediction_data(
    clim_df = clim_grid,
    year_value = surface_year,
    diet_value = diet_value,
    habitat_value = habitat_value,
    label_value = label_value
  )

  pred <- predict_rf_prevalence(tmp$pred_features)

  clim_grid %>%
    dplyr::mutate(
      Predicted_Risk = pred,
      Scenario_Label = label_value,
      Diet_Setting = ifelse(is.null(diet_value), get_train_mode_value("Diet"), diet_value),
      Habitat_Setting = ifelse(is.null(habitat_value), get_train_mode_value("Habitat"), habitat_value)
    )
}


temp_seq <- seq(
  climate_domain$temp_q025[1],
  climate_domain$temp_q975[1],
  length.out = n_temp
)

rain_seq <- seq(
  climate_domain$rain_q025[1],
  climate_domain$rain_q975[1],
  length.out = n_rain
)

clim_grid <- expand.grid(
  Temperature = temp_seq,
  Rainfall = rain_seq
) %>%
  as_tibble()


surface_overall <- predict_surface(
  clim_grid = clim_grid,
  label_value = "Overall"
)

write.csv(
  surface_overall,
  file.path(out_dir, "Response_surface_overall_predictions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

surface_summary <- surface_overall %>%
  dplyr::summarise(
    min_predicted_risk = min(Predicted_Risk, na.rm = TRUE),
    q25_predicted_risk = quantile(Predicted_Risk, 0.25, na.rm = TRUE),
    median_predicted_risk = median(Predicted_Risk, na.rm = TRUE),
    mean_predicted_risk = mean(Predicted_Risk, na.rm = TRUE),
    q75_predicted_risk = quantile(Predicted_Risk, 0.75, na.rm = TRUE),
    max_predicted_risk = max(Predicted_Risk, na.rm = TRUE)
  )

max_cell <- surface_overall %>%
  dplyr::slice_max(Predicted_Risk, n = 1, with_ties = FALSE) %>%
  dplyr::select(Temperature, Rainfall, Predicted_Risk)

high_threshold <- quantile(surface_overall$Predicted_Risk, 0.90, na.rm = TRUE)

high_risk_envelope <- surface_overall %>%
  dplyr::filter(Predicted_Risk >= high_threshold) %>%
  dplyr::summarise(
    high_risk_threshold = high_threshold,
    temp_min_high = min(Temperature, na.rm = TRUE),
    temp_median_high = median(Temperature, na.rm = TRUE),
    temp_max_high = max(Temperature, na.rm = TRUE),
    rain_min_high = min(Rainfall, na.rm = TRUE),
    rain_median_high = median(Rainfall, na.rm = TRUE),
    rain_max_high = max(Rainfall, na.rm = TRUE)
  )

summary_out <- dplyr::bind_cols(surface_summary, max_cell, high_risk_envelope)

write.csv(
  summary_out,
  file.path(out_dir, "Response_surface_overall_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nOverall response surface summary:\n")
print(summary_out)

training_climate_plot <- model_data %>%
  dplyr::select(Temperature, Rainfall, Prevalence)

p_surface <- ggplot2::ggplot(surface_overall, ggplot2::aes(x = Temperature, y = Rainfall)) +
  ggplot2::geom_raster(ggplot2::aes(fill = Predicted_Risk), interpolate = TRUE) +
  ggplot2::geom_contour(
    ggplot2::aes(z = Predicted_Risk),
    colour = "grey30",
    alpha = 0.55,
    linewidth = 0.25
  ) +
  ggplot2::geom_point(
    data = training_climate_plot,
    ggplot2::aes(x = Temperature, y = Rainfall),
    inherit.aes = FALSE,
    colour = "black",
    alpha = 0.18,
    size = 0.7
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c("#FFF7BC", "#FEC44F", "#FE9929", "#D95F0E", "#993404"),
    name = "Predicted\nrisk",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  ggplot2::labs(
    x = "Annual mean temperature (°C)",
    y = "Mean daily precipitation (mm/day)",
    title = "Temperature–precipitation response surface for predicted T. gondii infection risk",
    subtitle = "Non-climatic covariates held at representative training values; points show observed training climates"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(out_dir, "Response_surface_temperature_precipitation_overall.pdf"),
  p_surface,
  width = 9,
  height = 7,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)


available_diet <- sort(unique(na.omit(as.character(model_data$Diet))))
available_habitat <- sort(unique(na.omit(as.character(model_data$Habitat))))

cat("\nAvailable Diet levels:\n")
print(available_diet)

cat("\nAvailable Habitat levels:\n")
print(available_habitat)

diet_candidates <- available_diet[available_diet %in% c(
  "Carnivore", "Carnivorous", "Herbivore", "Herbivorous", "Omnivore", "Omnivorous"
)]

if (length(diet_candidates) == 0) {
  diet_candidates <- names(sort(table(model_data$Diet), decreasing = TRUE))[1:min(3, length(unique(model_data$Diet)))]
}

habitat_candidates <- available_habitat[available_habitat %in% c(
  "Freshwater", "Terrestrial", "Marine", "Semi-aquatic", "Aquatic"
)]

if (length(habitat_candidates) == 0) {
  habitat_candidates <- names(sort(table(model_data$Habitat), decreasing = TRUE))[1:min(3, length(unique(model_data$Habitat)))]
}

# 7.1 Surfaces by habitat, holding diet at modal value
surface_by_habitat <- dplyr::bind_rows(lapply(habitat_candidates, function(hab) {
  predict_surface(
    clim_grid = clim_grid,
    habitat_value = hab,
    label_value = paste0("Habitat: ", hab)
  )
}))

write.csv(
  surface_by_habitat,
  file.path(out_dir, "Response_surface_by_habitat_predictions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (nrow(surface_by_habitat) > 0) {
  p_habitat <- ggplot2::ggplot(surface_by_habitat, ggplot2::aes(x = Temperature, y = Rainfall)) +
    ggplot2::geom_raster(ggplot2::aes(fill = Predicted_Risk), interpolate = TRUE) +
    ggplot2::geom_contour(
      ggplot2::aes(z = Predicted_Risk),
      colour = "grey30",
      alpha = 0.45,
      linewidth = 0.2
    ) +
    ggplot2::facet_wrap(~ Habitat_Setting) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFF7BC", "#FEC44F", "#FE9929", "#D95F0E", "#993404"),
      name = "Predicted\nrisk",
      labels = scales::number_format(accuracy = 0.01)
    ) +
    ggplot2::labs(
      x = "Annual mean temperature (°C)",
      y = "Mean daily precipitation (mm/day)",
      title = "Temperature–precipitation response surfaces by habitat setting"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(out_dir, "Response_surface_temperature_precipitation_by_habitat.pdf"),
    p_habitat,
    width = 12,
    height = 7,
    dpi = 600,
    device = "pdf",
    useDingbats = FALSE
  )
}

# 7.2 Surfaces by diet, holding habitat at modal value
surface_by_diet <- dplyr::bind_rows(lapply(diet_candidates, function(diet) {
  predict_surface(
    clim_grid = clim_grid,
    diet_value = diet,
    label_value = paste0("Diet: ", diet)
  )
}))

write.csv(
  surface_by_diet,
  file.path(out_dir, "Response_surface_by_diet_predictions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (nrow(surface_by_diet) > 0) {
  p_diet <- ggplot2::ggplot(surface_by_diet, ggplot2::aes(x = Temperature, y = Rainfall)) +
    ggplot2::geom_raster(ggplot2::aes(fill = Predicted_Risk), interpolate = TRUE) +
    ggplot2::geom_contour(
      ggplot2::aes(z = Predicted_Risk),
      colour = "grey30",
      alpha = 0.45,
      linewidth = 0.2
    ) +
    ggplot2::facet_wrap(~ Diet_Setting) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFF7BC", "#FEC44F", "#FE9929", "#D95F0E", "#993404"),
      name = "Predicted\nrisk",
      labels = scales::number_format(accuracy = 0.01)
    ) +
    ggplot2::labs(
      x = "Annual mean temperature (°C)",
      y = "Mean daily precipitation (mm/day)",
      title = "Temperature–precipitation response surfaces by dietary guild setting"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      panel.grid = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(out_dir, "Response_surface_temperature_precipitation_by_diet.pdf"),
    p_diet,
    width = 12,
    height = 7,
    dpi = 600,
    device = "pdf",
    useDingbats = FALSE
  )
}


surface_zones <- surface_overall %>%
  dplyr::mutate(
    temp_zone = dplyr::case_when(
      Temperature < 5 ~ "Cold (<5°C)",
      Temperature >= 5 & Temperature < 15 ~ "Cool-temperate (5-15°C)",
      Temperature >= 15 & Temperature < 25 ~ "Warm (15-25°C)",
      Temperature >= 25 ~ "Hot (≥25°C)"
    ),
    rain_zone = dplyr::case_when(
      Rainfall < 1 ~ "Dry (<1 mm/day)",
      Rainfall >= 1 & Rainfall < 3 ~ "Moderate (1-3 mm/day)",
      Rainfall >= 3 ~ "Wet (≥3 mm/day)"
    )
  ) %>%
  dplyr::group_by(temp_zone, rain_zone) %>%
  dplyr::summarise(
    mean_predicted_risk = mean(Predicted_Risk, na.rm = TRUE),
    median_predicted_risk = median(Predicted_Risk, na.rm = TRUE),
    q25 = quantile(Predicted_Risk, 0.25, na.rm = TRUE),
    q75 = quantile(Predicted_Risk, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  surface_zones,
  file.path(out_dir, "Response_surface_climate_zone_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nClimate-zone response surface summary:\n")
print(surface_zones)
