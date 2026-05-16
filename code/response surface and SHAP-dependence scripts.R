
packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "caret",
  "randomForest", "fastshap", "tibble", "scales",
  "RColorBrewer", "forcats"
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

out_dir <- "D:/R/R/Tgondii_SHAP_interaction_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

shap_nsim <- 50

max_shap_n <- 500

n_grid <- 80


rf_model   <- readRDS(rf_model_path)
dummies    <- readRDS(dummies_path)
nzv        <- readRDS(nzv_path)
preProc    <- readRDS(preProc_path)
train_mode <- readRDS(train_mode_path)

data_raw <- readxl::read_excel(ml_file, sheet = "Sheet1")

if ("Sample Type" %in% names(data_raw) && !("Sample_Type" %in% names(data_raw))) {
  data_raw <- data_raw %>% dplyr::rename(Sample_Type = `Sample Type`)
}
if ("Rainfall (mm/day)" %in% names(data_raw) && !("Rainfall" %in% names(data_raw))) {
  data_raw <- data_raw %>% dplyr::rename(Rainfall = `Rainfall (mm/day)`)
}

required_raw_cols <- c(
  "Prevalence", "Year", "Rainfall", "Temperature", "Longitude", "Latitude",
  "Total", "Order", "Family", "Diet", "Habitat", "Sample_Type",
  "Methods", "Continent", "Country"
)

missing_raw_cols <- setdiff(required_raw_cols, names(data_raw))
if (length(missing_raw_cols) > 0) {
  stop("ml.xlsx is missing required columns: ", paste(missing_raw_cols, collapse = ", "))
}

model_data <- data_raw %>%
  dplyr::select(dplyr::all_of(required_raw_cols)) %>%
  dplyr::filter(
    !is.na(Prevalence),
    !is.na(Rainfall),
    !is.na(Temperature),
    !is.na(Latitude),
    Prevalence >= 0,
    Prevalence <= 1
  )

if (median(model_data$Rainfall, na.rm = TRUE) > 50) {
  stop("Rainfall appears to be annual precipitation. Expected mm/day.")
}

feature_vars <- colnames(model_data)[colnames(model_data) != "Prevalence"]

cat("\nModel data dimensions:\n")
print(dim(model_data))
cat("\nRainfall summary, expected mm/day:\n")
print(summary(model_data$Rainfall))


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

prepare_rf_newdata <- function(newdata_features) {
  newdata_features <- as.data.frame(newdata_features)

  if (!is.null(dummies$lvls)) {
    for (v in names(dummies$lvls)) {
      if (v %in% colnames(newdata_features)) {
        newdata_features[[v]] <- factor(
          as.character(newdata_features[[v]]),
          levels = dummies$lvls[[v]]
        )
      }
    }
  } else {
    newdata_features <- newdata_features %>%
      dplyr::mutate(dplyr::across(where(is.character), as.factor))
  }

  newdata_features$Prevalence <- 0

  dummy_mat <- predict(dummies, newdata = newdata_features)
  dummy_df <- as.data.frame(dummy_mat)
  dummy_df <- dummy_df[, !colnames(dummy_df) %in% "Prevalence", drop = FALSE]

  if (length(nzv) > 0) {
    if (is.numeric(nzv)) {
      idx <- nzv[nzv <= ncol(dummy_df)]
      if (length(idx) > 0) dummy_df <- dummy_df[, -idx, drop = FALSE]
    } else {
      dummy_df <- dummy_df[, setdiff(colnames(dummy_df), nzv), drop = FALSE]
    }
  }

  if (!is.null(preproc_predictors)) {
    dummy_df <- align_columns(dummy_df, preproc_predictors, fill_value = 0)
  }

  x <- as.data.frame(predict(preProc, dummy_df))
  x <- align_columns(x, required_predictors, fill_value = 0)

  x
}

predict_rf_prevalence <- function(newdata_features) {
  x <- prepare_rf_newdata(newdata_features)
  pred <- as.numeric(predict(rf_model, newdata = x))
  pred <- pmin(pmax(pred, 0), 1)
  pred
}

pred_wrapper <- function(object, newdata) {
  predict_rf_prevalence(newdata)
}


X <- model_data %>%
  dplyr::select(dplyr::all_of(feature_vars))

complete_idx <- stats::complete.cases(X)
X_complete <- X[complete_idx, , drop = FALSE]
metadata_complete <- model_data[complete_idx, , drop = FALSE]

if (nrow(X_complete) > max_shap_n) {
  set.seed(123)
  shap_rows <- sample(seq_len(nrow(X_complete)), max_shap_n)
  X_shap <- X_complete[shap_rows, , drop = FALSE]
  metadata_shap <- metadata_complete[shap_rows, , drop = FALSE]
} else {
  X_shap <- X_complete
  metadata_shap <- metadata_complete
}

cat("\nRows used for SHAP analysis:\n")
print(nrow(X_shap))


shap_features <- c("Rainfall", "Temperature")

shap_values <- fastshap::explain(
  object = rf_model,
  X = X_shap,
  pred_wrapper = pred_wrapper,
  nsim = shap_nsim,
  feature_names = shap_features,
  adjust = TRUE
)

shap_df <- as.data.frame(shap_values)

if (!all(shap_features %in% colnames(shap_df))) {
  stop(
    "SHAP output does not contain expected columns: ",
    paste(shap_features, collapse = ", "),
    ". Actual columns are: ",
    paste(colnames(shap_df), collapse = ", ")
  )
}

shap_wide <- metadata_shap %>%
  dplyr::select(
    Rainfall, Temperature, Latitude, Habitat, Diet, Prevalence,
    Continent, Country
  ) %>%
  dplyr::mutate(
    SHAP_Rainfall = shap_df$Rainfall,
    SHAP_Temperature = shap_df$Temperature
  )

shap_long <- shap_wide %>%
  tidyr::pivot_longer(
    cols = c(SHAP_Rainfall, SHAP_Temperature),
    names_to = "Feature",
    values_to = "SHAP_value"
  ) %>%
  dplyr::mutate(
    Feature = dplyr::recode(
      Feature,
      "SHAP_Rainfall" = "Rainfall",
      "SHAP_Temperature" = "Temperature"
    )
  )

write.csv(
  shap_long,
  file.path(out_dir, "SHAP_values_Rainfall_Temperature_long.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  shap_wide,
  file.path(out_dir, "SHAP_values_Rainfall_Temperature_wide.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)



available_habitat <- sort(unique(na.omit(as.character(shap_wide$Habitat))))
cat("\nAvailable Habitat levels:\n")
print(available_habitat)

freshwater_levels <- available_habitat[
  grepl("fresh|aquatic|water|river|wetland", available_habitat, ignore.case = TRUE)
]

if (length(freshwater_levels) == 0) {
  message("No freshwater-like Habitat level detected. Please inspect available Habitat levels.")
}

shap_wide <- shap_wide %>%
  dplyr::mutate(
    Freshwater_Associated = ifelse(
      as.character(Habitat) %in% freshwater_levels,
      "Freshwater-associated",
      "Other habitats"
    ),
    Freshwater_Associated = factor(
      Freshwater_Associated,
      levels = c("Other habitats", "Freshwater-associated")
    )
  )

rainfall_habitat_summary <- shap_wide %>%
  dplyr::group_by(Freshwater_Associated) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_rainfall = mean(Rainfall, na.rm = TRUE),
    mean_SHAP_Rainfall = mean(SHAP_Rainfall, na.rm = TRUE),
    median_SHAP_Rainfall = median(SHAP_Rainfall, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  rainfall_habitat_summary,
  file.path(out_dir, "Interaction_summary_SHAP_Rainfall_by_FreshwaterHabitat.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_rain_hab <- ggplot2::ggplot(
  shap_wide,
  ggplot2::aes(x = Rainfall, y = SHAP_Rainfall, colour = Freshwater_Associated)
) +
  ggplot2::geom_point(alpha = 0.75, size = 1.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  ggplot2::labs(
    x = "Mean daily precipitation (mm/day)",
    y = "SHAP value for precipitation",
    colour = NULL,
    title = "SHAP dependence: precipitation effect by freshwater-associated habitat"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "SHAP_dependence_Rainfall_by_FreshwaterHabitat.pdf"),
  p_rain_hab,
  width = 8,
  height = 6,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)


shap_wide <- shap_wide %>%
  dplyr::mutate(
    abs_lat = abs(Latitude),
    latitude_band = dplyr::case_when(
      abs_lat < 23.5 ~ "Tropical",
      abs_lat >= 23.5 & abs_lat < 45 ~ "Subtropical",
      abs_lat >= 45 & abs_lat < 66.5 ~ "Mid-high latitude",
      abs_lat >= 66.5 ~ "Polar",
      TRUE ~ NA_character_
    ),
    latitude_band = factor(
      latitude_band,
      levels = c("Tropical", "Subtropical", "Mid-high latitude", "Polar")
    )
  )

temperature_latitude_summary <- shap_wide %>%
  dplyr::group_by(latitude_band) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_temperature = mean(Temperature, na.rm = TRUE),
    mean_abs_latitude = mean(abs_lat, na.rm = TRUE),
    mean_SHAP_Temperature = mean(SHAP_Temperature, na.rm = TRUE),
    median_SHAP_Temperature = median(SHAP_Temperature, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  temperature_latitude_summary,
  file.path(out_dir, "Interaction_summary_SHAP_Temperature_by_LatitudeBand.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_temp_lat <- ggplot2::ggplot(
  shap_wide,
  ggplot2::aes(x = Temperature, y = SHAP_Temperature, colour = latitude_band)
) +
  ggplot2::geom_point(alpha = 0.75, size = 1.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  ggplot2::labs(
    x = "Annual mean temperature (°C)",
    y = "SHAP value for temperature",
    colour = "Latitude band",
    title = "SHAP dependence: temperature effect by latitude band"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "SHAP_dependence_Temperature_by_LatitudeBand.pdf"),
  p_temp_lat,
  width = 8,
  height = 6,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)


available_diet <- sort(unique(na.omit(as.character(shap_wide$Diet))))
cat("\nAvailable Diet levels:\n")
print(available_diet)

# Collapse rare diet categories if necessary
diet_counts <- shap_wide %>%
  dplyr::count(Diet, sort = TRUE)

major_diets <- diet_counts %>%
  dplyr::filter(n >= 10) %>%
  dplyr::pull(Diet) %>%
  as.character()

shap_wide <- shap_wide %>%
  dplyr::mutate(
    Diet_Grouped = ifelse(as.character(Diet) %in% major_diets, as.character(Diet), "Other"),
    Diet_Grouped = factor(Diet_Grouped)
  )

diet_rainfall_summary <- shap_wide %>%
  dplyr::group_by(Diet_Grouped) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_rainfall = mean(Rainfall, na.rm = TRUE),
    mean_SHAP_Rainfall = mean(SHAP_Rainfall, na.rm = TRUE),
    median_SHAP_Rainfall = median(SHAP_Rainfall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(mean_SHAP_Rainfall))

write.csv(
  diet_rainfall_summary,
  file.path(out_dir, "Interaction_summary_SHAP_Rainfall_by_Diet.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_diet_rain <- ggplot2::ggplot(
  shap_wide,
  ggplot2::aes(x = Rainfall, y = SHAP_Rainfall, colour = Diet_Grouped)
) +
  ggplot2::geom_point(alpha = 0.7, size = 1.5) +
  ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  ggplot2::labs(
    x = "Mean daily precipitation (mm/day)",
    y = "SHAP value for precipitation",
    colour = "Dietary guild",
    title = "SHAP dependence: precipitation effect by dietary guild"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "SHAP_dependence_Rainfall_by_Diet.pdf"),
  p_diet_rain,
  width = 8,
  height = 6,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)


make_reference_data <- function(n_rows) {
  ref <- tibble::tibble(
    Year = rep(2020, n_rows),
    Rainfall = rep(median(model_data$Rainfall, na.rm = TRUE), n_rows),
    Temperature = rep(median(model_data$Temperature, na.rm = TRUE), n_rows),
    Longitude = rep(median(model_data$Longitude, na.rm = TRUE), n_rows),
    Latitude = rep(median(model_data$Latitude, na.rm = TRUE), n_rows),
    Total = rep(median(model_data$Total, na.rm = TRUE), n_rows),
    Order = rep(names(sort(table(model_data$Order), decreasing = TRUE))[1], n_rows),
    Family = rep(names(sort(table(model_data$Family), decreasing = TRUE))[1], n_rows),
    Diet = rep(names(sort(table(model_data$Diet), decreasing = TRUE))[1], n_rows),
    Habitat = rep(names(sort(table(model_data$Habitat), decreasing = TRUE))[1], n_rows),
    Sample_Type = rep(names(sort(table(model_data$Sample_Type), decreasing = TRUE))[1], n_rows),
    Methods = rep(names(sort(table(model_data$Methods), decreasing = TRUE))[1], n_rows),
    Continent = rep(names(sort(table(model_data$Continent), decreasing = TRUE))[1], n_rows),
    Country = rep(names(sort(table(model_data$Country), decreasing = TRUE))[1], n_rows)
  )

  ref %>% dplyr::select(dplyr::all_of(feature_vars))
}

rain_grid <- seq(
  quantile(model_data$Rainfall, 0.025, na.rm = TRUE),
  quantile(model_data$Rainfall, 0.975, na.rm = TRUE),
  length.out = n_grid
)

temp_grid <- seq(
  quantile(model_data$Temperature, 0.025, na.rm = TRUE),
  quantile(model_data$Temperature, 0.975, na.rm = TRUE),
  length.out = n_grid
)

habitat_levels_for_pd <- c("Other habitats", "Freshwater-associated")

pd_rain_hab <- dplyr::bind_rows(lapply(habitat_levels_for_pd, function(hab_group) {
  ref <- make_reference_data(length(rain_grid))
  ref$Rainfall <- rain_grid

  if (hab_group == "Freshwater-associated" && length(freshwater_levels) > 0) {
    ref$Habitat <- freshwater_levels[1]
  } else {
    modal_hab <- names(sort(table(model_data$Habitat), decreasing = TRUE))[1]
    if (modal_hab %in% freshwater_levels) {
      other_levels <- setdiff(available_habitat, freshwater_levels)
      if (length(other_levels) > 0) modal_hab <- other_levels[1]
    }
    ref$Habitat <- modal_hab
  }

  tibble::tibble(
    Rainfall = rain_grid,
    Predicted_Risk = predict_rf_prevalence(ref),
    Group = hab_group
  )
}))

write.csv(
  pd_rain_hab,
  file.path(out_dir, "Partial_dependence_Rainfall_by_FreshwaterHabitat.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_pd_rain_hab <- ggplot2::ggplot(
  pd_rain_hab,
  ggplot2::aes(x = Rainfall, y = Predicted_Risk, colour = Group)
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::labs(
    x = "Mean daily precipitation (mm/day)",
    y = "Predicted infection risk",
    colour = NULL,
    title = "Partial dependence: precipitation by freshwater-associated habitat"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "Partial_dependence_Rainfall_by_FreshwaterHabitat.pdf"),
  p_pd_rain_hab,
  width = 8,
  height = 6,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)

lat_values <- tibble::tibble(
  latitude_band = c("Tropical", "Subtropical", "Mid-high latitude"),
  Latitude = c(0, 35, 55)
)

pd_temp_lat <- dplyr::bind_rows(lapply(seq_len(nrow(lat_values)), function(i) {
  ref <- make_reference_data(length(temp_grid))
  ref$Temperature <- temp_grid
  ref$Latitude <- lat_values$Latitude[i]

  tibble::tibble(
    Temperature = temp_grid,
    Predicted_Risk = predict_rf_prevalence(ref),
    latitude_band = lat_values$latitude_band[i]
  )
}))

write.csv(
  pd_temp_lat,
  file.path(out_dir, "Partial_dependence_Temperature_by_LatitudeBand.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_pd_temp_lat <- ggplot2::ggplot(
  pd_temp_lat,
  ggplot2::aes(x = Temperature, y = Predicted_Risk, colour = latitude_band)
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::labs(
    x = "Annual mean temperature (°C)",
    y = "Predicted infection risk",
    colour = "Latitude band",
    title = "Partial dependence: temperature by latitude setting"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "Partial_dependence_Temperature_by_LatitudeBand.pdf"),
  p_pd_temp_lat,
  width = 8,
  height = 6,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)

# 9.3 Rainfall curves by diet
diet_levels_for_pd <- major_diets
if (length(diet_levels_for_pd) == 0) {
  diet_levels_for_pd <- names(sort(table(model_data$Diet), decreasing = TRUE))[1]
}

pd_rain_diet <- dplyr::bind_rows(lapply(diet_levels_for_pd, function(diet_value) {
  ref <- make_reference_data(length(rain_grid))
  ref$Rainfall <- rain_grid
  ref$Diet <- diet_value

  tibble::tibble(
    Rainfall = rain_grid,
    Predicted_Risk = predict_rf_prevalence(ref),
    Diet_Grouped = diet_value
  )
}))

write.csv(
  pd_rain_diet,
  file.path(out_dir, "Partial_dependence_Rainfall_by_Diet.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

p_pd_rain_diet <- ggplot2::ggplot(
  pd_rain_diet,
  ggplot2::aes(x = Rainfall, y = Predicted_Risk, colour = Diet_Grouped)
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::labs(
    x = "Mean daily precipitation (mm/day)",
    y = "Predicted infection risk",
    colour = "Dietary guild",
    title = "Partial dependence: precipitation by dietary guild"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(out_dir, "Partial_dependence_Rainfall_by_Diet.pdf"),
  p_pd_rain_diet,
  width = 8,
  height = 6,
  dpi = 600,
  device = "pdf",
  useDingbats = FALSE
)

