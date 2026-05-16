# ==========================================================
# Spatial block cross-validation for RF model
# Based on ml.xlsx
# ==========================================================

library(readxl)
library(dplyr)
library(caret)
library(randomForest)
library(tibble)
library(purrr)

set.seed(123)

data_raw <- readxl::read_excel("ml.xlsx", sheet = "Sheet1")

model_data <- data_raw %>%
  dplyr::rename(
    Event = Event,
    Total = Total,
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

str(model_data)
summary(model_data$Prevalence)

block_size <- 10

model_data <- model_data %>%
  dplyr::mutate(
    lon_block = floor((Longitude + 180) / block_size),
    lat_block = floor((Latitude + 90) / block_size),
    spatial_block = paste(lon_block, lat_block, sep = "_")
  )

block_counts <- model_data %>%
  dplyr::count(spatial_block, name = "n") %>%
  dplyr::arrange(desc(n))

print(block_counts)

valid_blocks <- block_counts %>%
  dplyr::filter(n >= 3) %>%
  dplyr::pull(spatial_block)

ml_data <- model_data %>%
  dplyr::filter(spatial_block %in% valid_blocks)

cat("Number of records after removing sparse blocks:", nrow(ml_data), "\n")
cat("Number of spatial blocks:", length(unique(ml_data$spatial_block)), "\n")

k <- 5

unique_blocks <- unique(ml_data$spatial_block)

block_fold_table <- tibble(
  spatial_block = unique_blocks,
  fold = sample(rep(1:k, length.out = length(unique_blocks)))
)

ml_data <- ml_data %>%
  dplyr::left_join(block_fold_table, by = "spatial_block")

print(table(ml_data$fold))

get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

calc_metrics <- function(obs, pred) {
  obs <- as.numeric(obs)
  pred <- as.numeric(pred)
  
  # If obs has no variance in a fold, R2 cannot be calculated
  r2 <- ifelse(
    sd(obs, na.rm = TRUE) == 0,
    NA_real_,
    cor(obs, pred, use = "complete.obs")^2
  )
  
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae <- mean(abs(obs - pred), na.rm = TRUE)
  
  tibble(
    R2 = r2,
    RMSE = rmse,
    MAE = mae
  )
}

# Variables used for model fitting
feature_vars <- c(
  "Year",
  "Rainfall",
  "Temperature",
  "Longitude",
  "Latitude",
  "Total",
  "Order",
  "Family",
  "Diet",
  "Habitat",
  "Sample_Type",
  "Methods",
  "Continent",
  "Country"
)

# ----------------------------------------------------------
# 5. Spatial block CV loop
# ----------------------------------------------------------

spatial_cv_results <- list()
spatial_cv_predictions <- list()
spatial_cv_importance <- list()

for (i in 1:k) {
  
  message("Running spatial block fold ", i, " of ", k)
  
  train_df <- ml_data %>%
    dplyr::filter(fold != i)
  
  test_df <- ml_data %>%
    dplyr::filter(fold == i)
  
  if (nrow(test_df) < 10 || nrow(train_df) < 30) {
    warning("Fold ", i, " skipped because train or test set is too small.")
    next
  }
  
  train_model_df <- train_df %>%
    dplyr::select(Prevalence, all_of(feature_vars))
  
  test_model_df <- test_df %>%
    dplyr::select(Prevalence, all_of(feature_vars))
  
  # Convert character variables to factors
  train_model_df <- train_model_df %>%
    dplyr::mutate(across(where(is.character), as.factor))
  
  test_model_df <- test_model_df %>%
    dplyr::mutate(across(where(is.character), as.factor))
  
  # Match test factor levels to training levels
  for (v in names(train_model_df)) {
    if (is.factor(train_model_df[[v]]) && v %in% names(test_model_df)) {
      test_model_df[[v]] <- factor(
        as.character(test_model_df[[v]]),
        levels = levels(train_model_df[[v]])
      )
    }
  }
  
  # Impute numerical variables using training medians
  for (v in names(train_model_df)) {
    if (is.numeric(train_model_df[[v]])) {
      med <- median(train_model_df[[v]], na.rm = TRUE)
      
      train_model_df[[v]][is.na(train_model_df[[v]])] <- med
      test_model_df[[v]][is.na(test_model_df[[v]])] <- med
    }
  }
  
  # Impute categorical variables using training modes
  for (v in names(train_model_df)) {
    if (is.factor(train_model_df[[v]])) {
      mode_v <- get_mode(train_model_df[[v]])
      
      train_model_df[[v]][is.na(train_model_df[[v]])] <- mode_v
      test_model_df[[v]][is.na(test_model_df[[v]])] <- mode_v
      
      test_model_df[[v]] <- factor(
        as.character(test_model_df[[v]]),
        levels = levels(train_model_df[[v]])
      )
    }
  }
  
  dummies_fold <- caret::dummyVars(
    Prevalence ~ .,
    data = train_model_df,
    fullRank = FALSE
  )
  
  x_train <- predict(dummies_fold, newdata = train_model_df) %>%
    as.data.frame()
  
  x_test <- predict(dummies_fold, newdata = test_model_df) %>%
    as.data.frame()
  
  # Remove response column if it appears
  x_train <- x_train[, !colnames(x_train) %in% "Prevalence", drop = FALSE]
  x_test <- x_test[, !colnames(x_test) %in% "Prevalence", drop = FALSE]
  
  # Align columns in test to training
  missing_cols <- setdiff(colnames(x_train), colnames(x_test))
  
  if (length(missing_cols) > 0) {
    for (mc in missing_cols) {
      x_test[[mc]] <- 0
    }
  }
  
  extra_cols <- setdiff(colnames(x_test), colnames(x_train))
  
  if (length(extra_cols) > 0) {
    x_test <- x_test[, setdiff(colnames(x_test), extra_cols), drop = FALSE]
  }
  
  x_test <- x_test[, colnames(x_train), drop = FALSE]
  

  nzv_fold <- caret::nearZeroVar(x_train)
  
  if (length(nzv_fold) > 0) {
    x_train <- x_train[, -nzv_fold, drop = FALSE]
    x_test <- x_test[, colnames(x_train), drop = FALSE]
  }
  
  preproc_fold <- caret::preProcess(
    x_train,
    method = c("center", "scale")
  )
  
  x_train_pp <- predict(preproc_fold, x_train)
  x_test_pp <- predict(preproc_fold, x_test)
  
  y_train <- train_model_df$Prevalence
  y_test <- test_model_df$Prevalence
  
  rf_fold <- randomForest::randomForest(
    x = x_train_pp,
    y = y_train,
    ntree = 500,
    importance = TRUE
  )
  
  pred_test <- predict(rf_fold, newdata = x_test_pp)
  pred_test <- pmin(pmax(pred_test, 0), 1)
  
  metrics <- calc_metrics(y_test, pred_test) %>%
    dplyr::mutate(
      fold = i,
      n_train = nrow(train_df),
      n_test = nrow(test_df)
    )
  
  spatial_cv_results[[i]] <- metrics
  
  spatial_cv_predictions[[i]] <- test_df %>%
    dplyr::select(
      Longitude,
      Latitude,
      Prevalence,
      spatial_block,
      fold,
      Continent,
      Country
    ) %>%
    dplyr::mutate(
      Predicted = pred_test,
      Residual = Prevalence - Predicted
    )
  
  imp <- randomForest::importance(rf_fold) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Variable") %>%
    dplyr::mutate(fold = i)
  
  spatial_cv_importance[[i]] <- imp
}

spatial_cv_metrics <- dplyr::bind_rows(spatial_cv_results)

spatial_cv_summary <- spatial_cv_metrics %>%
  dplyr::summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_RMSE = sd(RMSE, na.rm = TRUE),
    mean_MAE = mean(MAE, na.rm = TRUE),
    sd_MAE = sd(MAE, na.rm = TRUE),
    total_test_n = sum(n_test, na.rm = TRUE)
  )

spatial_cv_pred <- dplyr::bind_rows(spatial_cv_predictions)
spatial_cv_imp <- dplyr::bind_rows(spatial_cv_importance)

print(spatial_cv_metrics)
print(spatial_cv_summary)

# Save outputs
write.csv(
  spatial_cv_metrics,
  "Spatial_block_CV_fold_metrics.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  spatial_cv_summary,
  "Spatial_block_CV_summary.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  spatial_cv_pred,
  "Spatial_block_CV_predictions.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  spatial_cv_imp,
  "Spatial_block_CV_variable_importance.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

library(dplyr)
library(caret)
library(randomForest)
library(tibble)

continents <- unique(ml_data$Continent)
continents <- continents[!is.na(continents)]

continent_cv_results <- list()

for (ct in continents) {
  
  message("Leaving out continent: ", ct)
  
  train_df <- ml_data %>%
    dplyr::filter(Continent != ct)
  
  test_df <- ml_data %>%
    dplyr::filter(Continent == ct)
  
  if (nrow(test_df) < 10 || nrow(train_df) < 30) {
    warning("Skipping ", ct, " because sample size is too small.")
    next
  }
  
  train_model_df <- train_df %>%
    dplyr::select(Prevalence, all_of(feature_vars))
  
  test_model_df <- test_df %>%
    dplyr::select(Prevalence, all_of(feature_vars))
  
  train_model_df <- train_model_df %>%
    dplyr::mutate(across(where(is.character), as.factor))
  
  test_model_df <- test_model_df %>%
    dplyr::mutate(across(where(is.character), as.factor))
  
  for (v in names(train_model_df)) {
    if (is.factor(train_model_df[[v]]) && v %in% names(test_model_df)) {
      test_model_df[[v]] <- factor(
        as.character(test_model_df[[v]]),
        levels = levels(train_model_df[[v]])
      )
    }
  }
  
 
  for (v in names(train_model_df)) {
    if (is.numeric(train_model_df[[v]])) {
      med <- median(train_model_df[[v]], na.rm = TRUE)
      train_model_df[[v]][is.na(train_model_df[[v]])] <- med
      test_model_df[[v]][is.na(test_model_df[[v]])] <- med
    }
  }
  

  get_mode <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    names(sort(table(x), decreasing = TRUE))[1]
  }
  
  for (v in names(train_model_df)) {
    if (is.factor(train_model_df[[v]])) {
      mode_v <- get_mode(train_model_df[[v]])
      train_model_df[[v]][is.na(train_model_df[[v]])] <- mode_v
      test_model_df[[v]][is.na(test_model_df[[v]])] <- mode_v
      test_model_df[[v]] <- factor(
        as.character(test_model_df[[v]]),
        levels = levels(train_model_df[[v]])
      )
    }
  }
  
  dummies_ct <- caret::dummyVars(
    Prevalence ~ .,
    data = train_model_df,
    fullRank = FALSE
  )
  
  x_train <- predict(dummies_ct, newdata = train_model_df) %>%
    as.data.frame()
  
  x_test <- predict(dummies_ct, newdata = test_model_df) %>%
    as.data.frame()
  
  x_train <- x_train[, !colnames(x_train) %in% "Prevalence", drop = FALSE]
  x_test <- x_test[, !colnames(x_test) %in% "Prevalence", drop = FALSE]
  

  missing_cols <- setdiff(colnames(x_train), colnames(x_test))
  if (length(missing_cols) > 0) {
    for (mc in missing_cols) {
      x_test[[mc]] <- 0
    }
  }
  
  extra_cols <- setdiff(colnames(x_test), colnames(x_train))
  if (length(extra_cols) > 0) {
    x_test <- x_test[, setdiff(colnames(x_test), extra_cols), drop = FALSE]
  }
  
  x_test <- x_test[, colnames(x_train), drop = FALSE]
  
  nzv_ct <- caret::nearZeroVar(x_train)
  
  if (length(nzv_ct) > 0) {
    x_train <- x_train[, -nzv_ct, drop = FALSE]
    x_test <- x_test[, colnames(x_train), drop = FALSE]
  }
  
  preproc_ct <- caret::preProcess(
    x_train,
    method = c("center", "scale")
  )
  
  x_train_pp <- predict(preproc_ct, x_train)
  x_test_pp <- predict(preproc_ct, x_test)
  
  rf_ct <- randomForest::randomForest(
    x = x_train_pp,
    y = train_model_df$Prevalence,
    ntree = 500,
    importance = TRUE
  )
  
  pred_ct <- predict(rf_ct, newdata = x_test_pp)
  pred_ct <- pmin(pmax(pred_ct, 0), 1)
  
  metrics_ct <- calc_metrics(test_model_df$Prevalence, pred_ct) %>%
    dplyr::mutate(
      left_out_continent = ct,
      n_train = nrow(train_df),
      n_test = nrow(test_df)
    )
  
  continent_cv_results[[ct]] <- metrics_ct
}

continent_cv_metrics <- dplyr::bind_rows(continent_cv_results)

continent_cv_summary <- continent_cv_metrics %>%
  dplyr::summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_RMSE = sd(RMSE, na.rm = TRUE),
    mean_MAE = mean(MAE, na.rm = TRUE),
    sd_MAE = sd(MAE, na.rm = TRUE)
  )

print(continent_cv_metrics)
print(continent_cv_summary)

write.csv(
  continent_cv_metrics,
  "Leave_one_continent_out_CV_metrics.csv",
  row.names = FALSE
)

write.csv(
  continent_cv_summary,
  "Leave_one_continent_out_CV_summary.csv",
  row.names = FALSE
)