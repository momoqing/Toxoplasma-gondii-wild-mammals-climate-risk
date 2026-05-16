
library(readxl)
library(dplyr)
library(caret)
library(randomForest)
library(xgboost)
library(e1071)
library(glmnet)
library(rpart)
library(ggplot2)

set.seed(1234)
options(scipen = 999)

data <- read_excel("ml.xlsx", sheet = "Sheet1")
model_data <- data %>%
  select(Prevalence, Year, Rainfall, Temperature, Longitude, Latitude, Total,
         Order, Family, Diet, Habitat, Sample_Type, Methods, Continent, Country) %>%
  mutate(across(where(is.numeric), ~replace_na(., median(., na.rm = TRUE))),
         across(where(is.character), ~replace_na(., names(which.max(table(.))))),
         across(where(is.character), as.factor)) %>%
  filter(Prevalence >= 0 & Prevalence <= 1)

train_index <- createDataPartition(model_data$Prevalence, p = 0.75, list = FALSE)
train_data <- model_data[train_index, ]
test_data  <- m
dummies <- dummyVars(Prevalence ~ ., data = train_data, fullRank = TRUE)
x_train_dum <- predict(dummies, train_data)
x_test_dum  <- predict(dummies, test_data)

nzv <- nearZeroVar(x_train_dum)
if(length(nzv) > 0){
  x_train_dum <- x_train_dum[, -nzv]
  x_test_dum  <- x_test_dum[, -nzv]
}

preProc <- preProcess(x_train_dum, method = c("center", "scale"))
x_train <- predict(preProc, x_train_dum)
x_test  <- predict(preProc, x_test_dum)
f
y_train <- train_data$Prevalence
y_test  <- test_data$Prevalence

ctrl <- trainControl(method = "cv", number = 10, verboseIter = FALSE)
metric_fun <- function(actual, pred) {
  data.frame(
    RMSE = RMSE(pred, actual),
    MAE = MAE(pred, actual),
    R2 = R2(pred, actual)
  )
}
model_results <- list()

cat("\n=== 训练 线性回归 模型 ===\n")
set.seed(1234)
lm_model <- train(x_train, y_train, method = "lm", trControl = ctrl, metric = "RMSE")
lm_pred <- predict(lm_model, x_test)
lm_metrics <- metric_fun(y_test, lm_pred)
model_results[["Linear_Regression"]] <- lm_metrics

cat("\n=== 训练 KNN 模型 ===\n")
set.seed(1234)
knn_grid <- expand.grid(k = seq(3, 21, 2))
knn_model <- train(x_train, y_train, method = "knn", trControl = ctrl, tuneGrid = knn_grid, metric = "RMSE")
knn_pred <- predict(knn_model, x_test)
knn_metrics <- metric_fun(y_test, knn_pred)
model_results[["KNN"]] <- knn_metrics

cat("\n=== 训练 SVM 模型 ===\n")
svm_grid <- expand.grid(
  sigma = c(0.02, 0.05, 0.08),  # 核函数宽度
  C = c(0.8, 1, 1.5)           # 正则化强度
)
set.seed(1234)
svm_model <- train(
  x = x_train,
  y = y_train,
  method = "svmRadial",
  trControl = ctrl,
  tuneGrid = svm_grid,
  metric = "RMSE"
)
svm_pred <- predict(svm_model, x_test)
svm_metrics <- metric_fun(y_test, svm_pred)
model_results[["SVM"]] <- svm_metrics

cat("\n=== 训练 Random Forest 模型 ===\n")
rf_grid <- expand.grid(mtry = seq(5, min(30, length(feature_names)), by = 5))
set.seed(1234)
rf_model <- train(
  x = x_train,
  y = y_train,
  method = "rf",
  trControl = ctrl,
  tuneGrid = rf_grid,
  metric = "RMSE",
  ntree = 500,         
  maxnodes = 15,        
  importance = TRUE
)
rf_pred <- predict(rf_model, x_test)
rf_metrics <- metric_fun(y_test, rf_pred)
model_results[["Random_Forest"]] <- rf_metrics

cat("\n=== 训练 XGBoost ===")
xgb_model <- xgboost(
  data = x_train,
  label = y_train,
  nrounds = 100,
  max_depth = 3,
  eta = 0.1,
  objective = "reg:squarederror",
  verbose = 0
)
xgb_pred <- predict(xgb_model, x_test)
model_results[["XGBoost"]] <- metric_fun(y_test, xgb_pred)

cat("\n=== 训练 决策树 模型 ===\n")
set.seed(1234)
dt_grid <- expand.grid(cp = seq(0.001,0.05,0.005))
dt_model <- train(x_train, y_train, method = "rpart", trControl = ctrl, tuneGrid = dt_grid, metric = "RMSE")
dt_pred <- predict(dt_model, x_test)
dt_metrics <- metric_fun(y_test, dt_pred)
model_results[["Decision_Tree"]] <- dt_metrics

# 4.7 GLMNET
cat("\n=== 训练 GLMNET 模型 ===\n")
set.seed(1234)
glm_grid <- expand.grid(
  alpha = c(0.3, 0.5, 0.7),   # 弹性网络最优区间
  lambda = c(0.005, 0.01, 0.02)
)
glm_model <- train(as.matrix(x_train), y_train, method = "glmnet", trControl = ctrl, tuneGrid = glm_grid, metric = "RMSE")
glm_pred <- predict(glm_model, x_test)
glm_metrics <- metric_fun(y_test, glm_pred)
model_results[["GLMNET"]] <- glm_metrics

plot_df <- rbind(
  data.frame(Observed = y_train, Predicted = rf_train_pred, Group = "Training"),
  data.frame(Observed = y_test, Predicted = rf_pred, Group = "Test")
)

fig_A <- ggplot(plot_df, aes(x = Observed, y = Predicted)) +
  geom_point(aes(shape = Group, color = Group), size = 2.5) +
  scale_shape_manual(values = c(15, 16)) +  # 15=方形(训练) 16=圆点(测试)
  scale_color_manual(values = c("#2E86AB", "#E63946")) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1) +
  labs(
    title = "(A) Measured versus predicted effect size using the RF model",
    x = "Measured Prevalence", y = "Predicted Prevalence",
    caption = paste0(
      "Train: R²=", round(rf_train_metrics$R2,3), " RMSE=", round(rf_train_metrics$RMSE,3), " MAE=", round(rf_train_metrics$MAE,3),
      "\nTest:  R²=", round(rf_metrics$R2,3), " RMSE=", round(rf_metrics$RMSE,3), " MAE=", round(rf_metrics$MAE,3)
    )
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.caption = element_text(hjust = 0))

final_results <- bind_rows(model_results, .id = "Model") %>% arrange(desc(R2))
fig_S3 <- ggplot(final_results, aes(x = reorder(Model, -R2), y = R2, fill = Model)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = round(R2,3)), vjust = -0.3) +
  labs(title = "Fig. S3 Performance of machine learning models",
       x = "Model", y = expression(R^2)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave("FigA_RF_Performance.pdf", fig_A, width = 10, height = 7, 
       device = "pdf", dpi = 300, useDingbats = FALSE)
ggsave("FigS3_Model_Comparison.pdf", fig_S3, width = 11, height = 6, 
       device = "pdf", dpi = 300, useDingbats = FALSE)


ggsave("FigA_RF_Performance.pdf", fig_A, width = 10, height = 7, 
       device = "pdf", dpi = 300, useDingbats = FALSE)
ggsave("FigS3_Model_Comparison.pdf", fig_S3, width = 11, height = 6, 
       device = "pdf", dpi = 300, useDingbats = FALSE)

saveRDS(rf_model, "Best_RF_Model.rds")
