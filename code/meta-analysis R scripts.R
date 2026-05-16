# 加载必要的包
if (!require("metafor")) install.packages("metafor")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("readxl")) install.packages("readxl")
if (!require("tidyr")) install.packages("tidyr")
if (!require("scales")) install.packages("scales")
if (!require("patchwork")) install.packages("patchwork")

library(metafor)
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(scales)
library(patchwork)
getwd()
setwd("D:/R/R")


results_dir <- file.path(getwd(), "metafor_results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cat("正在读取数据...\n")
data <- read_excel("1.xlsx", sheet = 1, col_names = FALSE)
colnames(data) <- c("Protozoa", "Subgroup_type", "Group", "Event", "Total")


data <- data %>%
  mutate(Protozoa = ifelse(Protozoa == "", NA, Protozoa)) %>%
  fill(Protozoa, .direction = "down") %>%
  filter(Event != "Event" | is.na(Event)) %>%
  filter(Total != "Total" | is.na(Total)) %>%
  mutate(
    Event = as.numeric(Event),
    Total = as.numeric(Total)
  ) %>%
  filter(
    !is.na(Event), !is.na(Total),
    Event <= Total, Total > 0, !is.na(Protozoa)
  )

cat("数据清理完成！\n")
cat("原虫种类：", unique(data$Protozoa), "\n")

run_metafor_analysis <- function(p_data, p_name) {
  
  cat(paste("\n=== 分析：", p_name, " ===\n"))
  cat("样本量：", nrow(p_data), "\n")
  
  if (nrow(p_data) < 2) {
    cat("样本不足，跳过\n")
    return(NULL)
  }
  
  # 创建安全的文件名
  safe_name <- gsub("[^a-zA-Z0-9]", "_", p_name)
  
  # (1) 计算效应量
  dat <- escalc(
    measure = "Loght", 
    xi = Event, 
    ni = Total, 
    data = p_data
  )
  
  res <- tryCatch({
    rma(yi = yi, vi = vi, data = dat, method = "REML")
  }, error = function(e) {
    cat("总体模型拟合失败：", e$message, "\n")
    return(NULL)
  })
  
  res_subgroup <- tryCatch({
    n_subgroups <- length(unique(p_data$Subgroup_type))
    if (n_subgroups >= nrow(p_data)) {
      stop("亚组数量大于等于样本量，无法分析")
    }
    rma(yi = yi, vi = vi, mods = ~ Subgroup_type - 1, data = dat, method = "REML")
  }, error = function(e) {
    cat("亚组分析跳过：", e$message, "\n")
    return(NULL)
  })
  

  plot_data <- p_data %>%
    mutate(
      yi = log(Event / (Total - Event)),
      vi = 1/Event + 1/(Total - Event),
      ci.lb = yi - 1.96 * sqrt(vi),
      ci.ub = yi + 1.96 * sqrt(vi),
      weight = 1 / vi,
      weight = weight / sum(weight) * 100,  # 归一化权重为百分比
      Subgroup = Subgroup_type,
      Group = Group,
      Study = paste(Subgroup_type, Group, sep = " - "),
      slab = paste(Subgroup_type, Group, sep = " - ")
    ) %>%
    arrange(Subgroup, desc(weight))
  
  plot_data <- plot_data %>%
    mutate(
      yi = ifelse(is.infinite(yi), 0, yi),
      vi = ifelse(is.infinite(vi) | is.na(vi), 1, vi)
    )
  
  plot_data <- plot_data %>%
    mutate(
      prop = plogis(yi),
      prop.lb = plogis(ci.lb),
      prop.ub = plogis(ci.ub),
      Prop_CI = paste0(round(prop * 100, 1), "% (", 
                       round(prop.lb * 100, 1), "%, ", 
                       round(prop.ub * 100, 1), "%)"),
      Event_Total = paste0(Event, "/", Total),
      Weight = paste0(round(weight, 1), "%")
    )
  
  plot_data$y_order <- factor(plot_data$Study, levels = rev(plot_data$Study))
  
  p_forest <- ggplot(plot_data, aes(y = y_order, x = prop)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = weight), color = "#2c3e50") +
    geom_errorbar(
      aes(xmin = prop.lb, xmax = prop.ub), 
      width = 0.2, 
      color = "#34495e",
      orientation = "y"
    ) +
    facet_wrap(~ Subgroup, scales = "free_y", ncol = 1) +
    scale_x_continuous(labels = percent_format(accuracy = 1), 
                       limits = c(0, min(max(plot_data$prop.ub, na.rm = TRUE) * 1.1, 1))) +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    labs(
      title = paste(p_name, "Subgroup Analysis Forest Plot"),
      x = "Prevalence (95% CI)",
      y = ""
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#ecf0f1"),
      strip.text = element_text(face = "bold", size = 9),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    )
  
  table_data <- plot_data %>%
    select(y_order, Study, Event_Total, Weight, Prop_CI) %>%
    pivot_longer(cols = -y_order, names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = c("Study", "Event_Total", "Weight", "Prop_CI"),
                             labels = c("Study", "Event/Total", "Weight", "Prevalence (95% CI)")))
  
  p_table <- ggplot(table_data, aes(x = variable, y = y_order)) +
    geom_text(aes(label = value), hjust = 0, size = 3) +
    facet_wrap(~ plot_data$Subgroup[match(y_order, plot_data$y_order)], 
               scales = "free_y", ncol = 1) +
    scale_x_discrete(position = "top", expand = c(0, 0.5)) +
    labs(x = "", y = "") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      plot.margin = margin(5.5, 0, 5.5, 5.5)
    )
  
  p_combined <- p_table + p_forest + plot_layout(widths = c(2, 1.5))
  
  pdf_file <- file.path(results_dir, paste0(safe_name, "_forest_plot.pdf"))
  ggsave(
    filename = pdf_file,
    plot = p_combined,
    width = 14,
    height = 6 + nrow(plot_data) * 0.2,
    device = "pdf"
  )
  
  cat("森林图（可编辑PDF）已保存：", pdf_file, "\n")
  
  result_file <- file.path(results_dir, paste0(safe_name, "_statistics.txt"))
  sink(result_file)
  cat("Protozoan:", p_name, "\n")
  cat("Measure: PLN (Logit transformation)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  if (!is.null(res)) {
    cat("1. Overall Meta-Analysis (Random Effects Model)\n")
    print(res)
    cat("\n\n3. Heterogeneity Statistics:\n")
    cat(paste("I² =", round(res$I2, 2), "%\n"))
    cat(paste("τ² =", round(res$tau2, 4), "\n"))
    cat(paste("Q-test p-value =", round(res$pval, 6), "\n"))
  } else {
    cat("1. Overall Meta-Analysis: Failed to fit model\n")
  }
  
  if (!is.null(res_subgroup)) {
    cat("\n\n2. Subgroup Analysis (Meta-Regression)\n")
    print(res_subgroup)
  } else {
    cat("\n\n2. Subgroup Analysis: Not performed (insufficient data)\n")
  }
  
  sink()
  cat("统计结果已保存：", result_file, "\n")
  
  return(list(model = res, subgroup = res_subgroup, plot = p_combined))
}

protozoa_list <- unique(data$Protozoa)
all_results <- list()

for (p in protozoa_list) {
  p_data <- data %>% filter(Protozoa == p)
  all_results[[p]] <- run_metafor_analysis(p_data, p)
}

cat("\n===== 全部分析完成！=====\n")
cat("结果保存在：", results_dir, "\n")
cat("生成的文件：\n")
print(list.files(results_dir))
