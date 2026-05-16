
pacman::p_load(
  readxl, dplyr, tidyr, stringr,
  terra, sf,
  ggplot2, ggspatial, viridis, patchwork,
  maps
)

sf::sf_use_s2(FALSE)

normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
period_labels <- c(
  "Pre-2000" = "2000年以前",
  "2000-2010" = "2000-2010年",
  "Post-2010" = "2010年以后"
)

df_raw <- read_excel("zon.xlsx", sheet = 1)
df_clean <- df_raw %>%
  select(
    Author, Year, Event = Event...3, Total = Total...4,
    Species, Order, Family, Longitude, Latitude
  ) %>%
  drop_na(Longitude, Latitude, Total, Species, Order) %>%
  filter(Total > 0, Total >= 5) %>%
  mutate(
    Positive_Rate = Event / Total,
    Period = case_when(
      Year < 2000 ~ "Pre-2000",
      Year >= 2000 & Year <= 2010 ~ "2000-2010",
      Year > 2010 ~ "Post-2010"
    ),
    Species = str_squish(Species) %>% str_to_sentence(),
    Order = str_squish(Order) %>% str_to_sentence()
  )

aoh_folder <- "D:/R/R/zonation_output/10974868/"
aoh_files_all <- list.files(aoh_folder, pattern = "\\.rds$", full.names = TRUE)
aoh_filenames <- basename(aoh_files_all)
aoh_orders <- str_remove(aoh_filenames, "\\.rds$") %>% str_squish()

target_orders <- unique(df_clean$Order)
matched_indices <- which(aoh_orders %in% target_orders)
target_aoh_files <- aoh_files_all[matched_indices]
target_aoh_orders <- aoh_orders[matched_indices]

aoh_list <- list()
for (i in 1:length(target_aoh_files)) {
  aoh_sv <- readRDS(target_aoh_files[i])
  aoh_sv_subset <- aoh_sv[, "sciname", drop = FALSE]
  values(aoh_sv_subset)$sciname <- values(aoh_sv_subset)$sciname %>%
    str_squish() %>% 
    str_to_sentence() %>% 
    str_remove_all("\\.|\\?|\\*|\\d+") %>%  # 去除特殊符号/数字
    str_trim()
  names(aoh_sv_subset) <- "Species"
  aoh_list[[target_aoh_orders[i]]] <- aoh_sv_subset
}

aoh_combined <- aoh_list[[1]]
if (length(aoh_list) > 1) {
  for (i in 2:length(aoh_list)) {
    aoh_combined <- rbind(aoh_combined, aoh_list[[i]])
  }
}

target_species <- unique(df_clean$Species)
aoh_species <- values(aoh_combined)$Species
matched_species_indices <- which(aoh_species %in% target_species)
aoh_matched <- aoh_combined[matched_species_indices, ]

global_template <- rast(
  xmin = -180, xmax = 180, ymin = -60, ymax = 90,
  resolution = 0.1, crs = "EPSG:4326"
)
values(global_template) <- 1

matched_species_unique <- unique(values(aoh_matched)$Species)
species_rasters <- list()
for (sp in matched_species_unique) {
  sp_indices <- which(values(aoh_matched)$Species == sp)
  sp_sv <- aoh_matched[sp_indices, ]
  sp_raster <- rasterize(sp_sv, global_template, field = 1, background = NA, touches = TRUE)
  species_rasters[[sp]] <- sp_raster
}
aoh_stack <- rast(species_rasters)

period_species <- df_clean %>%
  filter(Event > 0) %>%
  group_by(Period) %>%
  summarise(Species_List = list(unique(Species)), .groups = "drop")

period_risk_rasters <- list()
period_labels <- c("Pre-2000"="2000年以前", "2000-2010"="2000-2010年", "Post-2010"="2010年以后")

for (i in 1:nrow(period_species)) {
  p <- period_species$Period[i]
  sp_list <- period_species$Species_List[[i]]
  sp_in_aoh <- intersect(sp_list, names(aoh_stack))
  
  if (length(sp_in_aoh) > 0) {
    period_stack <- aoh_stack[[sp_in_aoh]]
    risk_raster <- sum(period_stack, na.rm = TRUE)
    period_risk_rasters[[p]] <- risk_raster
  }
}

world_map <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>%
  st_transform(crs = "EPSG:4326")

if(length(period_risk_rasters) == 0){
  stop("❌ 没有生成任何风险栅格，请检查AOH数据与物种匹配！")
}
all_risk_values <- unlist(lapply(period_risk_rasters, function(x) na.omit(values(x))))
max_risk <- ceiling(max(all_risk_values, na.rm=TRUE))
min_risk <- 0



plot_dfs <- list()
for (p in names(period_risk_rasters)) {
  df <- as.data.frame(period_risk_rasters[[p]], xy = TRUE) %>%
    drop_na()
  colnames(df)[3] <- "Positive_Species_Count"
  df$Period <- period_labels[p]
  plot_dfs[[p]] <- df
}

plot_risk_map <- function(plot_df, period_name, max_risk, world_map) {
  ggplot() +
    geom_rect(aes(xmin=-180, xmax=180, ymin=-60, ymax=90), fill="#d0e8f2", color=NA) +
    geom_tile(data = plot_df, aes(x = x, y = y, fill = Positive_Species_Count)) +
    geom_sf(data = world_map, fill = NA, color = "#000000", linewidth = 0.4) +
    scale_fill_viridis_c(
      option = "inferno",        
      limits = c(min_risk, max_risk),  
      breaks = seq(min_risk, max_risk, by=1),
      name = NULL,               
      na.value = "transparent",
      guide = guide_colorbar(
        direction = "vertical",
        barwidth = unit(0.5, "cm"),
        barheight = unit(10, "cm"),
        label.position = "right",
        label.hjust = 0,
        ticks = TRUE,
        ticks.colour = "black",
        ticks.linewidth = 1
      )
    ) +
    labs(title = period_name, subtitle = "全球野生哺乳动物弓形虫患病风险地图") +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 90), expand = FALSE) +
    theme_minimal(base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(face = "bold", size = 12, colour = "black"),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    annotation_scale(location = "bl") +
    annotation_north_arrow(location = "tr")
}
cat("\n正在生成PDF地图...\n")
for (p in names(period_risk_rasters)) {
  final_map <- plot_risk_map(plot_dfs[[p]], period_labels[[p]], max_risk, world_map)
  ggsave(paste0("弓形虫风险地图_", p, ".pdf"), final_map, width=16, height=9, dpi=600, bg="white")
  cat("✅ 已保存：", p, "地图\n")
}

combined_map <- wrap_plots(
  plot_risk_map(plot_dfs[["Pre-2000"]], period_labels[["Pre-2000"]], max_risk, world_map),
  plot_risk_map(plot_dfs[["2000-2010"]], period_labels[["2000-2010"]], max_risk, world_map),
  plot_risk_map(plot_dfs[["Post-2010"]], period_labels[["Post-2010"]], max_risk, world_map),
  ncol=3
) + plot_annotation(title = "全球野生哺乳动物弓形虫风险分时段对比") & theme(legend.position="none")

ggsave("弓形虫风险地图_三时段合并版.pdf", combined_map, width=24, height=8, dpi=1200, bg="white")

