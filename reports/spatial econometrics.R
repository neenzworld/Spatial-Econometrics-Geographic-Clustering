# ==============================================================================
# SPATIAL ECONOMETRICS & GEOGRAPHIC CLUSTERING
# Author: Neena Abimbola Albarus
# Description: Cross-sectional spatial analysis of OWID global life expectancy data.
#              Includes spatial joins, contiguity weight matrices, Moran's I,
#              and Local Indicators of Spatial Association (LISA).
# ==============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(spdep)
library(RColorBrewer)
library(viridis)

# ------------------------------------------------------------------------------
# DYNAMIC DATA PATH RESOLUTION
# ------------------------------------------------------------------------------
# Check potential paths so it works whether working directory is root or /reports
possible_paths <- c(
  "data/life_expectancy_data.csv",
  "../data/life_expectancy_data.csv",
  "/Users/yawdgyal/Desktop/R for github/Spatial Econometrics & Geographic Clustering/data/life_expectancy_data.csv"
)

data_path <- possible_paths[file.exists(possible_paths)][1]

if (is.na(data_path)) {
  stop("Could not locate 'life_expectancy_data.csv'. Please verify it exists in your project's data/ folder.")
}

# Ingest raw OWID panel dataset
owid_raw <- read_csv(data_path)

# Standardize column names and filter to 2019 cross-section
owid_2019 <- owid_raw %>%
  rename(
    country  = Entity,
    code     = Code,
    year     = Year,
    life_exp = `Life expectancy`
  ) %>%
  filter(year == 2019, !is.na(code)) %>%
  select(code, country, life_exp)

# Fetch spatial geometries from Natural Earth
world_polygons <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(iso_a3, name, continent, geometry) %>%
  filter(continent != "Antarctica")

# Join attribute data with vector boundaries on ISO 3-letter codes
spatial_df <- world_polygons %>%
  inner_join(owid_2019, by = c("iso_a3" = "code")) %>%
  filter(!st_is_empty(geometry))

# Figure 1: Global Continuous Choropleth Map
fig1 <- ggplot(data = spatial_df) +
  geom_sf(aes(fill = life_exp), color = "white", size = 0.15) +
  scale_fill_viridis_c(
    option = "magma", 
    name = "Life Exp.\n(Years)",
    na.value = "grey80"
  ) +
  theme_minimal() +
  labs(
    title = "Figure 1: Global Distribution of National Life Expectancy (2019)",
    subtitle = "Continuous Spatial Projection Illustrating Global Health Disparities",
    caption = "Source: Our World in Data & Natural Earth"
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.title = element_text(face = "bold", size = 13)
  )

print(fig1)

# Figure 2: Continental Density & Distribution Disparities
fig2 <- ggplot(spatial_df, aes(x = reorder(continent, life_exp, FUN = median), y = life_exp, fill = continent)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.2, outlier.color = "red", alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5, color = "grey20") +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Figure 2: Continental Distribution & Intra-Regional Variance",
    subtitle = "Violin and Boxplot Decompositions of Health Metrics (2019)",
    x = "Continent",
    y = "Life Expectancy at Birth (Years)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13)
  )

print(fig2)

# Spatial Weight Matrix & Spatial Lag Construction
neighbors_list <- poly2nb(spatial_df, queen = TRUE)
spatial_weights <- nb2listw(neighbors_list, style = "W", zero.policy = TRUE)

spatial_df$spatial_lag_life_exp <- lag.listw(spatial_weights, spatial_df$life_exp, zero.policy = TRUE)

spatial_df <- spatial_df %>%
  mutate(
    scale_life_exp = scale(life_exp)[,1],
    scale_lag_life_exp = scale(spatial_lag_life_exp)[,1]
  )

# Figure 3: Global Moran's Scatterplot
fig3 <- ggplot(spatial_df, aes(x = scale_life_exp, y = scale_lag_life_exp)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(color = "steelblue", alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
  theme_minimal() +
  labs(
    title = "Figure 3: Global Moran's Scatterplot of Health Metrics",
    subtitle = "Slope Equals the Global Moran's I Spatial Autocorrelation Coefficient",
    x = "Standardized National Life Expectancy (Z-Score)",
    y = "Standardized Spatial Lag / Neighbor Average (Z-Score)"
  ) +
  theme(plot.title = element_text(face = "bold", size = 13))

print(fig3)

# Global & Local Moran's I Calculations
moran_global <- moran.test(spatial_df$life_exp, listw = spatial_weights, zero.policy = TRUE)
print(moran_global)

lisa_results <- localmoran(spatial_df$life_exp, listw = spatial_weights, zero.policy = TRUE)

spatial_df <- spatial_df %>%
  mutate(
    local_moran_i = lisa_results[, "Ii"],
    lisa_p_value  = lisa_results[, "Pr(z != E(Ii))"],
    quadrant = case_when(
      scale_life_exp > 0 & scale_lag_life_exp > 0 ~ "High-High",
      scale_life_exp < 0 & scale_lag_life_exp < 0 ~ "Low-Low",
      scale_life_exp < 0 & scale_lag_life_exp > 0 ~ "Low-High",
      scale_life_exp > 0 & scale_lag_life_exp < 0 ~ "High-Low",
      TRUE ~ "Not Significant"
    ),
    cluster_category = if_else(lisa_p_value < 0.05, quadrant, "Not Significant")
  )

# Figure 4: LISA Cluster Map
lisa_colors <- c(
  "High-High"       = "#d7191c",
  "Low-Low"         = "#2b83ba",
  "Low-High"        = "#fdae61",
  "High-Low"        = "#abdda4",
  "Not Significant" = "grey85"
)

fig4 <- ggplot(data = spatial_df) +
  geom_sf(aes(fill = cluster_category), color = "white", size = 0.15) +
  scale_fill_manual(values = lisa_colors, name = "LISA Cluster Type\n(p < 0.05)") +
  theme_minimal() +
  labs(
    title = "Figure 4: Local Indicators of Spatial Association (LISA) Cluster Map",
    subtitle = "Statistically Significant Health Hotspots (High-High) & Coldspots (Low-Low)"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13)
  )

print(fig4)

# Figure 5: Continental Distribution of Local Moran Scores
fig5 <- ggplot(spatial_df, aes(x = continent, y = local_moran_i, fill = continent)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "red") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(
    title = "Figure 5: Magnitude of Local Spatial Clustering Across Continents",
    subtitle = "Higher Positive Local Moran Statistics Indicate Stronger Local Spatial Dependence",
    x = "Continent",
    y = "Local Moran's I Statistic (Ii)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 13)
  )

print(fig5)