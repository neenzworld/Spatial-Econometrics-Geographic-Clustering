# Spatial Econometrics & Geographic Clustering

**Cross-Sectional Spatial Analysis of Global Health Disparities, Contiguity Weight Matrices and LISA Clusters**

Author: **Neena Abimbola Albarus**

---

## Overview

This workbook applies cross-sectional spatial econometrics to global health data, using Queen contiguity weight matrices, spatial lag construction, Global Moran's I tests, and Local Indicators of Spatial Association (LISA) to detect statistically significant life expectancy hotspots, coldspots, and spatial outliers. It bridges choropleth cartography with formal spatial dependence diagnostics to demonstrate why OLS independence assumptions fail in comparative public health research.

---

## Prerequisites

Before attempting this workbook, you should already know how to:

- Perform data wrangling with **`dplyr`** and **`tidyr`**
- Build layered visualizations with **`ggplot2`**
- Join attribute tables to **`sf`** geometry objects
- Interpret basic **OLS regression output** (coefficients, residuals, t-statistics)

---

## Learning Objectives

By completing this workbook, you will develop competence in:

1. Executing spatial joins between cross-sectional attribute data and vector boundary polygons using `sf` and `rnaturalearth`
2. Constructing spatial neighbor linkages and row-standardized spatial weight matrices (**W**) using `spdep`
3. Calculating spatial lag variables (**Wy**) to systematically quantify neighborhood context effects
4. Running **Global Moran's I** hypothesis tests to evaluate global spatial autocorrelation
5. Deriving and mapping **Local Indicators of Spatial Association (LISA)** to detect statistically significant health hotspots, coldspots, and spatial outliers

---

## Repository Structure

```
.
├── README.md
├── spatial_econometrics.qmd        # Main Quarto workbook
├── data/
│   └── life_expectancy_data.csv    # OWID cross-sectional life expectancy (2019)
└── reports/                        # Rendered HTML output (optional)
    └── spatial_econometrics.html
```

---

## Data Sources

| Source | Description | Link |
|---|---|---|
| **Our World in Data (OWID)** | Cross-sectional national life expectancy records (2019) | https://ourworldindata.org/life-expectancy |
| **Natural Earth** | Medium-scale vector country boundary polygons (ISO A3 codes) | https://www.naturalearthdata.com/ |

The workbook joins these two sources on ISO 3-letter country codes (`iso_a3`) to produce a spatial data frame suitable for spatial econometric analysis.

---

## Required R Packages

```r
install.packages(c(
  "tidyverse",
  "sf",
  "rnaturalearth",
  "spdep",
  "RColorBrewer",
  "viridis"
))
```

> **Note:** `sf` may require system-level geospatial libraries (GDAL, GEOS, PROJ). On macOS, install via Homebrew (`brew install gdal geos proj`) or use the `r2u`/`pak` binary mirrors on Linux.

---

## Workflow Summary

| Module | Purpose |
|---|---|
| **Module 1** | Dynamic path resolution, OWID ingestion, and spatial join with Natural Earth polygons |
| **Module 2** | Choropleth mapping and continental distributional visualizations |
| **Module 3** | Queen contiguity weight matrix construction and spatial lag computation |
| **Module 4** | Global Moran's I test and LISA cluster diagnostics |
| **Lab Exercises** | KNN sensitivity analysis, distance-threshold weights, and spatial outlier extraction |

---

## How to Render

From the project root:

```bash
quarto render spatial_econometrics.qmd
```

Or open the `.qmd` file in **RStudio** and click **Render**.

The workbook expects `life_expectancy_data.csv` to be located in either `data/` or `../data/` relative to the working directory. A hardcoded fallback path is included for local development but should be removed or generalized for distribution.

---

## Key Concepts

- **Tobler's First Law of Geography** — *"Everything is related to everything else, but near things are more related than distant things."*
- **Queen Contiguity** — Neighbor definition based on shared boundary edges *or* vertices
- **Row Standardization (`style = "W"`)** — Each row of the spatial weight matrix sums to 1, so the spatial lag becomes a weighted average of neighbors
- **Global Moran's I** — A single scalar index of overall spatial autocorrelation ($I > 0$ clustering, $I < 0$ dispersion, $I \approx 0$ randomness)
- **LISA Quadrants** — High-High (hotspot), Low-Low (coldspot), High-Low and Low-High (spatial outliers)

---

## License

For educational use. Please cite Our World in Data and Natural Earth as original data providers.
