# Libraries needed
##### Set Working Directory
setwd("C:\\Users\\sasin\\Sync\\Sasindu\\Trent Metabolism\\Lake Ontario\\2024_chl study")

##### Upload the data file
raw_data <- read_csv("raw_data.csv")


library(tidyverse)
library(lubridate)
library(viridis)
library(patchwork)
library(oce)  # For oceanographic calculations

# ============================================================
# 1. FIX DATETIME ISSUES
# ============================================================

cat("Step 1: Fixing datetime issues...\n")

# Examine the datetime problem
cat("\nOriginal datetime format (first 10 rows):\n")
raw_data %>% 
  select(contains("datetime"), contains("epoch")) %>% 
  head(10) %>%
  print()

# Fix datetime using unix epoch (reliable source)
datetime_fixed <- raw_data %>%
  mutate(
    # Convert unix epoch to proper datetime
    datetime_utc = as_datetime(raw_datetime_unixepoch_utc, tz = "UTC"),
    sci_datetime_utc = as_datetime(raw_sci_datetime_unixepoch_utc, tz = "UTC"),
    ctd_datetime_utc = as_datetime(raw_ctd_timestamp_unixepoch_utc, tz = "UTC")
  ) %>%
  arrange(datetime_utc) %>%
  distinct(datetime_utc, .keep_all = TRUE)  # Remove duplicate timestamps

# Check for non-monotonic timestamps (time going backwards)
cat("\nChecking for non-monotonic timestamps...\n")
time_diffs <- diff(datetime_fixed$raw_datetime_unixepoch_utc)
non_monotonic <- which(time_diffs <= 0)

if(length(non_monotonic) > 0) {
  cat("WARNING: Found", length(non_monotonic), "non-monotonic time points\n")
  cat("First 20 locations:", head(non_monotonic, 20), "\n\n")
  
  # Show examples of problematic rows
  cat("Example of non-monotonic timestamps:\n")
  problem_rows <- datetime_fixed[c(non_monotonic[1]-1, non_monotonic[1], non_monotonic[1]+1), ]
  print(problem_rows %>% select(datetime_utc, raw_datetime_unixepoch_utc, depth.m))
} else {
  cat("✓ All timestamps are monotonic (no time reversals)\n")
}

# Summary of temporal coverage
cat("\nTemporal coverage:\n")
cat("Start:", format(min(datetime_fixed$datetime_utc, na.rm = TRUE)), "\n")
cat("End:", format(max(datetime_fixed$datetime_utc, na.rm = TRUE)), "\n")
cat("Duration:", difftime(max(datetime_fixed$datetime_utc, na.rm = TRUE), 
                          min(datetime_fixed$datetime_utc, na.rm = TRUE), 
                          units = "days"), "\n")

# ============================================================
# 2. CHLOROPHYLL & SENSOR CALIBRATION
# ============================================================

cat("\nStep 2: Applying sensor calibrations...\n")

# CALIBRATION PARAMETERS
# TODO: Replace these with your actual calibration values from instrument documentation
# These should be from the FLBBCD calibration sheet

# Chlorophyll fluorescence (FLBBCD)
chl_dark_counts <- 46        # Dark count offset
chl_scale_factor <- 0.0073   # μg/L per count

# CDOM fluorescence
cdom_dark_counts <- 50       # Dark count offset
cdom_scale_factor <- 0.091   # ppb per count

# Optical backscatter (700nm)
bb_dark_counts <- 50         # Dark count offset
bb_scale_factor <- 2.076e-06 # m^-1 sr^-1 per count

cat("\nCalibration parameters:\n")
cat("Chlorophyll: dark =", chl_dark_counts, ", scale =", chl_scale_factor, "μg/L per count\n")
cat("CDOM: dark =", cdom_dark_counts, ", scale =", cdom_scale_factor, "ppb per count\n")
cat("Backscatter: dark =", bb_dark_counts, ", scale =", bb_scale_factor, "m^-1 sr^-1 per count\n")

# ============================================================
# 3. OXYGEN SOLUBILITY FOR FRESHWATER (Lake Ontario)
# ============================================================

# Oxygen solubility function for freshwater (Benson & Krause approximation)
# Highly accurate for 0–40°C at 1 atm, suitable for Lake Ontario
o2_solubility_fresh_mgL <- function(temp_C) {
  # USGS / Benson-Krause based polynomial
  exp(
    -139.34411 + 
      (1.575701e5 / (temp_C + 273.15)) - 
      (6.642308e7 / (temp_C + 273.15)^2) + 
      (1.243800e10 / (temp_C + 273.15)^3) - 
      (8.621949e11 / (temp_C + 273.15)^4)
  )
}

# ============================================================
# 4. APPLY ALL CALIBRATIONS AND CONVERSIONS
# ============================================================

calibrated_data <- datetime_fixed %>%
  mutate(
    # ---------------------------------------------------------
    # Convert lat/lon from DMM (Degrees Decimal Minutes) to DD (Decimal Degrees)
    # Format: DDDMM.MMMM -> DDD + MM.MMMM/60
    # ---------------------------------------------------------
    lat_dd = floor(abs(lat.DMM) / 100) + (abs(lat.DMM) %% 100) / 60,
    lon_dd = -(floor(abs(lon.DMM) / 100) + (abs(lon.DMM) %% 100) / 60),  # Negative for Western Hemisphere
    # Handle NA propagation
    lat_dd = if_else(is.na(lat.DMM), NA_real_, lat_dd),
    lon_dd = if_else(is.na(lon.DMM), NA_real_, lon_dd),
    
    # ---------------------------------------------------------
    # Chlorophyll-a concentration (μg/L)
    # ---------------------------------------------------------
    chlorophyll_ug_L = chl_scale_factor * (flbbcd_chlor_sig.nodim - chl_dark_counts),
    # Set negative values to NA (physical impossibility)
    chlorophyll_ug_L = if_else(chlorophyll_ug_L < 0, NA_real_, chlorophyll_ug_L),
    
    # ---------------------------------------------------------
    # CDOM (Colored Dissolved Organic Matter) in ppb
    # ---------------------------------------------------------
    cdom_ppb = cdom_scale_factor * (flbbcd_cdom_sig.nodim - cdom_dark_counts),
    cdom_ppb = if_else(cdom_ppb < 0, NA_real_, cdom_ppb),
    
    # ---------------------------------------------------------
    # Optical backscatter coefficient (m^-1 sr^-1)
    # ---------------------------------------------------------
    backscatter_m_sr = bb_scale_factor * (flbbcd_bb_sig.nodim - bb_dark_counts),
    backscatter_m_sr = if_else(backscatter_m_sr < 0, NA_real_, backscatter_m_sr),
    
    # ---------------------------------------------------------
    # Salinity (PSU) using Practical Salinity Scale
    # Lake Ontario typically 0.1-0.3 PSU, but we calculate anyway
    # ---------------------------------------------------------
    salinity_PSU = swSCTp(
      conductivity = sci_water_cond.S.m,
      temperature = sci_water_temp.degC,
      pressure = sci_water_pressure.bar * 10,  # Convert bar to dbar (1 bar = 10 dbar)
      conductivityUnit = "S/m"
    ),
    
    # ---------------------------------------------------------
    # Density (kg/m³) from T, S, P
    # ---------------------------------------------------------
    density_kg_m3 = swRho(
      salinity = salinity_PSU,
      temperature = sci_water_temp.degC,
      pressure = sci_water_pressure.bar * 10
    ),
    
    # ---------------------------------------------------------
    # Calculate depth from pressure (more accurate than depth.m)
    # Uses latitude for gravity correction
    # ---------------------------------------------------------
    depth_from_pressure_m = swDepth(
      pressure = sci_water_pressure.bar * 10,
      latitude = lat_dd  # Use converted decimal degrees
    ),
    
    # Use pressure-derived depth if available, otherwise fall back to depth.m
    depth_best_m = coalesce(depth_from_pressure_m, depth.m),
    
    # ---------------------------------------------------------
    # Dissolved Oxygen: Convert from % saturation to mg/L
    # Using freshwater approximation (appropriate for Lake Ontario)
    # ---------------------------------------------------------
    # Calculate saturation concentration at given temperature
    O2_sat_mg_L = o2_solubility_fresh_mgL(rinkoii_temp.degC),
    # Convert % saturation to actual concentration
    DO_mg_L = O2_sat_mg_L * (rinkoii_DO_perc.sat / 100),
    # Flag invalid values
    DO_mg_L = if_else(
      rinkoii_DO_perc.sat <= 0 | 
        is.na(rinkoii_temp.degC) | 
        rinkoii_temp.degC < -2 | 
        rinkoii_temp.degC > 40,
      NA_real_, 
      DO_mg_L
    ),
    
    # ---------------------------------------------------------
    # Add temporal variables for analysis
    # ---------------------------------------------------------
    date = as_date(datetime_utc),
    day = date,
    hour = hour(datetime_utc),
    day_of_year = yday(datetime_utc),
    
    # ---------------------------------------------------------
    # Create depth bins for stratified analysis
    # ---------------------------------------------------------
    depth_bin = cut(
      depth_best_m, 
      breaks = c(0, 5, 10, 20, 30, 50, 100, Inf),
      labels = c("0-5m", "5-10m", "10-20m", "20-30m", "30-50m", "50-100m", ">100m"),
      include.lowest = TRUE
    )
  )

# ============================================================
# 5. DATA QUALITY CHECKS AND SUMMARY
# ============================================================

cat("\nStep 3: Data quality summary...\n\n")

# Count valid observations for key variables
cat("Valid observations:\n")
cat("  Total rows:", nrow(calibrated_data), "\n")
cat("  Valid chlorophyll:", sum(!is.na(calibrated_data$chlorophyll_ug_L)), 
    sprintf("(%.1f%%)", 100 * mean(!is.na(calibrated_data$chlorophyll_ug_L))), "\n")
cat("  Valid positions (lat/lon):", sum(!is.na(calibrated_data$lat_dd) & !is.na(calibrated_data$lon_dd)),
    sprintf("(%.1f%%)", 100 * mean(!is.na(calibrated_data$lat_dd) & !is.na(calibrated_data$lon_dd))), "\n")
cat("  Valid depth:", sum(!is.na(calibrated_data$depth_best_m)),
    sprintf("(%.1f%%)", 100 * mean(!is.na(calibrated_data$depth_best_m))), "\n")
cat("  Valid temperature:", sum(!is.na(calibrated_data$sci_water_temp.degC)),
    sprintf("(%.1f%%)", 100 * mean(!is.na(calibrated_data$sci_water_temp.degC))), "\n")
cat("  Valid DO:", sum(!is.na(calibrated_data$DO_mg_L)),
    sprintf("(%.1f%%)", 100 * mean(!is.na(calibrated_data$DO_mg_L))), "\n")

# Chlorophyll statistics
cat("\nChlorophyll-a statistics (μg/L):\n")
print(summary(calibrated_data$chlorophyll_ug_L))

# Depth statistics
cat("\nDepth statistics (m):\n")
print(summary(calibrated_data$depth_best_m))

# Temperature statistics
cat("\nTemperature statistics (°C):\n")
print(summary(calibrated_data$sci_water_temp.degC))

# Geographic coverage
cat("\nGeographic coverage:\n")
cat("  Latitude range:", 
    sprintf("%.4f to %.4f", 
            min(calibrated_data$lat_dd, na.rm = TRUE),
            max(calibrated_data$lat_dd, na.rm = TRUE)), "\n")
cat("  Longitude range:", 
    sprintf("%.4f to %.4f", 
            min(calibrated_data$lon_dd, na.rm = TRUE),
            max(calibrated_data$lon_dd, na.rm = TRUE)), "\n")

# ============================================================
# 6. CREATE CLEAN DATASET FOR ANALYSIS
# ============================================================

# Filter for complete cases with valid chl, position, and depth
clean_data <- calibrated_data %>%
  filter(
    !is.na(chlorophyll_ug_L),
    !is.na(lat_dd),
    !is.na(lon_dd),
    !is.na(depth_best_m),
    chlorophyll_ug_L >= 0,  # Physical constraint
    depth_best_m >= 0       # Physical constraint
  )

cat("\nFinal clean dataset:\n")
cat("  Rows with complete data:", nrow(clean_data), 
    sprintf("(%.1f%% of total)", 100 * nrow(clean_data) / nrow(calibrated_data)), "\n")

# Check for potential outliers
cat("\nPotential outlier check:\n")
chl_q99 <- quantile(clean_data$chlorophyll_ug_L, 0.99, na.rm = TRUE)
chl_outliers <- sum(clean_data$chlorophyll_ug_L > chl_q99, na.rm = TRUE)
cat("  Chlorophyll > 99th percentile (", sprintf("%.2f", chl_q99), "μg/L):", 
    chl_outliers, "observations\n")

# Save diagnostic plot
cat("\nGenerating diagnostic plots...\n")

p_diag <- calibrated_data %>%
  select(datetime_utc, chlorophyll_ug_L, depth_best_m, sci_water_temp.degC, DO_mg_L) %>%
  pivot_longer(-datetime_utc, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = datetime_utc, y = value)) +
  geom_point(alpha = 0.3, size = 0.5) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  labs(title = "Time Series of Key Variables (After Calibration)",
       x = "Date/Time (UTC)", y = "Value") +
  theme_minimal()

print(p_diag)

cat("\n✓ Data preparation complete!\n")
cat("  Use 'calibrated_data' for full dataset\n")
cat("  Use 'clean_data' for analysis-ready data (complete cases only)\n")


# ============================================================
# 7. DIAGNOSTIC PLOTS
# ============================================================

cat("\nGenerating diagnostic plots...\n")

# Check if we have data to plot
vars_to_plot <- calibrated_data %>%
  select(datetime_utc, chlorophyll_ug_L, depth_best_m, sci_water_temp.degC, DO_mg_L) %>%
  pivot_longer(-datetime_utc, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))

if(nrow(vars_to_plot) > 0) {
  p_diag <- vars_to_plot %>%
    ggplot(aes(x = datetime_utc, y = value, color = variable)) +
    geom_point(alpha = 0.3, size = 0.5) +
    facet_wrap(~variable, scales = "free_y", ncol = 1,
               labeller = labeller(variable = c(
                 "chlorophyll_ug_L" = "Chlorophyll-a (μg/L)",
                 "depth_best_m" = "Depth (m)",
                 "sci_water_temp.degC" = "Temperature (°C)",
                 "DO_mg_L" = "Dissolved Oxygen (mg/L)"
               ))) +
    scale_color_manual(
      values = c(
        "chlorophyll_ug_L" = "#2ECC40",      # Green for chlorophyll
        "depth_best_m" = "#0074D9",          # Blue for depth
        "sci_water_temp.degC" = "#FF4136",   # Red for temperature
        "DO_mg_L" = "#FF851B"                # Orange for DO
      ),
      guide = "none"  # Hide legend since colors are obvious from facets
    ) +
    labs(title = "Time Series of Key Variables",
         x = "Date/Time (UTC)", y = "Value") +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey90", color = NA)
    )
  
  print(p_diag)
} else {
  cat("No valid data to plot\n")
}

cat("\n✓ Data preparation complete!\n")
cat("  Use 'calibrated_data' for full dataset\n")
cat("  Use 'clean_data' for analysis-ready data (complete cases only)\n")


# ============================================================
# 7. Offshore distance using polygon method
# ============================================================
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)  # Install with: install.packages("rnaturalearthhires", 
#              repos = "http://packages.ropensci.org", 
#              type = "source")
library(ggplot2)


cat("Step 1: Loading geographic data...\n")

# Get states/provinces for US and Canada
na_states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

# Create bounding box for Lake Ontario region (CORRECTED)
ontario_bbox <- st_bbox(c(xmin = -80, ymin = 43, xmax = -76, ymax = 44.5))
attr(ontario_bbox, "crs") <- st_crs(4326)

# Create polygon from bounding box
ontario_region <- st_as_sfc(ontario_bbox)

cat("Step 2: Extracting Lake Ontario region...\n")

# Clip to Lake Ontario region
ontario_land <- st_crop(na_states, ontario_region)

cat("Step 3: Converting glider data to spatial object...\n")

# Convert glider data to spatial
glider_sf <- clean_data %>%
  filter(!is.na(lat_dd), !is.na(lon_dd)) %>%
  st_as_sf(coords = c("lon_dd", "lat_dd"), crs = 4326)

cat("  Found", nrow(glider_sf), "valid observations with coordinates\n")

cat("Step 4: Transforming to UTM projection...\n")

# Transform to UTM 17N for accurate distance calculation (in meters)
glider_utm <- st_transform(glider_sf, crs = 32617)
ontario_land_utm <- st_transform(ontario_land, crs = 32617)

cat("Step 5: Extracting coastline...\n")

# Get the boundary (coastline) of the land
ontario_coast <- st_union(ontario_land_utm) %>%
  st_boundary()

cat("Step 6: Calculating distances to shore (this may take a moment)...\n")

# Calculate distance to shore for each point
distances <- st_distance(glider_utm, ontario_coast)
glider_utm$dist_to_shore_m <- as.numeric(distances)

# Convert to km for easier interpretation
glider_utm$dist_to_shore_km <- glider_utm$dist_to_shore_m / 1000

cat("Step 7: Adding distances back to original dataframe...\n")

# Add back to original dataframe
clean_data_with_distance <- clean_data %>%
  mutate(
    row_id = row_number()
  ) %>%
  left_join(
    glider_utm %>%
      st_drop_geometry() %>%
      mutate(row_id = row_number()) %>%
      select(row_id, dist_to_shore_m, dist_to_shore_km),
    by = "row_id"
  ) %>%
  select(-row_id)

# Summary
cat("\n✓ Distance calculation complete!\n\n")
cat("Distance to shore statistics (km):\n")
print(summary(clean_data_with_distance$dist_to_shore_km))

# Additional statistics
cat("\nDetailed statistics:\n")
cat("  Min distance:", sprintf("%.2f km", min(clean_data_with_distance$dist_to_shore_km, na.rm = TRUE)), "\n")
cat("  Max distance:", sprintf("%.2f km", max(clean_data_with_distance$dist_to_shore_km, na.rm = TRUE)), "\n")
cat("  Mean distance:", sprintf("%.2f km", mean(clean_data_with_distance$dist_to_shore_km, na.rm = TRUE)), "\n")
cat("  Median distance:", sprintf("%.2f km", median(clean_data_with_distance$dist_to_shore_km, na.rm = TRUE)), "\n")

# Quick visualization of the map
cat("\nGenerating map...\n")

p_map <- ggplot() +
  geom_sf(data = ontario_land, fill = "grey80", color = "grey50", linewidth = 0.5) +
  geom_point(data = clean_data_with_distance %>% filter(!is.na(dist_to_shore_km)), 
             aes(x = lon_dd, y = lat_dd, color = dist_to_shore_km),
             alpha = 0.6, size = 1.5) +
  scale_color_viridis_c(name = "Distance to\nShore (km)", option = "plasma") +
  labs(title = "Glider Track - Distance to Shore",
       subtitle = "Lake Ontario, April 2025",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  coord_sf(xlim = c(-80, -76), ylim = c(43, 44.5))

print(p_map)

cat("\n✓ Use 'clean_data_with_distance' for analysis with distance to shore\n")



# =============
# PLOTS
# ============

library(ggplot2)
library(viridis)
library(akima)
library(dplyr)


# ============================================================
# Prepare data
# ============================================================

plot_data <- clean_data_with_distance %>%
  filter(!is.na(dist_to_shore_km),
         !is.na(depth_best_m),
         !is.na(chlorophyll_ug_L),
         chlorophyll_ug_L >= 0,
         depth_best_m >= 0)

cat("Using", nrow(plot_data), "observations\n")

# ============================================================
# Interpolate to create smooth contours
# ============================================================

cat("Interpolating data for smooth contours...\n")

interp_data <- with(plot_data, 
                    interp(x = dist_to_shore_km, 
                           y = depth_best_m, 
                           z = chlorophyll_ug_L,
                           xo = seq(min(dist_to_shore_km), max(dist_to_shore_km), length = 100),
                           yo = seq(min(depth_best_m), max(depth_best_m), length = 100),
                           linear = TRUE,
                           extrap = FALSE))

# Convert to dataframe
interp_df <- expand.grid(
  dist_km = interp_data$x,
  depth_m = interp_data$y
) %>%
  mutate(chlorophyll = as.vector(interp_data$z))

# ============================================================
# PLOT 2: Clean contour plot (recommended)
# ============================================================

p2 <- ggplot() +
  # Filled contours (interpolated background)
  geom_raster(data = interp_df %>% filter(!is.na(chlorophyll)),
              aes(x = dist_km, y = depth_m, fill = chlorophyll),
              interpolate = TRUE) +
  # White contour lines
  geom_contour(data = interp_df %>% filter(!is.na(chlorophyll)),
               aes(x = dist_km, y = depth_m, z = chlorophyll),
               color = "white", alpha = 0.6, linewidth = 0.5,
               bins = 10) +
  # Color scale
  scale_fill_viridis_c(
    name = "Chlorophyll-a\n(μg/L)",
    option = "viridis",
    na.value = "grey90"
  ) +
  # Flip y-axis so depth increases downward
  scale_y_reverse() +
  # Labels
  labs(
    title = "Chlorophyll-a Distribution: Offshore Distance vs Depth",
    subtitle = "Lake Ontario, April 2025",
    x = "Distance from Shore (km)",
    y = "Depth (m)"
  ) +
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p2)

# ============================================================
# PLOT 3: With actual data points overlaid
# ============================================================

p3 <- ggplot() +
  # Filled contours
  geom_raster(data = interp_df %>% filter(!is.na(chlorophyll)),
              aes(x = dist_km, y = depth_m, fill = chlorophyll),
              interpolate = TRUE) +
  # Contour lines
  geom_contour(data = interp_df %>% filter(!is.na(chlorophyll)),
               aes(x = dist_km, y = depth_m, z = chlorophyll),
               color = "white", alpha = 0.5, linewidth = 0.4,
               bins = 10) +
  # Actual data points
  geom_point(data = plot_data,
             aes(x = dist_to_shore_km, y = depth_best_m),
             color = "white", alpha = 0.15, size = 0.5) +
  # Color scale
  scale_fill_viridis_c(
    name = "Chlorophyll-a\n(μg/L)",
    option = "viridis",
    na.value = "grey90"
  ) +
  scale_y_reverse() +
  labs(
    title = "Chlorophyll-a Distribution: Offshore Distance vs Depth",
    subtitle = "White dots show glider measurement locations",
    x = "Distance from Shore (km)",
    y = "Depth (m)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

print(p3)


# ============================================================
# PLOT 4: Filled contour (alternative style)
# ============================================================

p4 <- ggplot(interp_df %>% filter(!is.na(chlorophyll)), 
             aes(x = dist_km, y = depth_m, z = chlorophyll)) +
  geom_contour_filled(bins = 12) +
  scale_y_reverse() +
  labs(
    title = "Chlorophyll-a Distribution (Discrete Contours - Interpolated data)",
    x = "Distance from Shore (km)",
    y = "Depth (m)",
    fill = "Chlorophyll-a\n(μg/L)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

print(p4)


# ============================================================
# TEMPERATURE: Distance vs Depth
# ============================================================

cat("\n1. Processing Temperature data...\n")

# Prepare temperature data
temp_data <- clean_data_with_distance %>%
  filter(!is.na(dist_to_shore_km),
         !is.na(depth_best_m),
         !is.na(sci_water_temp.degC),
         depth_best_m >= 0)

cat("  Using", nrow(temp_data), "temperature observations\n")

# Interpolate temperature
temp_interp <- with(temp_data, 
                    interp(x = dist_to_shore_km, 
                           y = depth_best_m, 
                           z = sci_water_temp.degC,
                           xo = seq(min(dist_to_shore_km), max(dist_to_shore_km), length = 100),
                           yo = seq(min(depth_best_m), max(depth_best_m), length = 100),
                           linear = TRUE,
                           extrap = FALSE))

# Convert to dataframe
temp_interp_df <- expand.grid(
  dist_km = temp_interp$x,
  depth_m = temp_interp$y
) %>%
  mutate(temperature = as.vector(temp_interp$z))

# Plot temperature
p_temp <- ggplot(temp_interp_df %>% filter(!is.na(temperature)), 
                 aes(x = dist_km, y = depth_m, z = temperature)) +
  geom_contour_filled(bins = 12) +
  scale_y_reverse() +
  labs(
    title = "Water Temperature Distribution: Offshore Distance vs Depth",
    subtitle = "Lake Ontario, April 2025",
    x = "Distance from Shore (km)",
    y = "Depth (m)",
    fill = "Temperature\n(°C)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right"
  )

print(p_temp)

# ============================================================
# DISSOLVED OXYGEN: Distance vs Depth
# ============================================================

cat("\n2. Processing Dissolved Oxygen data...\n")

# Prepare DO data
do_data <- clean_data_with_distance %>%
  filter(!is.na(dist_to_shore_km),
         !is.na(depth_best_m),
         !is.na(DO_mg_L),
         depth_best_m >= 0,
         DO_mg_L > 0)

cat("  Using", nrow(do_data), "dissolved oxygen observations\n")

# Interpolate DO
do_interp <- with(do_data, 
                  interp(x = dist_to_shore_km, 
                         y = depth_best_m, 
                         z = DO_mg_L,
                         xo = seq(min(dist_to_shore_km), max(dist_to_shore_km), length = 100),
                         yo = seq(min(depth_best_m), max(depth_best_m), length = 100),
                         linear = TRUE,
                         extrap = FALSE))

# Convert to dataframe
do_interp_df <- expand.grid(
  dist_km = do_interp$x,
  depth_m = do_interp$y
) %>%
  mutate(dissolved_oxygen = as.vector(do_interp$z))

# Plot DO
p_do <- ggplot(do_interp_df %>% filter(!is.na(dissolved_oxygen)), 
               aes(x = dist_km, y = depth_m, z = dissolved_oxygen)) +
  geom_contour_filled(bins = 12) +
  scale_y_reverse() +
  labs(
    title = "Dissolved Oxygen Distribution: Offshore Distance vs Depth",
    subtitle = "Lake Ontario, April 2025",
    x = "Distance from Shore (km)",
    y = "Depth (m)",
    fill = "Dissolved\nOxygen\n(mg/L)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right"
  )

print(p_do)



# ============================================================
# Prepare data for Chl temporal plot
# ============================================================

temporal_chl_data <- clean_data_with_distance %>%
  filter(!is.na(day),
         !is.na(depth_best_m),
         !is.na(chlorophyll_ug_L),
         chlorophyll_ug_L >= 0,
         depth_best_m >= 0)

cat("Using", nrow(temporal_chl_data), "observations\n")
cat("Date range:", as.character(min(temporal_chl_data$day)), "to", 
    as.character(max(temporal_chl_data$day)), "\n")

# ============================================================
# PLOT 1: Line plot - Depth vs Date colored by Chlorophyll
# ============================================================

cat("\n1. Creating line plot (depth profile over time)...\n")

p_chl_time <- ggplot(temporal_chl_data, 
                 aes(x = day, y = depth_best_m, color = chlorophyll_ug_L)) +
  geom_point(alpha = 0.6, size = 1) +
    scale_color_viridis_c(
    name = "Chlorophyll-a\n(μg/L)",
    option = "viridis"
  ) +
  scale_y_reverse() +
  labs(
    title = "Chlorophyll-a Temporal Variation with Depth",
    subtitle = "Each vertical line represents one day's depth profile",
    x = "Date",
    y = "Depth (m)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_chl_time)

# ============================================================
# PLOT 2: Alternative line plot - showing depth zones over time
# ============================================================

cat("\n2. Creating line plot by depth zones...\n")

# Create depth zones
temporal_binned <- temporal_chl_data %>%
  mutate(
    depth_zone = cut(depth_best_m, 
                     breaks = c(0, 10, 20, 30, 50, 100, Inf),
                     labels = c("0-10m", "10-20m", "20-30m", "30-50m", "50-100m", ">100m"))
  ) %>%
  group_by(day, depth_zone) %>%
  summarise(
    mean_chl = mean(chlorophyll_ug_L, na.rm = TRUE),
    sd_chl = sd(chlorophyll_ug_L, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3)

p_line_zones <- ggplot(temporal_binned, aes(x = day, y = mean_chl, color = depth_zone)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_viridis_d(name = "Depth Zone", option = "plasma") +
  labs(
    title = "Chlorophyll-a Temporal Trends by Depth Zone",
    subtitle = "Mean chlorophyll concentration over time at different depths",
    x = "Date",
    y = "Mean Chlorophyll-a (μg/L)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_line_zones)

# ============================================================
# PLOT 3: Contour plot - Date vs Depth with Chlorophyll
# ============================================================

cat("\n3. Creating contour plot (time vs depth)...\n")

# Convert date to numeric for interpolation
temporal_chl_data_numeric <- temporal_chl_data %>%
  mutate(date_numeric = as.numeric(day))

# Interpolate
temporal_interp <- with(temporal_chl_data_numeric, 
                        interp(x = date_numeric, 
                               y = depth_best_m, 
                               z = chlorophyll_ug_L,
                               xo = seq(min(date_numeric), max(date_numeric), length = 100),
                               yo = seq(min(depth_best_m), max(depth_best_m), length = 100),
                               linear = TRUE,
                               extrap = FALSE))

# Convert back to dataframe
temporal_interp_df <- expand.grid(
  date_numeric = temporal_interp$x,
  depth_m = temporal_interp$y
) %>%
  mutate(
    chlorophyll = as.vector(temporal_interp$z),
    date = as.Date(date_numeric, origin = "1970-01-01")
  )

# Create contour plot
p_contour <- ggplot(temporal_interp_df %>% filter(!is.na(chlorophyll)), 
                    aes(x = date, y = depth_m, z = chlorophyll)) +
  geom_contour_filled(bins = 12) +
  scale_y_reverse() +
  labs(
    title = "Chlorophyll-a Temporal Distribution: Date vs Depth",
    subtitle = "Lake Ontario, April 2025",
    x = "Date",
    y = "Depth (m)",
    fill = "Chlorophyll-a\n(μg/L)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_contour)

# ============================================================
# PLOT 4: Smooth raster contour (alternative style)
# ============================================================

cat("\n4. Creating smooth raster contour plot...\n")

p_contour_smooth <- ggplot() +
  geom_raster(data = temporal_interp_df %>% filter(!is.na(chlorophyll)),
              aes(x = date, y = depth_m, fill = chlorophyll),
              interpolate = TRUE) +
  geom_contour(data = temporal_interp_df %>% filter(!is.na(chlorophyll)),
               aes(x = date, y = depth_m, z = chlorophyll),
               color = "white", alpha = 0.4, linewidth = 0.3,
               bins = 10) +
  scale_fill_viridis_c(
    name = "Chlorophyll-a\n(μg/L)",
    option = "viridis",
    na.value = "grey90"
  ) +
  scale_y_reverse() +
  labs(
    title = "Chlorophyll-a Temporal Distribution (Smooth)",
    subtitle = "Date vs Depth with white contour lines",
    x = "Date",
    y = "Depth (m)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_contour_smooth)

# ============================================================
# PLOT 5: Faceted by depth - showing temporal trends
# ============================================================

cat("\n5. Creating faceted depth profile plots...\n")

# Select specific depths to highlight
depth_layers <- temporal_chl_data %>%
  mutate(
    depth_layer = case_when(
      depth_best_m <= 10 ~ "Surface (0-10m)",
      depth_best_m > 10 & depth_best_m <= 20 ~ "Epilimnion (10-20m)",
      depth_best_m > 20 & depth_best_m <= 40 ~ "Metalimnion (20-40m)",
      depth_best_m > 40 ~ "Hypolimnion (>40m)"
    )
  ) %>%
  filter(!is.na(depth_layer)) %>%
  group_by(day, depth_layer) %>%
  summarise(
    mean_chl = mean(chlorophyll_ug_L, na.rm = TRUE),
    sd_chl = sd(chlorophyll_ug_L, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3)

p_faceted <- ggplot(depth_layers, aes(x = day, y = mean_chl)) +
  geom_ribbon(aes(ymin = mean_chl - sd_chl, ymax = mean_chl + sd_chl),
              alpha = 0.2, fill = "darkgreen") +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 2) +
  facet_wrap(~depth_layer, ncol = 1, scales = "free_y") +
  labs(
    title = "Chlorophyll-a Temporal Trends by Water Column Layer",
    subtitle = "Mean ± SD for each depth layer over time",
    x = "Date",
    y = "Chlorophyll-a (μg/L)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_faceted)



