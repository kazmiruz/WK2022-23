library(tidyverse)

library(geosphere)
library(sf)

#####################
### CTD compiling ###

#
## putting together metadata for our microbial samples based on all the info available in the reports and adding genohub_id

meta.wk.genohub <- read.csv("WK_DNA/meta.wk.genohub.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  select(-Date, -Code, -comments) %>%
  mutate(across(Cruise, ~ gsub("_", "", as.character(.))),                    # remove underscore from Cruise
         across(genohub_id, ~ gsub("_", "-", as.character(.))),               # adjust genohub_id
         plate.geno = toupper(paste0(plate, "_", genohub_id)),                # Make plate.geno
         join_var = ifelse(Cruise == "WK2022", `Depth [m]`, Target)) %>%      # assign Depth for WK2022, Target for WK2023
  select(-`Depth [m]`, -Target)

meta.wk.report <- read.csv("WK_DNA/meta.wk.report.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(join_var = ifelse(Cruise == "WK2022", `Depth [m]`, Target))          # assign Depth for WK2022, Target for WK2023

meta.wk <- meta.wk.genohub %>%
  inner_join(meta.wk.report, by = c("Cruise", "Station", "join_var")) %>%     # 1 sample whose origin I couldn't figure out is out
  separate(Date, into = c("Month", "Day", "Year"), sep = "/", convert = TRUE) %>%
  mutate (`Depth [m]` = as.numeric(`Depth [m]`)) %>%
  select(-join_var)


#
## wrangling CTD and TSG data from CanWin

ctd.wk <- rbind(
  
  # Autonomous 2022
  read.csv("WK_CTD/excel_2022_autonomous_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
    select(-Notes) %>%
    rename_with(~"Depth [m]", contains("Depth [salt water, m]")) %>%
    mutate (Cast = NA,
            across(contains("SBE 43"), ~ NA)),                 # in 2022, the oxygen data from the autonomous CTD was corrupted
  
  # Rosette 2022
  # Before using the csv, delete unnecessary some columns:
  # have to remove readings of the second `Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]` instrument;
  # since `PAR/Irradiance, Biospherical/Licor [umol photons/m^2/sec]`is enough,
  # have to remove readings of `SPAR, Biospherical/Licor [umol photons/m^2/sec]` and `CPAR/Corrected Irradiance [%]`
  read.csv("WK_CTD/excel_2022_rosette_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
    select(-Notes) %>%
    rename_with(~"Depth [m]", contains("Depth [salt water, m]")) %>%
    mutate(`Latitude Out [Degrees North]` = NA,
           `Longitude Out [Degrees East]` = NA,
           `Beam Transmission, WET Labs C-Star [%]` = NA,
           `Beam Attenuation, WET Labs C-Star [1/m]` = NA) %>%
    dplyr::rename("Latitude In [Degrees North]" = "Latitude [Degrees North]",
                  "Longitude In [Degrees East]" = "Longitude [Degrees East]"),
  
  # Autonomous 2023  
  read.csv("WK_CTD/2023_autonomous_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
    select(-`Logbook Notes`) %>%
    rename_with(~"Depth [m]", contains("Depth [salt water, m]")) %>%
    mutate (Cast = NA),
  
  # Rosette 2023
  # Before using the csv, delete unnecessary some columns:
  # since `PAR/Irradiance, Biospherical/Licor [umol photons/m^2/sec]`is enough,
  # have to remove readings of `SPAR, Biospherical/Licor [umol photons/m^2/sec]` and `CPAR/Corrected Irradiance [%]`
  read.csv("WK_CTD/excel_2023_rosette_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
    select(-Notes) %>%
    rename_with(~"Depth [m]", contains("Depth [salt water, m]")) %>%
    mutate(`Fluorescence, WET Labs CDOM [mg/m^3]` = NA,
           `Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]` = NA,
                         `Latitude Out [Degrees North]` = NA,
                         `Longitude Out [Degrees East]` = NA,
                         `Beam Transmission, WET Labs C-Star [%]` = NA,
                         `Beam Attenuation, WET Labs C-Star [1/m]` = NA) %>%
    dplyr::rename("Latitude In [Degrees North]" = "Latitude [Degrees North]",
                  "Longitude In [Degrees East]" = "Longitude [Degrees East]")) %>%
  
  mutate(`Station` = as.character(Station),                                                 # force Station to character
         `Latitude [deg]` = coalesce(`Latitude In [Degrees North]`,
                                     `Latitude Out [Degrees North]`),                       # Lat in, fall back on Lat out
         `Longitude [deg]` = coalesce(`Longitude In [Degrees East]`,
                                      `Longitude Out [Degrees East]`),                      # Lon in, fall back on Lon out
         Cruise = ifelse(grepl("2022", Cruise), "WK2022",
                         ifelse(grepl("2023", Cruise), "WK2023", Cruise))) %>%              # shorten Cruise name
  select(-contains(" In"), -contains(" Out"),                                               # remove Latitude/Longitude In/Out
         -contains("Pressure"),                                                             # remove Pressure b/c we have Depth
         -contains("SBE 43 [ml/l]"), -contains("SBE 43 [% saturation]"),                    # remove extra Oxygen values
         -contains("Conductivity"), -contains("Volume Anomaly"), -`Timestamp (UTC)`) %>%    # remove unnecessary/transient variables
  rename_with(~ "CDOM [mg/m^3]", contains("CDOM")) %>%
  rename_with(~ "Sigma-theta [kg/m^3]", contains("sigma-theta")) %>%
  rename_with(~ gsub(" \\(UTC\\)", "", .x)) %>%
  rename_with(~ gsub("ITS-90, ", "", .x)) %>%                          # Remove "ITS-90, "
  rename_with(~ gsub("/Irradiance", "", .x)) %>%                       # Remove "/Irradiance"
  rename_with(~ gsub(", [^\\[]*", " ", .x)) %>%                        # Remove everything after comma but before bracket `[`
  rename_with(~ trimws(.x)) %>%                                        # Remove leading and trailing spaces from column names
  filter(!is.na(Year) & !is.na(Month) & !is.na(Day) & 
           !is.na(`Depth [m]`) & !is.na(`Latitude [deg]`) & !is.na(`Longitude [deg]`)) # ensure there are no NA values in key columns


library(readxl)
# --- WK2022 ---
tsg_wk2022 <- read.csv("WK_CTD/tsg_2022_merged.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(Cruise = "WK2022",
         Timestamp = ymd_hm(Timestamp))  # Parse WK2022 format

# --- WK2023 ---
tsg_wk2023 <- list.files("WK_CTD/WK TSG data 2023", pattern = "\\.xlsx$", full.names = TRUE) %>%
  map_dfr(read_xlsx) %>%
  select(-timeQ, -`c0S/m`) %>%
  mutate(Timestamp = ymd_hms(date_time)) %>%  # Parse WK2023 format
  select(-date_time) %>%
  dplyr::rename(`Latitude [deg]` = latitude,
                `Longitude [deg]` = longitude,
                `Salinity [PSU]` = sal00,
                `Temperature [deg C]` = t090C) %>%
  mutate(Cruise = "WK2023")

# --- Combine and standardize ---
tsg.wk <- bind_rows(tsg_wk2022, tsg_wk2023) %>%
  mutate(`Depth [m]` = 2,
         Station = "TSG",
         Cast = NA,
         Year = year(Timestamp),
         Month = month(Timestamp),
         Day = day(Timestamp),
         Hour = hour(Timestamp),
         Minute = minute(Timestamp),
         Second = 0) %>%
  select(-Timestamp)


ctd.tsg.wk <- bind_rows(ctd.wk, tsg.wk) %>%
  distinct(`Latitude [deg]`, `Longitude [deg]`, `Depth [m]`, Year, Month, Day, .keep_all = TRUE)


#
## Add CTD data to the metadata

meta_ctd_candidates <- dplyr::left_join(
  meta.wk,
  ctd.tsg.wk %>%
    dplyr::rename(Station_CTD = Station) %>%
    dplyr::select(-Cruise, -Hour, -Minute, -Second),
  by = c("Year", "Month", "Day")
)


meta.wk.incomplete <- meta_ctd_candidates %>%
  
  # Calculate lat/lon/depth diffs
  dplyr::mutate(
    lat_diff = abs(`Latitude [deg].y` - `Latitude [deg].x`),
    lon_diff = abs(`Longitude [deg].y` - `Longitude [deg].x`),
    depth_diff = abs(`Depth [m].x` - `Depth [m].y`),
    geo_diff = lat_diff + lon_diff,
    
    # Count non-NA CTD variables (y columns)
    non_na_count = rowSums(!is.na(select(., ends_with(".y"))))
  ) %>%
  
  # Apply xyz dimension filters (0.020°Lat ≈ 2.2 km; 0.015°Lon ≈ 0.9 km; wait for 2023 TSG)
  dplyr::filter(lat_diff <= 0.020 & lon_diff <= 0.015 & depth_diff <= 1) %>%  
  
  # Grouping
  dplyr::group_by(Year, Month, Day, Station, `Depth [m].x`) %>%
  
  # Since there are multiple legitimate matches, select the ones w/ the most complete cases and then w/ the least Lat/Lon dif-ce
  dplyr::arrange(desc(non_na_count), geo_diff, depth_diff) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  
  # Clean columns
  dplyr::rename_with(~ sub("\\.x$", "", .x)) %>%
  dplyr::select(-ends_with(".y"),
                -lat_diff, -lon_diff, -depth_diff, -geo_diff, -non_na_count,
                -plate, -Sample, -Note) %>%
  
  # Organize final columns
  dplyr::select ("Cruise", "Year", "Month", "Day", "Time",
                 "Station", "Station_CTD", "Cast", "Latitude [deg]", "Longitude [deg]", "Bot. Depth [m]",
                 "Depth [m]", "Code", "Target",
                 "genohub_id", "plate.geno", everything())

# hist((meta.wk.complete %>% filter(lat_diff < 0.05))$lat_diff)
# hist((meta.wk.complete %>% filter(lon_diff < 0.05))$lon_diff)

meta.wk.complete <- bind_rows(
  meta.wk.incomplete,
  meta.wk %>% 
    anti_join(meta.wk.incomplete, by = c("Year", "Month", "Day", "Station", "Depth [m]")) %>%
    tibble::add_column(!!!setNames(rep(list(NA), 
                                       length(setdiff(colnames(meta.wk.incomplete), colnames(meta.wk)))), 
                                   setdiff(colnames(meta.wk.incomplete), colnames(meta.wk))))) %>%
  select(-plate, -Sample, -Note) %>%
  
  # Distance to coast
  st_as_sf(coords = c("Longitude [deg]", "Latitude [deg]"), crs = 4326, remove = FALSE) %>%
  mutate("Distance to coast [km]" = st_distance(., st_read("ne_110m_coastline/ne_110m_coastline.shp")) %>%
           apply(1, min) %>% as.numeric() / 1000) %>%
  st_drop_geometry()


#
## save to csv
write.csv(meta.wk.complete, "OUTPUT/meta.wk.complete.csv", row.names = FALSE, fileEncoding = "UTF-8")

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################


##########################
### EXPLORE and OUTPUT ###

# Open a PDF to save the histograms
pdf("OUTPUT/wk_parameters_hist.pdf", width = 8.5, height = 12)

# Set up 4 plots per page
par(mfrow = c(3, 2), mar = c(5, 6, 4, 3))  # 3 rows and 2 columns per page; wider left margin

# Loop through numeric columns
for (col_name in names(meta.wk.complete)) {
  # Check if the column is numeric
  if (is.numeric(meta.wk.complete[[col_name]])) {
    # Extract the column, removing NA values
    column_data <- na.omit(meta.wk.complete[[col_name]])
    
    # Skip if the column is empty after removing NAs
    if (length(column_data) == 0) next
    
    # Calculate basic statistics
    col_min <- min(column_data)
    col_max <- max(column_data)
    col_mean <- mean(column_data)
    col_sd <- sd(column_data)
    
    # Generate x-axis range based on the data
    x_range <- range(column_data, na.rm = TRUE)
    
    # Create a simple histogram
    hist(column_data,
         main = NULL,
         xlab = bquote(bold(.(col_name))),
         xlim = x_range,  # Use calculated range for x-axis
         col = "lightblue",
         border = "black")
    
    # Add summary statistics above the graph
    mtext(paste("Min:", round(col_min, 2)), side = 3, line = 2, adj = 0, cex = 0.8)
    mtext(paste("Max:", round(col_max, 2)), side = 3, line = 1, adj = 0, cex = 0.8)
    mtext(paste("Mean ± SD:", round(col_mean, 2), "±", round(col_sd, 2)), side = 3, line = 0, adj = 0, cex = 0.8)
  }
}

# Close the PDF
dev.off()
