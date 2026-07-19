# WK_1---------------------------------------------------------------------
# Purpose:
#   Compile all environmental data (CTD, CastAway, TSG) and merge with sample metadata to produce a complete metadata table for
#   the WK microbial samples. Derived parameters (light attenuation, mixed‑layer depth, distance to coast) are also calculated.
#
# Required external files:
#   WK_CTD/excel_2022_autonomous_ctd_downcasts.csv
#   WK_CTD/2023_autonomous_ctd_downcasts.csv
#   WK_CTD/excel_2022_rosette_ctd_downcasts.csv
#   WK_CTD/excel_2023_rosette_ctd_downcasts.csv
#   WK_CTD/2022_rosette_bottle_data.csv
#   WK_CTD/2023_rosette_bottle_data.csv
#   WK_CTD/tsg_2022_merged.csv
#   WK_CTD/WK_TSG_data_2023/  (folder containing *.xlsx files)
#   WK_CTD/JBE_2022_CastAway_CC1744004_Data/ (folder containing *.csv files)
#   WK_CTD/JBE_2023_CastAway_CC1744004_Data/ (folder containing *.csv files)
#   WK_DNA/meta.wk.genohub.csv
#   WK_DNA/meta.wk.report.approx.csv
#   ne_10m/ne_10m_coastline.shp (Natural Earth shapefile; other NE layers optional)
#
# Output:
#   OUTPUT/meta.wk.complete.csv     (final metadata table)
#   OUTPUT/wk_parameters_hist.pdf   (histograms of all numeric variables)
# -------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(sf)

##########################
### Compiling CTD data ###
# CTD_AUTO + CTD_DOWN + CTD_BOT + TSG + CastAway 

#
## Autonomous CTD
ctd.auto <- bind_rows(read.csv("WK_CTD/excel_2022_autonomous_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                        rename_with(~"Depth [m]", contains("Depth [salt water, m]")),
                      read.csv("WK_CTD/2023_autonomous_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                        rename_with(~"Depth [m]", contains("Depth [salt water, m]"))) %>%
  
  # Fix the variables
  mutate(Code = "ROS-A", Code_ODV = 1,
         Station = as.character(Station),                                                   # force Station to character
         Cruise = ifelse(grepl("2022", Cruise), "WK2022",
                         ifelse(grepl("2023", Cruise), "WK2023", Cruise)),                  # shorten Cruise name
         `Latitude [deg]` = coalesce(`Latitude In [Degrees North]`,
                                     `Latitude Out [Degrees North]`),                       # Lat in, fall back on Lat out
         `Longitude [deg]` = coalesce(`Longitude In [Degrees East]`,
                                      `Longitude Out [Degrees East]`)) %>%                  # Lon in, fall back on Lon out
  
  select(-contains("Notes"),
         -`Beam Transmission, WET Labs C-Star [%]`,                                         # remove Transmission, leave only Attenuation (the 1/m one)
         -contains(" In"), -contains(" Out"),                                               # remove Latitude/Longitude In/Out
         -contains("Pressure"),                                                             # remove Pressure b/c we have Depth
         -contains("SBE 43 [ml/l]"), -contains("SBE 43 [% saturation]"),                    # remove extra Oxygen values
         -contains("Conductivity"), -contains("Volume Anomaly"),                            # remove unnecessary/transient variables
         -contains("Time")) %>%                                                             # remove UTC Time
  
  rename_with(~ "Chlorophyll-a fluorescence [mg/m^3]", contains("ECO-AFL/FL")) %>%
  rename_with(~ "CDOM fluorescence [mg/m^3]", contains("CDOM")) %>%
  rename_with(~ "Beam attenuation [1/m]", contains ("Beam attenuation")) %>%
  rename_with(~ "Sigma-theta [kg/m^3]", contains("sigma-theta")) %>%
  rename_with(~ gsub(" \\(UTC\\)", "", .x)) %>%
  rename_with(~ gsub("ITS-90, ", "", .x)) %>%                          # Remove "ITS-90, "
  rename_with(~ gsub("/Irradiance", "", .x)) %>%                       # Remove "/Irradiance"
  rename_with(~ gsub(", [^\\[]*", " ", .x)) %>%                        # Remove everything after comma but before bracket `[`
  rename_with(~ trimws(.x)) %>%                                        # Remove leading and trailing spaces from column names
  filter(!is.na(Year) & !is.na(Month) & !is.na(Day) & 
           !is.na(`Depth [m]`) & !is.na(`Latitude [deg]`) & !is.na(`Longitude [deg]`)) %>%  # ensure there are no NA values in key columns
  mutate(`Oxygen [umol/kg]` = if_else(Cruise == "WK2022", NA_real_, `Oxygen [umol/kg]`))    # in 2022, the O2 data from the auto CTD is wrong (the CTD was dropped)


#
## Rosette downcast CTD
ctd.ros_down <- bind_rows(read.csv("WK_CTD/excel_2022_rosette_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE),
            read.csv("WK_CTD/excel_2023_rosette_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE)) %>%
  
  # Fix the variables
  mutate(Code = "ROS-D", Code_ODV = 2,
         Station = as.character(Station),                                               # force Station to character
         Cruise = ifelse(grepl("2022", Cruise), "WK2022",
                         ifelse(grepl("2023", Cruise), "WK2023", Cruise))) %>%          # shorten Cruise name

  select(-Notes, - `Fluorescence, WET Labs ECO-AFL/FL, 2 [mg/m^3]`,                     # 2nd instrument - there are negative values and a lot are 0s
         -contains("Pressure"),                                                         # remove Pressure b/c we have Depth
         -contains("SBE 43 [ml/l]"), -contains("SBE 43 [% saturation]"),                # remove extra Oxygen values
         -contains("Conductivity"),                                                     # remove transient variable Conductivity
         -contains("Volume Anomaly"),                                                   # remove Specific Volume Anomaly
         -contains("SPAR"), -contains ("CPAR"),                                         # remove extra PAR
         -contains("Time")) %>%                                                         # remove UTC Time
  
  rename_with(~ "Latitude [deg]", contains("Lat")) %>%
  rename_with(~ "Longitude [deg]", contains("Lon")) %>%
  rename_with(~"Depth [m]", contains("Depth [salt water, m]")) %>%
  rename_with(~ "Chlorophyll-a fluorescence [mg/m^3]", contains("ECO-AFL/FL")) %>%
  rename_with(~ "CDOM fluorescence [mg/m^3]", contains("CDOM")) %>%
  rename_with(~ "Sigma-theta [kg/m^3]", contains("sigma-theta")) %>%
  rename_with(~ gsub(" \\(UTC\\)", "", .x)) %>%
  rename_with(~ gsub("ITS-90, ", "", .x)) %>%                          # Remove "ITS-90, "
  rename_with(~ gsub("/Irradiance", "", .x)) %>%                       # Remove "/Irradiance"
  rename_with(~ gsub(", [^\\[]*", " ", .x)) %>%                        # Remove everything after comma but before bracket `[`
  rename_with(~ trimws(.x)) %>%                                        # Remove leading and trailing spaces from column names

  filter(!is.na(Year) & !is.na(Month) & !is.na(Day) & 
           !is.na(`Depth [m]`) & !is.na(`Latitude [deg]`) & !is.na(`Longitude [deg]`)) # ensure there are no NA values in key columns


#
## Upcast CTD for Rosette bottle
ctd.ros_bot <- bind_rows(read.csv("WK_CTD/2022_rosette_bottle_data.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                           mutate(Cruise = "WK2022"),
                          read.csv("WK_CTD/2023_rosette_bottle_data.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                           mutate(Cruise = "WK2023")) %>%
  
  # Fix the variables
  mutate(Code = "ROS-B", Code_ODV = 3,
         Station = as.character(Station)) %>%                                           # force Station to character
  separate(`Time (UTC)`, into = c("Hour", "Minute", "Second"), sep = ":", convert = TRUE) %>%
  
  select(-`Fluorescence, WET Labs ECO-AFL/FL, 2 [mg/m^3]`,                              # 2nd instrument - there are negative values and a lot are 0s
         -contains("sdev"),                                                             # remove sdev, min, max
         -contains("Pressure"),                                                         # remove Pressure b/c we have Depth
         -contains("SBE 43 [ml/l]"), -contains("SBE 43 [% saturation]"),                # remove extra Oxygen values
         -contains("Conductivity"),                                                     # remove transient variable Conductivity
         -contains("Volume Anomaly"),                                                   # remove Specific Volume Anomaly
         -contains("SPAR"), -contains ("CPAR"),                                         # remove extra PAR
         -contains("Time"),                                                             # remove UTC Time
         -contains("Descent")) %>%                                                      # remove Descent Rate
  
  rename_with(~ "Latitude [deg]", contains("Lat")) %>%
  rename_with(~ "Longitude [deg]", contains("Lon")) %>%
  rename_with(~ "Depth [m]", contains("Depth [salt water, m]")) %>%
  rename_with(~ "Chlorophyll-a fluorescence [mg/m^3]", contains("ECO-AFL/FL")) %>%
  rename_with(~ "CDOM fluorescence [mg/m^3]", contains("CDOM")) %>%
  rename_with(~ "Sigma-theta [kg/m^3]", contains("sigma-theta")) %>%   # rename this "density" to Sigma-theta
  rename_with(~ "Density [kg/m^3]", contains("density")) %>%           # rename this "density" to Density
  rename_with(~ gsub(" \\(UTC\\)", "", .x)) %>%
  rename_with(~ gsub("ITS-90, ", "", .x)) %>%                          # Remove "ITS-90, "
  rename_with(~ gsub("/Irradiance", "", .x)) %>%                       # Remove "/Irradiance"
  rename_with(~ gsub(", [^\\[]*", " ", .x)) %>%                        # Remove everything after comma but before bracket `[`
  rename_with(~ trimws(.x)) %>%                                        # Remove leading and trailing spaces from column names
  
  filter(!is.na(Year) & !is.na(Month) & !is.na(Day) & 
           !is.na(`Depth [m]`) & !is.na(`Latitude [deg]`) & !is.na(`Longitude [deg]`))  # ensure there are no NA values in key columns


#
## CastAway for Niskin (acknowledge Zou Zou Kuzyk)

# Folders
folders <- list(WK2022 = "WK_CTD/JBE_2022_CastAway_CC1744004_Data",
                WK2023 = "WK_CTD/JBE_2023_CastAway_CC1744004_Data")

# Function
process_castaway <- function(file, cruise) {
  lines <- readLines(file, warn = FALSE)
  header_idx <- which(!grepl("^%", lines))[1]
  metadata <- lines[1:(header_idx - 1)]
  data_text <- paste(lines[header_idx:length(lines)], collapse = "\n")
  
  # Extract metadata
  get_val <- function(pat) strsplit(grep(pat, metadata, value = TRUE), ",")[[1]][2]
  cast_time <- ymd_hms(get_val("^% Cast time \\(UTC\\)"))
  lat <- as.numeric(get_val("^% Start latitude"))
  lon <- as.numeric(get_val("^% Start longitude"))
  
  # Parse data
  df <- read.csv(text = data_text, header = TRUE, check.names = FALSE) %>%
    select(`Depth [m]` = `Depth (Meter)`,
           `Temperature [deg C]` = `Temperature (Celsius)`,
           `Salinity [PSU]` = `Salinity (Practical Salinity Scale)`,
           `Density [kg/m^3]` = `Density (Kilograms per Cubic Meter)`) %>%
    mutate(`Latitude [deg]` = lat,
           `Longitude [deg]` = lon,
           Cruise = cruise,
           Station = "CastAway",
           Code = "NISK", Code_ODV = 4,
           Year = year(cast_time), Month = month(cast_time), Day = day(cast_time),
           Hour = hour(cast_time), Minute = minute(cast_time), Second = as.numeric(format(cast_time, "%S")))
}

# Applying the function and binding (some CastAway don't have Lat/Long but this is OK b/c we know where they were based on time)
castaway <- bind_rows(lapply(list.files(folders$WK2022, "\\.csv$", full.names = TRUE), process_castaway, cruise = "WK2022"),
                      lapply(list.files(folders$WK2023, "\\.csv$", full.names = TRUE), process_castaway, cruise = "WK2023"))


#
## TSG for FT
tsg <- bind_rows(read.csv("WK_CTD/tsg_2022_merged.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                   mutate(Cruise = "WK2022", date_time = ymd_hms(date_time)),
                 list.files("WK_CTD/WK_TSG_data_2023", pattern = "\\.xlsx$", full.names = TRUE) %>% map_dfr(read_xlsx) %>% mutate(Cruise = "WK2023")) %>%
  rename(`Latitude [deg]` = latitude,
         `Longitude [deg]` = longitude,
         `Salinity [PSU]` = sal00,
         `Temperature [deg C]` = t090C) %>%
  mutate(`Depth [m]` = 2,
         Station = "TSG",
         Code = "FT", Code_ODV = 5,
         Year = year(date_time),
         Month = month(date_time),
         Day = day(date_time),
         Hour = hour(date_time),
         Minute = minute(date_time),
         Second = as.numeric(format(date_time, "%S"))) %>% select(-timeQ, -date_time, -`c0S/m`) %>%
  filter(`Salinity [PSU]` < 35)                                                         # some Salinity is off


#
## Put together all the oceanographic data and derive light attenuation and MLD

ctd <- bind_rows(ctd.auto, ctd.ros_down, ctd.ros_bot, castaway, tsg) %>% filter(!is.na(`Salinity [PSU]`)) %>% # join instrumental readings and rmv rows w/out S
  group_by(Year, Cruise, Month, Day, Station, Code, Cast) %>%     # group per PROFILE (date + station + cast) so Kd/MLD reference the correct cast, never crossing dates/casts

  # Calculate light attenuation coefficient: Kirk (2011) and Weiskerger et al. (2018).
  mutate(
    .par_valid = is.finite(`Depth [m]`) & is.finite(`PAR [umol photons/m^2/sec]`) & `PAR [umol photons/m^2/sec]` > 0, # rows usable as a surface reference
    .E0 = if (any(.par_valid)) `PAR [umol photons/m^2/sec]`[.par_valid][which.min(`Depth [m]`[.par_valid])] else NA_real_, # E0 = PAR at shallowest VALID (finite, >0) depth of this profile; NA if none
    `Light attenuation coefficient [1/m]` = if_else(
      `Depth [m]` > 0 & `PAR [umol photons/m^2/sec]` > 0 & !is.na(.E0),  # Apply formula only if D > 0, valid PAR, and a valid surface reference exists
      (1 / `Depth [m]`) *                                            # Kd = 1/z*ln(E0/Ez)
        log(.E0 / pmax(`PAR [umol photons/m^2/sec]`, 1e-6)),         # ensure there is no division by 0
      NA_real_) |> (\(x) if_else(x < 0, 0, x))()) %>%                # Remove negative Kd values (when Ez>E0 aka artifact)
  select(-.par_valid, -.E0) %>%

  # Calculate MLD only for CTD stations (ROS-A, ROS-D, ROS-B) and where sigma-theta is available: Peralta-Ferriz and Woodgate (2015)
  mutate(`Mixed Layer Depth [m]` = if (Code[1] %in% c("ROS-A", "ROS-D")) {   # do it only for autonomous and downcast CTD data (upcasts don't have all the depth)
    sigma_ref <- `Sigma-theta [kg/m^3]`[which.min(`Depth [m]`)]                  # σθ at the shallowest measured depth aka reference depth
    mixed_depths <- `Depth [m]`[`Sigma-theta [kg/m^3]` <= (sigma_ref + 0.1)]     # find all depths where σθ ≤ (surface σθ + 0.1 kg/m³) aka w/in threshold
    if (length(mixed_depths) > 0) max(mixed_depths) else max(`Depth [m]`)} else {NA_real_}) %>% ungroup() # get deepest MLD; none exceed threshold - MLD = max depth

# with(ctd %>% filter(`Depth [m]` >= 1 & `Depth [m]` <= 40) %>% drop_na(`Beam attenuation [1/m]`, `Light attenuation coefficient [1/m]`), {
#   plot(`Beam attenuation [1/m]`, `Light attenuation coefficient [1/m]`)
#   lines(lowess(`Beam attenuation [1/m]`, `Light attenuation coefficient [1/m]`), col = "red", lwd = 2)
# })

#
## putting together metadata for our microbial samples based on all the info available in the reports and adding genohub_id
meta.wk.genohub <- read.csv("WK_DNA/meta.wk.genohub.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(across(Cruise, ~ gsub("_", "", as.character(.))),                       # remove underscore from Cruise
         across(genohub_id, ~ gsub("_", "-", as.character(.))),                  # adjust genohub_id
         plate.geno = toupper(paste0(plate, "_", genohub_id)),                   # Make plate.geno
         join_var = ifelse(Cruise == "WK2022", `Depth [m]`, Target)) %>%         # join variable: Depth for WK2022, Target for WK2023
  select(-`Depth [m]`, -Target)

meta.wk.report <- read.csv("WK_DNA/meta.wk.report.approx.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(join_var = ifelse(Cruise == "WK2022", `Depth [m]`, Target))             # join variable: Depth for WK2022, Target for WK2023

meta.wk <- meta.wk.genohub %>%
  inner_join(meta.wk.report, by = c("Cruise", "Station", "join_var")) %>%
  separate(Date, into = c("Month", "Day", "Year"), sep = "/", convert = TRUE) %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%    # use the HMS time from your own meta.wk record
  mutate (Second  = 00) %>%
  select(-join_var) %>%
  mutate(`Depth [m]` = as.numeric(`Depth [m]`))

#############

# #
# ## save to csv for ODV to assess proximity of your sampled stations & the stations you use to approximate the metadata
# # Code = 0 for your sampled stations & Code 1-4 for CTD and CastAway; ignore TSG here
# write.csv(bind_rows(ctd %>% filter(Code_ODV < 5), meta.wk %>% mutate(Code_ODV = 0)), "OUTPUT/ctd.odv.csv", row.names = FALSE, fileEncoding = "UTF-8")

#############

#
## Add metadata for samples

### add BOTTLE (upcast) CTD ###

# build “ROS-A” (AUTO CTD) lookup table - will use later
ctd_ref <- ctd %>%
  filter(Code == "ROS-A") %>%
  transmute(Cruise, Year, Month,
            Day_ref = Day, Station_ref = Station, Code_ref = Code,
            Lon_ref = `Longitude [deg]`, Lat_ref = `Latitude [deg]`, Dep_ref = `Depth [m]`,
            CDOM_ref = `CDOM fluorescence [mg/m^3]`, Fluo_ref = `Chlorophyll-a fluorescence [mg/m^3]`, BA_ref = `Beam attenuation [1/m]`)

# --- Step 1: get original ROS-B data (before any approximation) ---
rosb_original <- meta.wk %>% filter(`Approx Code` == "ROS-B") %>%
  left_join(ctd %>% filter(Code == "ROS-B") %>%
              select(-`Bot. Depth [m]`, -`Mixed Layer Depth [m]`, -Hour, -Minute, -Second),
            by = join_by(Code, Cruise, Station, Cast, `Bottle Position`, Year, Month, Day),
            suffix = c(".meta", ".ctd")) %>%
  mutate(`Latitude [deg]`  = `Latitude [deg].ctd`,
         `Longitude [deg]` = `Longitude [deg].ctd`,
         `Depth [m]` = `Depth [m].ctd`) %>%
  select(-ends_with(".meta"), -ends_with(".ctd"))

# --- Step 2: find best valid matches (same as before, but without final select) ---
valid_matches <- rosb_original %>%
  left_join(ctd_ref, by = c("Cruise", "Year", "Month"), relationship = "many-to-many") %>%
  mutate(day_diff = abs(Day - Day_ref),
         dep_diff = abs(`Depth [m]` - Dep_ref),
         lat_diff = abs(`Latitude [deg]` - Lat_ref), 
         lon_diff = abs(`Longitude [deg]` - Lon_ref),
         dist_diff = sqrt((lat_diff * 111.3)^2 + (lon_diff * 111.3 * cos(`Latitude [deg]` * pi / 180))^2)) %>%
  filter(dist_diff <= 2 & dep_diff <= 0.25 & day_diff == 0) %>%
  group_by(Cruise, Station, Cast, `Bottle Position`, Year, Month) %>%
  slice_min(order_by = dist_diff + dep_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Cruise, Station, Cast, `Bottle Position`, Year, Month, Day, `Depth [m]`,
         lat_diff, lon_diff, dist_diff, dep_diff, day_diff,
         CDOM_ref, Fluo_ref, BA_ref, Station_ref, Code_ref, Dep_ref)

# --- Step 3: left join back to original, keeping all rows, and coalesce ---
meta.wk.bot <- rosb_original %>%
  left_join(valid_matches,
            by = c("Cruise", "Station", "Cast", "Bottle Position", "Year", "Month", "Day", "Depth [m]")) %>%
  mutate(`CDOM fluorescence [mg/m^3]` = if_else(Year == 2023 & !is.na(CDOM_ref), CDOM_ref, `CDOM fluorescence [mg/m^3]`),
         `Chlorophyll-a fluorescence [mg/m^3]` = if_else(Year == 2023 & !is.na(Fluo_ref), Fluo_ref, `Chlorophyll-a fluorescence [mg/m^3]`),
         `Beam attenuation [1/m]` = if_else(!is.na(BA_ref), BA_ref, `Beam attenuation [1/m]`),
         `Approx Station` = coalesce(Station_ref, `Approx Station`),
         `Approx Code` = coalesce(Code_ref, `Approx Code`),
         `Approx Depth` = coalesce(Dep_ref, `Approx Depth`)) %>%
  select(-CDOM_ref, -Fluo_ref, -BA_ref, -Station_ref, -Code_ref, -Dep_ref)

# range(meta.wk.bot$dist_diff, na.rm = TRUE); median(meta.wk.bot$dist_diff, na.rm = TRUE) # all
# range(meta.wk.bot$dep_diff, na.rm = TRUE); median(meta.wk.bot$dep_diff, na.rm = TRUE) # all
# 
# # run after WK_2_TAX
# meta.wk.bot %>% filter(plate.geno %in% meta.wk.small$plate.geno) %>% pull(dist_diff) %>% {list(range = range(., na.rm = TRUE), median = median(., na.rm = TRUE))} # I use
# meta.wk.bot %>% filter(plate.geno %in% meta.wk.small$plate.geno) %>% pull(dep_diff) %>% {list(range = range(., na.rm = TRUE), median = median(., na.rm = TRUE))} # I use

##

### many FT and two special samples (BI-06 (75 m) in 2023 and STN-13-Z2 (1 m) in 2022) were paired w/ AUTO CTD ###

# BI-06 (75 m) is a non-rosette Niskin sample, but matches well with a standalone CTD cast by the same name made on the same day in close proximity,
# BI-06 standalone CTD cast however reached only 63 m

# At first, I picked manually the ROS-A stations that were sampled on the same DAY as & in close proximity to my FT samples (Approx Station),
# as well as specified the exact DEPTH I want my sample to match w/ in ROS-A data (e.g., Depth = 2; Approx Depth = 1.982)
# However, later, I also double checked my pick using the commented out code to make sure the distance difference is < 1.5 km
meta.wk.auto <- meta.wk %>%
  filter(`Approx Code` == "ROS-A") %>%
  left_join(ctd %>%
              filter(Code == "ROS-A") %>% 
              mutate(`Approx Code` = Code, `Approx Station` = Station, `Approx Depth` = `Depth [m]`) %>%
              select(-Code, -Station, -`Depth [m]`,
                     -`Bottle Position`, -Cast, -`Bot. Depth [m]`, -`Mixed Layer Depth [m]`,
                     -Hour, -Minute, -Second,
                     ROSA_Lat = `Latitude [deg]`, ROSA_Lon = `Longitude [deg]`),
            by = join_by(`Approx Code`, Cruise, `Approx Station`, `Approx Depth`, Year, Month, Day)) %>%
  # report difference in distance and depth for sanity
  mutate(dep_diff = abs(`Depth [m]` - `Approx Depth`),
         dist_diff = sqrt((abs(`Latitude [deg]` - ROSA_Lat) * 111.3)^2 +
                          (abs(`Longitude [deg]` - ROSA_Lon) * 111.3 * cos(`Latitude [deg]` * pi / 180))^2))

# median(meta.wk.auto %>% filter(Code == "FT") %>% .$dist_diff) # all
# 
# # run after WK_2_TAX
# meta.wk.auto %>% filter(plate.geno %in% meta.wk.small$plate.geno, Code == "FT") %>% pull(dist_diff) %>% median(na.rm = TRUE) # I use
# meta.wk.auto %>% filter(plate.geno %in% meta.wk.small$plate.geno, Code == "FT") %>% pull(dist_diff) %>% range(na.rm = TRUE) # I use

# ### This code basically checks if any ROS-A stations are close enough to the FT samples to be used for approximation of parameters ###
# ### This is to validate the choice of the Approx Station used above, which were originally picked manually ###
# ft_approx <- meta.wk %>% filter(`Approx Code` == "ROS-A") %>%
#   left_join(ctd %>% filter(Code == "ROS-A") %>% select(Cruise, Year, Month, Day, `Approx Station` = Station, Depth_roza = `Depth [m]`,
#                                                        Lat_roza = `Latitude [deg]`, Lon_roza = `Longitude [deg]`),
#             by = c("Cruise", "Year", "Month", "Day"), relationship = "many-to-many") %>%
#   mutate(lat_diff = abs(`Latitude [deg]` - Lat_roza),
#          lon_diff = abs(`Longitude [deg]` - Lon_roza),
#          dep_diff = abs(`Depth [m]` - Depth_roza),
#          dist_km = sqrt((lat_diff * 111.3)^2 + (lon_diff * 111.3 * cos(`Latitude [deg]` * pi / 180))^2)) %>%
#   # Keep only matches within your thresholds
#   filter(dist_km <= 2) %>%
#   # For each FT sample, keep the closest match (smallest sum of distance + depth diff)
#   group_by(Cruise, Station, Cast, `Bottle Position`, Year, Month, Day, Hour, Minute) %>%
#   slice_min(order_by = dist_km + dep_diff, n = 1, with_ties = FALSE) %>%
#   ungroup()
# identical(ft_approx$`Approx Station.x`, ft_approx$`Approx Station.y`) # should be true

##

### match "left-over" NISK w/ CastAway CTD (these are only zodiac samples; they also could not find an appropriate match w/ any close-by AUTO CTD) ###
meta.wk.nisk <- meta.wk %>%
  filter(`Approx Code` == "NISK") %>% mutate (`Approx Depth` = round(`Approx Depth`, 4)) %>%
  left_join(ctd %>%
              filter(Code == "NISK") %>% mutate(`Approx Code` = Code, `Approx Station` = Station, `Approx Depth` = round(`Depth [m]`, 4)) %>%
              select (- Code, -Station, -`Depth [m]`,
                      -`Bottle Position`, -Cast, -`Latitude [deg]`, - `Longitude [deg]`, -`Bot. Depth [m]`, -`Mixed Layer Depth [m]`,
                      - Minute, -Second),
            by = join_by(`Approx Code`, Cruise, `Approx Station`, `Approx Depth`, Year, Month, Day, Hour))

# In all cases, CastAway measurements were taken at the exact same HOUR and at the exact same Lat/Long the left-over NISK samples were taken;
# the MRM-associated CastAway lacks Lat/Long, but we know that it was taken on Aug 15, 2023 at ~12:30 pm (low tide at the Moose River).

# Upon review of the CastAway data, I found the most optimal way to encode the match between the left-over NISK samples and CastAway to be by
# using precise timing (DAY and HOUR) & specifying the exact DEPTH I want my sample to match w/ in CastAway data (e.g., Depth = 1; Approx Depth = 1.194960).

##

### approximate "left-over" FT w/ TSG (the ones that could not find a match w/ AUTO CTD) ###

meta.wk.ft <- meta.wk %>% filter(`Approx Code` == "FT")              # samples requiring FT data
ctd_ft <- ctd %>% filter(Code == "FT") %>%                           # TSG readings at FT
  mutate (`Approx Code` = Code, `Approx Station` = Station, `Approx Depth` = `Depth [m]`) %>%
  select(`Approx Code`, Cruise, `Approx Station`, `Approx Depth`,
         # FT has high resolution - readings are taken every 10 sec, so technically should work with Month/Day/Hour precision beside Lat/Long
         Year, Month, Day, Hour,
         `Latitude [deg]`, `Longitude [deg]`, `Depth [m]`,
         `Salinity [PSU]`, `Temperature [deg C]`)
meta.wk.ft <- merge(meta.wk.ft, ctd_ft,                              # merge the two and take the closest match
                    by = c("Approx Code", "Cruise", "Approx Station", "Approx Depth", "Year", "Month", "Day", "Hour")) %>%
  mutate(lon_diff = abs(`Longitude [deg].x` - `Longitude [deg].y`),  # depth is the same in all FT, so no need to calc. the dif-ce
         lat_diff = abs(`Latitude [deg].x` - `Latitude [deg].y`),
         dist_diff = sqrt((lat_diff * 111.3)^2 + (lon_diff * 111.3 * cos(`Latitude [deg].x` * pi / 180))^2)) %>%
  group_by(Cruise, Station, Year, Month, Day, Hour) %>%                        # don't need to group by YMD as we added by YMD, but let it be for consistency
  slice_min(order_by = dist_diff, n = 1, with_ties = FALSE) %>% ungroup() %>%  # a slight dif-ce in distance is due to bad Lat/Long precision in the sampling loc-n
  # result: recorded at practically the same time and loc-n: same day and hour and w/in 6-75 m of the sampling loc-n
  mutate(`Latitude [deg]` = `Latitude [deg].y`, `Longitude [deg]` = `Longitude [deg].y`, `Depth [m]` = `Depth [m].y`) %>%
  select(-contains(".x"), -contains(".y"))

##

### complete data ###
meta.wk.complete <- bind_rows(meta.wk.bot,
                              meta.wk.auto,
                              meta.wk.nisk,
                              meta.wk.ft) %>%
  transmute(Cruise, Year, Month, Day, Hour, Minute, Second,
            Code, `Approx Code`,
            Station, `Approx Station`, `Latitude [deg]`, `Longitude [deg]`,
            `Depth [m]`, `Bot. Depth [m]`,
            genohub_id, plate.geno, dna_conc,
            `Salinity [PSU]`, `Temperature [deg C]`, `Potential Temperature [deg C]`, `Density [kg/m^3]`, `Sigma-theta [kg/m^3]`,
            `Oxygen [umol/kg]`,
            `PAR [umol photons/m^2/sec]`, `Light attenuation coefficient [1/m]`, `Beam attenuation [1/m]`,
            `CDOM fluorescence [mg/m^3]`, `Chlorophyll-a fluorescence [mg/m^3]`) %>%
  st_as_sf(coords = c("Longitude [deg]", "Latitude [deg]"), crs = 4326, remove = FALSE) %>%
  mutate("Distance to coast [km]" = round(st_distance(., st_read("ne_10m/ne_10m_coastline.shp")) %>%
           apply(1, min) %>% as.numeric() / 1000, 0)) %>%
  st_drop_geometry()


# save to csv to upload downstream
write.csv(meta.wk.complete, "OUTPUT/meta.wk.complete.csv", row.names = FALSE, fileEncoding = "UTF-8")


#################################################################################
### Explore the distribution of parameters to make sure there are no oddities ###

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

# # Targeted exploration
# salinity change at CTD_STN_13_Z2_B
# with(ctd %>% filter(Station == "CastAway", Year == 2022, Month == 8, Day == 11), plot(`Depth [m]`*(-1) ~ `Salinity [PSU]`))
# with(ctd %>% filter(Station == "CTD_STN_13_Z2_B"), points(`Depth [m]`*(-1) ~ `Salinity [PSU]`, col = "red"))
# ctd %>% filter(Station == "CastAway", Year == 2022, Month == 8, Day == 11) %>% View()
# 
# ctd %>% filter(Code == c("ROS-B", "ROS-A", "ROS-D")) %>% View()
# 
# correlation of fluorescence with other parameters - none
# par(mfrow = c(2, 2))
# with(meta.wk.complete, plot(`CDOM fluorescence [mg/m^3]` ~ `Beam attenuation [1/m]`))
# with(meta.wk.complete, plot(`CDOM fluorescence [mg/m^3]` ~ `Light attenuation coefficient [1/m]`))
# with(meta.wk.complete, plot(`Chlorophyll-a fluorescence [mg/m^3]` ~ `Beam attenuation [1/m]`))
# with(meta.wk.complete, plot(`Chlorophyll-a fluorescence [mg/m^3]` ~ `Light attenuation coefficient [1/m]`))
# par(mfrow = c(1, 1))
# 
# View(castaway %>% filter (Year == 2022, Month == 8, Day %in% c(11, 17))) # the Maquatua River at low (A11, 10 pm) and high (A17, 9 pm) tides; sampled at low tide
# View(castaway %>% filter (Year == 2023, Month == 8, Day == 15))          # near the Moose River at low (noon) and high (5 pm) tides; sampled at low tide
