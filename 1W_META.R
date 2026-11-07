# 1------------------------------------------------------------------------
# Purpose:
#   Compile metadata.
#   Compile all environmental data (CTD, TSG, CastAway) and sampling/DNA-work information to produce a table
#   that includes all available ("complete") metadata for the DNA samples collected during the WK2022 and WK2023 expeditions.
#   Derived parameters (mixed‑layer depth, distance to coast, and geo_loc_name for MIMARKS) are also calculated.
#
# Required external files:
#   CTD/excel_2022_autonomous_ctd_downcasts.csv (2022 standalone CTD)
#   CTD/2023_autonomous_ctd_downcasts.csv       (2023 standalone CTD)
#   CTD/excel_2022_rosette_ctd_downcasts.csv    (2022 downcast of CTD+bottles)
#   CTD/excel_2023_rosette_ctd_downcasts.csv    (2023 downcast of CTD+bottles)
#   CTD/2022_rosette_bottle_data.csv            (2022 upcast of CTD+bottles)
#   CTD/2023_rosette_bottle_data.csv            (2023 upcast of CTD+bottles)
#   CTD/tsg_2022_merged.csv                     (2022 thermosalinograph)
#   CTD/WK_TSG_data_2023/                       (2023 thermosalinograph)
#   CTD/JBE_2022_CastAway_CC1744004_Data/       (2022 CastAway)
#   CTD/JBE_2023_CastAway_CC1744004_Data/       (2023 CastAway)
#   DNA/meta.wk.genohub.csv                     (information on DNA extraction/analysis w/ genohub ID )
#   DNA/meta.wk.report.approx.csv               (reviewed and cross-referenced information on the samples collected by the microbial team)
#   EXT/ne_10m/                                 (includes the coastline shape file and associated files)
#
# Output:
#   EXT/meta.wk.complete.rds                    (final metadata table; input for 2W_TAX)
#   OUT/1W_META/ctd.odv.csv                     (exploration: CTD/CastAway profiles + DNA sample positions for Ocean Data View, to check station pairing)
#   OUT/1W_META/wk_parameters_hist.pdf          (exploration: histograms of all numeric variables to check for oddities)
# -------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(sf)
library(lubridate)
library(gsw)

##########################
### Compiling CTD data ###
# CTD_AUTO + CTD_DOWN + CTD_BOT + TSG + CastAway 

#
## Autonomous CTD
ctd.auto <- bind_rows(read.csv("CTD/excel_2022_autonomous_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                        rename_with(~"Depth [m]", contains("Depth [salt water, m]")),
                      read.csv("CTD/2023_autonomous_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                        rename_with(~"Depth [m]", contains("Depth [salt water, m]"))) %>%
  
  # Fix the variables
  mutate(Sampling = "ROS-A", Code = 1,
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
  mutate(`Oxygen [umol/kg]` = if_else(Cruise == "WK2022", NA_real_, `Oxygen [umol/kg]`))    # 2022: the O2 data from the auto CTD is wrong (CTD was dropped)


#
## Rosette downcast CTD
ctd.ros_down <- bind_rows(read.csv("CTD/excel_2022_rosette_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE),
                          read.csv("CTD/excel_2023_rosette_ctd_downcasts.csv", stringsAsFactors = FALSE, check.names = FALSE)) %>%
  
  # Fix the variables
  mutate(Sampling = "ROS-D", Code = 2,
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
ctd.ros_bot <- bind_rows(read.csv("CTD/2022_rosette_bottle_data.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                           mutate(Cruise = "WK2022"),
                         read.csv("CTD/2023_rosette_bottle_data.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                           mutate(Cruise = "WK2023")) %>%
  
  # Fix the variables
  mutate(Sampling = "ROS-B", Code = 3,
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
folders <- list(WK2022 = "CTD/JBE_2022_CastAway_CC1744004_Data",
                WK2023 = "CTD/JBE_2023_CastAway_CC1744004_Data")

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
           Sampling = "NISK", Code = 4,
           Year = year(cast_time), Month = month(cast_time), Day = day(cast_time),
           Hour = hour(cast_time), Minute = minute(cast_time), Second = as.numeric(format(cast_time, "%S")))}

# Applying the function and binding (some CastAway don't have Lat/Long but this is OK b/c we know where they were based on time)
castaway <- bind_rows(lapply(list.files(folders$WK2022, "\\.csv$", full.names = TRUE), process_castaway, cruise = "WK2022"),
                      lapply(list.files(folders$WK2023, "\\.csv$", full.names = TRUE), process_castaway, cruise = "WK2023"))


#
## TSG for FT
tsg <- bind_rows(read.csv("CTD/tsg_2022_merged.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
                   mutate(Cruise = "WK2022", date_time = ymd_hms(date_time)),
                 list.files("CTD/WK_TSG_data_2023", pattern = "\\.xlsx$", full.names = TRUE) %>% map_dfr(read_xlsx) %>% mutate(Cruise = "WK2023")) %>%
  rename_with(~ "Latitude [deg]",  .cols = matches("^lat", ignore.case = TRUE)) %>%
  rename_with(~ "Longitude [deg]", .cols = matches("^lon", ignore.case = TRUE)) %>%
  rename_with(~ "Salinity [PSU]",  .cols = matches("^sal", ignore.case = TRUE)) %>%
  rename_with(~ "Temperature [deg C]", .cols = matches("^t0", ignore.case = TRUE)) %>%
  mutate(`Depth [m]` = 2,
         Station = "TSG",
         Sampling = "FT", Code = 5,
         Year = year(date_time),
         Month = month(date_time),
         Day = day(date_time),
         Hour = hour(date_time),
         Minute = minute(date_time),
         Second = as.numeric(format(date_time, "%S"))) %>% select(-timeQ, -date_time, -`c0S/m`) %>%
  filter(`Salinity [PSU]` < 35)                                                         # some Salinity is off


#
## Put together all the oceanographic data and derive MLD from the CTD profiles (autonomous and downcast)
ctd <- bind_rows(ctd.auto, ctd.ros_down, ctd.ros_bot, castaway, tsg) %>% filter(!is.na(`Salinity [PSU]`)) # join instrumental readings and rmv rows w/out S


#
## putting together metadata for our microbial samples based on all the info available in the reports and adding genohub_id
meta.wk.genohub <- read.csv("DNA/meta.wk.genohub.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(across(Cruise, ~ gsub("_", "", as.character(.))),                       # remove underscore from Cruise
         across(genohub_id, ~ gsub("_", "-", as.character(.))),                  # adjust genohub_id
         plate.geno = toupper(paste0(plate, "_", genohub_id)),                   # Make plate.geno
         join_var = ifelse(Cruise == "WK2022", `Depth [m]`, Target)) %>%         # join variable: Depth for WK2022, Target for WK2023
  select(-`Depth [m]`, -Target)

meta.wk.report <- read.csv("DNA/meta.wk.report.approx.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(join_var = ifelse(Cruise == "WK2022", `Depth [m]`, Target))             # join variable: Depth for WK2022, Target for WK2023

meta.wk <- meta.wk.genohub %>%
  inner_join(meta.wk.report, by = c("Cruise", "Station", "join_var")) %>%
  separate(Date, into = c("Month", "Day", "Year"), sep = "/", convert = TRUE) %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%      # use the HMS time from your own meta.wk record
  mutate (Second = 00) %>%
  select(-join_var) %>%
  mutate(`Depth [m]` = as.numeric(`Depth [m]`))

#############

## Exploration: save to csv for Ocean Data View (ODV) to assess the proximity of the sampled stations
## and the CTD/CastAway stations used to approximate their metadata
# Code = 0 for the sampled stations & Code 1-4 for CTD and CastAway; TSG (Code 5) is ignored here
write.csv(bind_rows(ctd %>% filter(Code < 5), meta.wk %>% mutate(Code = 0)), "OUT/1W_META/ctd.odv.csv", row.names = FALSE, fileEncoding = "UTF-8")

#############

#
## Add metadata for samples

### add BOTTLE (upcast) CTD ###

# build “ROS-A” (AUTO CTD) lookup table - will use later
ctd_ref <- ctd %>%
  filter(Sampling == "ROS-A") %>%
  transmute(Cruise, Year, Month,
            Day_ref = Day, Station_ref = Station, Sampling_ref = Sampling,
            Lon_ref = `Longitude [deg]`, Lat_ref = `Latitude [deg]`, Dep_ref = `Depth [m]`,
            CDOM_ref = `CDOM fluorescence [mg/m^3]`, Fluo_ref = `Chlorophyll-a fluorescence [mg/m^3]`, BA_ref = `Beam attenuation [1/m]`)

# --- Step 1: get original ROS-B data (before any approximation) ---
rosb_original <- meta.wk %>% filter(`Pairing` == "ROS-B") %>%
  left_join(ctd %>% filter(Sampling == "ROS-B") %>% select(-`Bot. Depth [m]`, -Hour, -Minute, -Second),
            by = join_by(Sampling, Cruise, Station, Cast, `Bottle Position`, Year, Month, Day), suffix = c(".meta", ".ctd")) %>%
  mutate(`Latitude [deg]`  = `Latitude [deg].ctd`, `Longitude [deg]` = `Longitude [deg].ctd`, `Depth [m]` = `Depth [m].ctd`) %>%
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
         CDOM_ref, Fluo_ref, BA_ref, Station_ref, Sampling_ref, Dep_ref)

# --- Step 3: left join back to original, keeping all rows, and coalesce ---
meta.wk.bot <- rosb_original %>%
  left_join(valid_matches,
            by = c("Cruise", "Station", "Cast", "Bottle Position", "Year", "Month", "Day", "Depth [m]")) %>%
  mutate(`CDOM fluorescence [mg/m^3]` = if_else(Year == 2023 & !is.na(CDOM_ref), CDOM_ref, `CDOM fluorescence [mg/m^3]`),
         `Chlorophyll-a fluorescence [mg/m^3]` = if_else(Year == 2023 & !is.na(Fluo_ref), Fluo_ref, `Chlorophyll-a fluorescence [mg/m^3]`),
         `Beam attenuation [1/m]` = if_else(!is.na(BA_ref), BA_ref, `Beam attenuation [1/m]`),
         `Paired Station` = coalesce(Station_ref, `Paired Station`),
         `Pairing` = coalesce(Sampling_ref, `Pairing`),
         `Approximate Depth` = coalesce(Dep_ref, `Approximate Depth`)) %>%
  select(-CDOM_ref, -Fluo_ref, -BA_ref, -Station_ref, -Sampling_ref, -Dep_ref)



##

### many FT and two special samples (BI-06 (75 m) in 2023 and STN-13-Z2 (1 m) in 2022) were paired w/ AUTO CTD ###

# BI-06 (75 m) is a non-rosette Niskin sample, but matches well with a standalone CTD cast by the same name made on the same day in close proximity,
# BI-06 standalone CTD cast however reached only 63 m

# At first, I picked manually the ROS-A stations that were sampled on the same DAY as & in close proximity to my FT samples (Paired Station),
# as well as specified the exact DEPTH I want my sample to match w/ in ROS-A data (e.g., Depth = 2; Approximate Depth = 1.982)
# However, later, I also double checked my pick using the commented out code to make sure the distance difference is < 1.5 km
meta.wk.auto <- meta.wk %>%
  filter(`Pairing` == "ROS-A") %>%
  left_join(ctd %>%
              filter(Sampling == "ROS-A") %>% 
              mutate(Pairing = Sampling, `Paired Station` = Station, `Approximate Depth` = `Depth [m]`) %>%
              select(-Sampling, -Station, -`Depth [m]`, -`Bottle Position`, -Cast, -`Bot. Depth [m]`, -Hour, -Minute, -Second,
                     ROSA_Lat = `Latitude [deg]`, ROSA_Lon = `Longitude [deg]`),
            by = join_by(`Pairing`, Cruise, `Paired Station`, `Approximate Depth`, Year, Month, Day)) %>%
  # report difference in distance and depth for sanity
  mutate(dep_diff = abs(`Depth [m]` - `Approximate Depth`),
         dist_diff = sqrt((abs(`Latitude [deg]` - ROSA_Lat) * 111.3)^2 +
                            (abs(`Longitude [deg]` - ROSA_Lon) * 111.3 * cos(`Latitude [deg]` * pi / 180))^2))


# ### This code checks if any ROS-A stations are close enough to the FT samples to be used for approximation of parameters ###
# ### This is to validate the choice of the Paired Station used above, which were originally picked manually ###
# ft_approx <- meta.wk %>% filter(`Pairing` == "ROS-A") %>%
#   left_join(ctd %>% filter(Sampling == "ROS-A") %>% select(Cruise, Year, Month, Day, `Paired Station` = Station, Depth_roza = `Depth [m]`,
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
# identical(ft_approx$`Paired Station.x`, ft_approx$`Paired Station.y`) # should be true

##

### match "left-over" NISK w/ CastAway CTD (these are only zodiac samples; they also could not find an appropriate match w/ any close-by AUTO CTD) ###
meta.wk.nisk <- meta.wk %>%
  filter(`Pairing` == "NISK") %>% mutate (`Approximate Depth` = round(`Approximate Depth`, 4)) %>%
  left_join(ctd %>%
              filter(Sampling == "NISK") %>% mutate(Pairing = Sampling, `Paired Station` = Station, `Approximate Depth` = round(`Depth [m]`, 4)) %>%
              select (- Sampling, -Station, -`Depth [m]`, -`Bottle Position`, -Cast, -`Latitude [deg]`, - `Longitude [deg]`, -`Bot. Depth [m]`,
                      - Minute, -Second), by = join_by(`Pairing`, Cruise, `Paired Station`, `Approximate Depth`, Year, Month, Day, Hour))

# In all cases, CastAway measurements were taken at the exact same HOUR and at the exact same Lat/Long the left-over NISK samples were taken;
# the MRM-associated CastAway lacks Lat/Long, but we know that it was taken on Aug 15, 2023 at ~12:30 pm (low tide at the Moose River).

# For practicality, I matched left-over NISK samples and CastAway measurements by time (DAY and HOUR) & exact DEPTH, where "Approximate Depth" is
# the depth closest to the reported "Depth" the sample was taken at (e.g., reported Depth = 1; based on CastAway, Approximate Depth = 1.194960).
# > meta.wk.nisk$`Depth [m]`
# [1] 1 1 1 4
# > meta.wk.nisk$`Approximate Depth`
# [1] 1.1950 1.0502 1.0541 4.0511


### approximate "left-over" FT w/ TSG (the ones that could not find a match w/ AUTO CTD) ###
meta.wk.ft <- meta.wk %>% filter(`Pairing` == "FT")                      # samples requiring FT data
ctd_ft <- ctd %>% filter(Sampling == "FT") %>%                           # TSG readings at FT
  mutate (Pairing = Sampling, `Paired Station` = Station, `Approximate Depth` = `Depth [m]`) %>%
  select(Pairing, Cruise, `Paired Station`, `Approximate Depth`,
         # FT has high resolution - readings are taken every 10 sec, so technically should work with Month/Day/Hour precision beside Lat/Long
         Year, Month, Day, Hour,
         `Latitude [deg]`, `Longitude [deg]`, `Depth [m]`,
         `Salinity [PSU]`, `Temperature [deg C]`)
meta.wk.ft <- merge(meta.wk.ft, ctd_ft,                              # merge the two and take the closest match
                    by = c("Pairing", "Cruise", "Paired Station", "Approximate Depth", "Year", "Month", "Day", "Hour")) %>%
  mutate(lon_diff = abs(`Longitude [deg].x` - `Longitude [deg].y`),  # depth is the same in all FT, so no need to calc. the dif-ce
         lat_diff = abs(`Latitude [deg].x` - `Latitude [deg].y`),
         dist_diff = sqrt((lat_diff * 111.3)^2 + (lon_diff * 111.3 * cos(`Latitude [deg].x` * pi / 180))^2)) %>%
  group_by(Cruise, Station, Year, Month, Day, Hour) %>%                        # don't need to group by YMD as we added by YMD, but let it be for consistency
  slice_min(order_by = dist_diff, n = 1, with_ties = FALSE) %>% ungroup() %>%  # a slight dif-ce in distance due to bad Lat/Long precision in the samp. loc-n
  # result: recorded at practically the same time and loc-n: same day and hour and w/in 6-75 m of the sampling loc-n
  mutate(`Latitude [deg]` = `Latitude [deg].y`, `Longitude [deg]` = `Longitude [deg].y`, `Depth [m]` = `Depth [m].y`) %>%
  select(-contains(".x"), -contains(".y"))

##

### complete data ###

################################
### Final WK sample metadata ###
################################

sp_crs <- 3347
coast <- st_read("EXT/ne_10m/ne_10m_coastline.shp", quiet = TRUE) %>% st_transform(sp_crs)
marine <- st_read("EXT/ne_10m/ne_10m_geography_marine_polys.shp", quiet = TRUE) %>%
  select(marine = name) %>% st_transform(sp_crs)

meta.wk.complete <- bind_rows(meta.wk.bot, meta.wk.auto, meta.wk.nisk, meta.wk.ft) %>%
  
  # Standardize names that differ in the WK source files
  dplyr::rename(`Bottom depth [m]` = `Bot. Depth [m]`, `Potential temperature [deg C]` = `Potential Temperature [deg C]`,
                Salinity = `Salinity [PSU]`, Bottle = `Bottle Position`, `Sample label` = Sample, genohub.id = genohub_id) %>%
  
  mutate(across(any_of(c("Latitude [deg]", "Longitude [deg]", "Depth [m]", "Temperature [deg C]", "Salinity",
                         "Potential temperature [deg C]", "Density [kg/m^3]", "Sigma-theta [kg/m^3]",
                         "Oxygen [umol/kg]")), as.numeric),
         
         # Preserve the original WK sampling-station name before chronological renumbering
         `Original station` = as.character(Station),
         
         # Reconstruct sampling time and fractional day of year
         Time = sprintf("%02d:%02d:%02d", Hour, Minute, Second),
         DOY = yday(make_date(Year, Month, Day)) + (Hour * 3600 + Minute * 60 + Second) / 86400,
         
         # Temporary TEOS-10 variables; pressure is derived only for these calculations
         .P = gsw_p_from_z(-`Depth [m]`, `Latitude [deg]`),
         .SA = gsw_SA_from_SP(Salinity, .P, `Longitude [deg]`, `Latitude [deg]`),
         .CT = gsw_CT_from_t(.SA, `Temperature [deg C]`, .P),
         
         # Retain supplied values and fill missing values using TEOS-10
         `Potential temperature [deg C]` = coalesce(`Potential temperature [deg C]`, gsw_pt_from_CT(.SA, .CT)),
         `Density [kg/m^3]` = coalesce(`Density [kg/m^3]`, gsw_rho(.SA, .CT, .P)),
         `Sigma-theta [kg/m^3]` = coalesce(`Sigma-theta [kg/m^3]`, gsw_sigma0(.SA, .CT)),
         `Sigma-t [kg/m^3]` = gsw_rho(.SA, gsw_CT_from_t(.SA, `Temperature [deg C]`, 0), 0) - 1000,
         `Dissolved oxygen [uM]` = `Oxygen [umol/kg]` * `Density [kg/m^3]` / 1000,
         `AOU [umol/kg]` = gsw_O2sol(.SA, .CT, .P, `Longitude [deg]`, `Latitude [deg]`) - `Oxygen [umol/kg]`,

         # Known WK filtration conditions
         `Sequential filtration` = 0L, `Pore size [μm]` = 0.22) %>%
  
  # MLD from the complete rosette downcast profile (Peralta-Ferriz & Woodgate 2015 threshold, 0.1 kg/m^3 relative to the shallowest bin):
  # the maximum depth to which sigma-theta stays CONTINUOUSLY within the threshold, i.e. the bin above the first exceedance.
  # For monotonic profiles, this equals their "maximum depth within the threshold"; density inversions below the pycnocline
  # (profiles they rejected in QC) do not extend the mixed layer, e.g. MLD = 1.49 m at 2022 STN-17 cast 22 despite a deeper 1.5 kg/m^3 inversion.
  # If the threshold is never exceeded, the water column is mixed to the bottom of the cast (MLD = deepest bin).
  left_join(ctd %>%
              filter(Sampling == "ROS-D", !is.na(`Sigma-theta [kg/m^3]`)) %>%
              mutate(Cast = suppressWarnings(as.integer(as.character(Cast)))) %>%
              group_by(Cruise, Station, Cast, Year, Month, Day) %>%
              arrange(`Depth [m]`, .by_group = TRUE) %>%
              summarise(`Mixed layer depth [m]` = {
                exceed <- which(`Sigma-theta [kg/m^3]` > `Sigma-theta [kg/m^3]`[1] + 0.1)
                if (length(exceed)) `Depth [m]`[exceed[1] - 1] else max(`Depth [m]`)
              }, .groups = "drop"),
            by = c("Cruise", "Station", "Cast", "Year", "Month", "Day")) %>%
  
  # Standardized geographic variables
  st_as_sf(coords = c("Longitude [deg]", "Latitude [deg]"), crs = 4326, remove = FALSE) %>% st_transform(sp_crs) %>%
  mutate(`Distance to coast [km]` = {
    i <- st_nearest_feature(geometry, coast)
    round(as.numeric(st_distance(geometry, st_geometry(coast)[i], by_element = TRUE)) / 1000)
  },
  geo_loc_name = {
    i <- st_nearest_feature(geometry, marine)
    d <- as.numeric(st_distance(geometry, st_geometry(marine)[i], by_element = TRUE))
    if_else(d <= 5000, paste0("Canada: ", marine$marine[i]), "Canada")
  }) %>%
  st_drop_geometry() %>%
  
  # Number stations by their first sampling time across WK2022 and WK2023
  {
    b <- .
    station_key <- b %>%
      mutate(.sample_time = make_datetime(Year, Month, Day, Hour, Minute, Second, tz = "UTC")) %>%
      arrange(.sample_time, Cruise, `Original station`) %>%
      distinct(Cruise, `Original station`) %>%
      mutate(Station = as.character(row_number()))
    
    b %>% select(-Station) %>% left_join(station_key, by = c("Cruise", "Original station"))
  } %>%
  
  # Types shared with the other project metadata
  mutate(across(any_of(c("Cruise", "Station", "Original station", "Target", "Sample label", "Time",
                         "plate", "genohub.id", "plate.geno")), as.character),
         across(any_of(c("Cast", "Bottle", "Year", "Month", "Day", "Hour", "Minute", "Second", "Sequential filtration")),
                \(x) suppressWarnings(as.integer(as.character(x))))) %>%
  
  # Final WK output: unavailable columns are omitted rather than created as NA
  select(any_of(c("Cruise", "Station", "Original station", "Sampling", "Cast", "Latitude [deg]", "Longitude [deg]", "geo_loc_name",
                  "Bottom depth [m]", "Target", "Bottle", "Sample label", "Year", "Month", "Day", "Hour", "Minute", "Second",
                  "Time", "DOY", "Depth [m]", "Potential temperature [deg C]", "Temperature [deg C]", "Salinity",
                  "Bottle salinity", "Density [kg/m^3]", "Sigma-t [kg/m^3]", "Sigma-theta [kg/m^3]",
                  "Chlorophyll-a fluorescence [mg/m^3]", "CDOM fluorescence [mg/m^3]", "Beam attenuation [1/m]",
                  "Nitrate (CTD) [uM]", "Nitrate [uM]", "Nitrite [uM]", "Phosphate [uM]", "Silicate [uM]",
                  "Dissolved oxygen [uM]", "AOU [umol/kg]", "Mixed layer depth [m]", "Distance to coast [km]",
                  "Sequential filtration", "Prefiltration mesh pore size [μm]", "Pore size [μm]", "Volume filtered [L]",
                  "DNA [mg/L]", "plate", "genohub.id", "plate.geno")))

saveRDS(meta.wk.complete, "EXT/meta.wk.complete.rds")


##################################################################################
### Exploration: distribution of parameters to make sure there are no oddities ###

pdf("OUT/1W_META/wk_parameters_hist.pdf", width = 8.5, height = 12)    # Open PDF
par(mfrow = c(3, 2), mar = c(5, 6, 4, 3))                      # Six plots per page with wider margins
for (col_name in names(meta.wk.complete)[sapply(meta.wk.complete, is.numeric)]) {      # Loop through numeric columns
  x <- na.omit(meta.wk.complete[[col_name]])                                           # Remove missing values
  if (!length(x)) next                                                                 # Skip empty columns
  hist(x, main = NULL, xlab = bquote(bold(.(col_name))), xlim = range(x), col = "lightblue", border = "black")  # Plot histogram
  mtext(paste("Min:", round(min(x), 2)), side = 3, line = 2, adj = 0, cex = 0.8)                                # Add minimum
  mtext(paste("Max:", round(max(x), 2)), side = 3, line = 1, adj = 0, cex = 0.8)                                # Add maximum
  mtext(paste("Mean ± SD:", round(mean(x), 2), "±", round(sd(x), 2)), side = 3, line = 0, adj = 0, cex = 0.8)}  # Add mean and SD
dev.off()                                                      # Close PDF


# # Targeted exploration
# # salinity change at CTD_STN_13_Z2_B
# with(ctd %>% filter(Station == "CastAway", Year == 2022, Month == 8, Day == 11), plot(`Depth [m]`*(-1) ~ `Salinity`))
# with(ctd %>% filter(Station == "CTD_STN_13_Z2_B"), points(`Depth [m]`*(-1) ~ `Salinity`, col = "red"))
# ctd %>% filter(Station == "CastAway", Year == 2022, Month == 8, Day == 11) %>% View()
# 
# ctd %>% filter(Sampling == c("ROS-B", "ROS-A", "ROS-D")) %>% View()
# 
# # correlation of fluorescence with other parameters - none
# with(meta.wk.complete, plot(`CDOM fluorescence [mg/m^3]` ~ `Beam attenuation [1/m]`))
# with(meta.wk.complete, plot(`Chlorophyll-a fluorescence [mg/m^3]` ~ `Beam attenuation [1/m]`))
# 
# View(castaway %>% filter (Year == 2022, Month == 8, Day %in% c(11, 17))) # the Maquatua River at low (A11, 10pm) and high (A17, 9pm) tides; sampled at low
# View(castaway %>% filter (Year == 2023, Month == 8, Day == 15))          # near the Moose River at low (noon) and high (5 pm) tides; sampled at low
