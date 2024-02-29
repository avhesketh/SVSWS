# Amelia Hesketh January 2023
########################### 
# Reading in, cleaning, and visualizing iButton temperature data

# Load packages
pkgs <- c("tidyverse","lubridate","lme4", "car", "stringi", 
          "DHARMa", "patchwork", "glmmTMB", "plotrix", "emmeans")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# Read in the raw temperature data

# cd into the appropriate temperature data folder
setwd("./raw_data/temperature")

# pull all the file names
file_names <- list.files()

# two files have ~$ at the beginning of the names, so manually replace them
# so that read_csv works
file_names[1] <- "SVSWS_D_4_temp_20200819.csv"
file_names[2] <- "SVSWS_D_5_temp_20200819.csv"

# prepare an empty dataframe to which to add all the raw data
ibuttons <- data.frame()

# this loop reads the files and adds them to the 'ibutton' amalgamated dataframe
for (i in 1:length(file_names)){
  temp_data <- read_csv(file_names[i], skip = 14)
  temp_data$info <- file_names[i] # retain file names, which have information about the tile ID
  ibuttons <- rbind(ibuttons, temp_data)
}

setwd("../../")

# We need to account for when PDT/PST shifts occurred since the ibuttons don't know, and
# we will be matching these data up with tide height data later to determine when tiles
# are emersed. These are the relevant date when time changes occurred during experiment.
pdt_2019_end <- ymd_hms("2019-11-03 02:00:00")
pdt_2020_begin <- ymd_hms("2020-03-08 02:00:00")
pdt_2020_end <- ymd_hms("2020-11-01 02:00:00")

# Cleaning the raw data

ibuttons$`Date/Time` <- as.character(ibuttons$`Date/Time`)

ibutton_separate <- ibuttons %>%
  # change date/time to a parse-able format
  mutate(date_time = if_else(info %in% 
                                       c("SVSWS_D_4_temp_20200819.csv", "SVSWS_F_0_temp_20210224.csv",
                                         "SVSWS_D_5_temp_20200819.csv", "SVSWS_A_0_temp_20200814.csv",
                                         "SVSWS_C_0_temp_20200814.csv", "SVSWS_D_13_temp_20190814.csv",
                                         "SVSWS_D_0_temp_20200819.csv"),
                             dmy_hm(substr(`Date/Time`, 3, length(`Date/Time`))), dmy_hms(`Date/Time`))) %>% 
  # remove useless columns
  select(-Unit, -`Date/Time`) %>% 
  # separate out file information into useful columns
  separate(info, into = c("project","block","number","X", "collect"), sep = "_") %>% 
  select(-project, -X) %>% # remove project code and empty data
  rename(temp = Value) %>% 
  mutate(number = as.numeric(number),
         # clean up to create a collection date column
         collect_date = str_remove_all(collect, ".csv")) %>% 
  select(-collect) %>% 
  mutate(collect_date = ymd(collect_date)) %>% # collection date = when the data were read from ibutton
  # here, we take values where times are in pdt that should be in pst (and vice versa) and fix it
  mutate(date_time = if_else(collect_date == "2020-01-22" & date_time >= pdt_2019_end,
                             date_time - hours(1), 
                             if_else(collect_date == "2020-03-14" & date_time >= pdt_2020_begin,
                                     date_time + hours(1), 
                                     if_else(collect_date == "2021-02-24" & date_time >= pdt_2020_end,
                                             date_time - hours(1), date_time), date_time),
                             date_time)) %>% 
  mutate(record_date = date(date_time)) %>% 
  group_by(block, number, date_time, collect_date, record_date) %>% 
  summarize(temp = mean(temp))

# this file has info about where tiles were during initial settlement vs. during the experiment
# tiles were relocated due to wave/log damage to more sheltered areas of shoreline
block_design <- read_csv("./raw_data/design/SVSWS_tilesetup.csv")

# ibuttons were sometimes set up or data were collected before/after field work 
# so as to not bring the laptop out into adverse weather
# thus, we need to filter data out that were not actually collected in the field

ibutton_filter <- ibutton_separate %>% ungroup()  %>% 
  mutate(remove = ifelse((collect_date == "2019-10-18" & date_time > ymd_hms("2019-10-18 3:00:00"))|
                           (collect_date == "2020-03-14" & date_time < ymd_hms("2020-01-22 20:00:00"))|
                           (collect_date == "2020-09-14" & date_time > ymd_hms("2020-09-14 12:00:00"))|
                           (record_date == "2019-06-06" & (date_time > ymd_hms("2019-06-06 10:00:00") & 
                                                             date_time <= ymd_hms("2019-06-06 16:00:00"))) |
                           (date_time >= ymd_hms("2021-02-24 19:00:00")),
                         TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  select(!remove)

# Tiles were moved to new locations due to log damage 
# Here, we divide data into *before* tiles were moved and *after* tiles were moved
# This is to ensure that individual tiles are tracked properly through time.

ibutton_tiles <- ibutton_filter %>% 
  mutate(original_no = ifelse(collect_date <= "2019-06-06", number, NA),
         original_block = ifelse(collect_date <= "2019-06-06", block, NA),
         new_no = ifelse(collect_date > "2019-06-06", number, NA),
         new_block = ifelse(collect_date > "2019-06-06", block, NA))

ibutton_originals <- ibutton_tiles %>% 
  filter(is.na(new_no)) %>% 
  select(-new_no, -new_block) %>% 
  full_join(block_design)

ibutton_position <- ibutton_tiles %>% 
  filter(is.na(original_no)) %>% 
  select(-original_no, -original_block) %>% 
  full_join(block_design) %>%
  full_join(ibutton_originals) %>%
  filter(is.na(date_time)==F) 

# The data are cleaned and tile IDs have been corrected. Now we assign treatments
# rather than using tile colour for year 1.
ibutton_clean <- ibutton_position %>% 
  mutate(trt_y1 = case_when(is.na(colour_y1) ~ "Rock",
                            colour_y1 == "white" ~ "C",
                            colour_y1 == "black" ~ "W"),
         trt_y2 = if_else(treatment == "rock","Rock",
                          treatment)) %>% 
  mutate(hour = hour(date_time)) %>% 
  mutate(date_time = ymd_h(paste(record_date, hour, sep = " "))) %>% 
  rename(date = record_date)

# We also recorded specific tile information about tile angle and aspect, but these didn't
# turn out to be important since they were similar between tiles.
# These variables can be removed.
cols_to_remove <- c("original_angle","new_angle",
                    "new_angle","new_compass","colour_y1",
                    "colour_y2", "survived","original_herb_trt", "hour")

ibutton_clean <- ibutton_clean %>% select(!all_of(cols_to_remove)) %>% filter(is.na(temp) == F)

# Write the clean file...
write_csv(ibutton_clean, "./clean_data/SVSWS_temp_clean.csv")


########################### 
# Cursory examination of temperature data

ibutton_clean <- read_csv("./clean_data/SVSWS_temp_clean.csv")

# Here, we visualize ALL the mean temperature data for the whole experiment, divided into 
# current treatment (white surface = cool, black surface = warm)
temp_summary <- ibutton_clean %>% 
  mutate(current_trt = if_else(date_time < "2020-04-03 16:00:00", trt_y1, if_else(trt_y2 == "Rock", trt_y2,substr(trt_y2,2,2)))) %>% 
  group_by(current_trt, date_time) %>% 
  summarise(mean_temp = mean(temp, na.rm = TRUE))

trt_temps <- ggplot(aes(y = mean_temp, x = date_time, colour = current_trt), 
                    data = temp_summary %>% filter(current_trt %in% c("W", "C", "Rock"))) +
  geom_line(lwd = 0.2) +
  facet_wrap(~current_trt) +
  theme_bw() +
  labs(x = "Time", y = "Average hourly temperature (˚C)", colour = "Treatment") +
  scale_colour_manual(values = c("cadetblue","grey50","firebrick4")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ylim(c(-5,40))

trt_temps
# Note that some data are missing for rock temperatures between 19 July to 19 August 2020.
# iButtons were pre-programmed incorrectly, and started recording too early.

###########################
# Filtering the data for analyses

# We need to isolate times when tiles are out of the water 
# by matching the level of tiles with that of the tide.

# Read in tide data and tidy it up 
tides <- read_csv("./raw_data/design/SVSWS_tides.csv", col_names = FALSE) %>% 
  separate(X1, c("date","blank", "time", "time_zone", "tide_height"), sep = "([\\  \\ ])") %>% 
  select(-blank) %>% 
  unite(c("date", "time"), col = date_time, sep = " ") %>% 
  mutate(tide_height = as.numeric(tide_height), date_time = ymd_hm(date_time)) %>% 
  group_by(date_time) %>% # need to average tide levels for when time changes occur
  summarize(tide_height = mean(tide_height)) %>% ungroup() %>% 
  mutate(date = date(date_time)) %>%   group_by(date) %>% 
  # find the minimum tide height to use as an additional covariate later
  mutate(min_tide = min(tide_height)) %>% ungroup()


# For subsequent analysis of overall temperature means & maxima, need to isolate only data where
# tiles are above water and during daytime low tides (March through October)

# clean up Climate Canada data to get daily sunrise and sunset times
daylight_hours_pdt <- read_csv("./raw_data/design/SVSWS_daylight_hours_pst.rtf",
                               col_names = FALSE, skip = 9) %>% 
  select(1, 4, 6) %>% 
  rename(date = 1, sunrise = 2, sunset = 3) %>% 
  mutate(date = mdy(stri_sub(date, -12,-1))) %>% 
  na.omit() %>% 
  mutate(sunrise_datetime = ymd_hms(paste(date, sunrise, sep = " ")),
         sunset_datetime = ymd_hms(paste(date, sunset, " "))) %>% 
  mutate(sunrise_datetime = if_else((date < date(pdt_2019_end)) | 
                                      (date > pdt_2020_begin & date > pdt_2020_end), 
                                    sunrise_datetime + hours(1), sunrise_datetime),
         sunset_datetime = if_else((date < date(pdt_2019_end)) | 
                                     (date > pdt_2020_begin & date > pdt_2020_end), 
                                   sunset_datetime + hours(1), sunset_datetime)) %>% 
  select(-sunrise, -sunset)

# isolate only low tides where it is also daylight
daylight_tides <- tides %>% 
  left_join(daylight_hours_pdt) %>% 
  filter(date_time > sunrise_datetime & date_time < sunset_datetime) %>% 
  select(date_time, date, tide_height, min_tide)

# we approximated tile heights based on the tide height when temperatures drop 
# as tides rise during peak summertime low tides (averaged over three tides)
# here, we read these tile level data in and clean them up
tile_heights <- read_csv("./raw_data/design/SVSWS_tilelevels.csv") %>% 
  pivot_longer(cols = c(tile_height_1:tile_height_3), values_to = "shore_level_m") %>% 
  group_by(new_block, new_number) %>% 
  summarize(tile_level = mean(shore_level_m)) %>% 
  rename(new_no = new_number)

# define date when treatments swapped from y1 to y2 treatments
swap_date <- ymd_hms("2020-04-03 17:00:00")

# filter temperature data that satisfy our requirements
out_of_water <- ibutton_clean %>% 
  # if date is before treatment swap date, then the "active" treatment is that of y1
  # else it is that of y2
  mutate(treatment = if_else(date_time < swap_date, trt_y1, trt_y2)) %>% 
  mutate(period = if_else(date_time < swap_date, "Year 1", "Year 2")) %>% 
  # only retain those temperatures collected during daytime low tides
  right_join(daylight_tides) %>% 
  left_join(tile_heights) %>% 
  # and only retain temperatures collected when water levels falls below level of the tile
  filter(date > "2019-06-06" & tide_height < tile_level)


# In addition to tides and daylight, daytime cloud cover probably plays an important role
# since passive warming relies on solar irradiance reaching the tiles
# Here, we read in and clean the cloud cover data from NOAA.
cloud_cover <- read_csv("./raw_data/design/SVSWS_cloud_cover.csv") %>% 
  select(DATE, HourlySkyConditions) %>% 
  separate(2, into = c("clouds_layer1", "clouds_layer1_okta", "clouds_layer1_height",
                       "clouds_layer2","clouds_layer2_okta", "clouds_layer2_height",
                       "clouds_layer3", "clouds_layer3_okta", "clouds_layer3_height")) %>% 
  # only consider the highest cloud layer (e.g., if overcast with scattered clouds, consider overcast)
  mutate(cloud_total = if_else(is.na(clouds_layer3)==F, clouds_layer3_okta,
                               if_else(is.na(clouds_layer2) == F, clouds_layer2_okta,
                                       clouds_layer1_okta))) %>% 
  rename(date_time = DATE) %>% 
  select(date_time, cloud_total) %>% 
  mutate(date = date(date_time),
         cloud_total = as.integer(cloud_total)) %>% 
  na.omit() %>% 
  # only retain clouds during daylight hours for daily mean value
  left_join(daylight_hours_pdt) %>% 
  filter(date_time > sunrise_datetime & date_time < sunset_datetime) %>% 
  mutate(date_time = round_date(ymd_hms(date_time), "hour"),
         cloud_total = if_else(cloud_total == 9, NA, cloud_total)) %>% 
  unique() %>% 
  separate(date_time, c("date", "time"), sep = " ") %>% 
  mutate(date = ymd(date)) %>% 
  # create daily mean cloud cover
  group_by(date) %>% summarize(cloud_total = mean(cloud_total, na.rm = T)) %>% 
  ungroup()

# final filter ... join to cloud cover data and retain only data collected
# between 15 June and 31 August of each year for looking @ temperature treatment 
# effectiveness
daily_temps <- out_of_water %>% 
  filter(is.na(temp) == F)  %>% 
  left_join(cloud_cover) %>% 
  select(-tide_height) %>% 
  filter(date >= "2019-06-15" & date < "2019-09-01" |
           date >= "2020-06-15" & date < "2020-09-01")

###########################
# Visualizing temperatures in each treatment

# daily mean temperatures
mean_daily <- daily_temps %>% 
  group_by(treatment, date, period, block, number,
           min_tide, cloud_total) %>% 
  summarize(mean_dt = mean(temp, na.rm = T)) %>% ungroup()

# daily maximum temperatures
max_daily <- daily_temps %>% 
  group_by(treatment, date, period, block, number, 
           min_tide, cloud_total) %>% 
  summarize(mdt = max(temp, na.rm = T)) %>% ungroup()

# generate dataframe for mean daily maximum temperature plot
av_mdt <- max_daily %>% 
  group_by(treatment, period, block) %>% 
  summarize(av_mdt = mean(mdt, na.rm = T), se_mdt = std.error(mdt),
            sd_mdt = sd(mdt)) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "Summer 2019", "Summer 2020"))

# define palettes for colour and shape for plots
pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B", "grey50")
pch.trt <- c(1,1,16,16,16,16,16)

# maximum daily temperature plot
mdt_plot <- ggplot(aes(x = treatment, y = av_mdt, col = treatment, pch = treatment),
                   data = av_mdt) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(width = 0.25, size = 0.8) +
  theme(strip.text = element_text(size = 14)) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  scale_linetype_manual(values = lty.trt) +
  labs(y = "Mean daily maximum temperature (ºC)", 
       x = "Treatment",
       col = "Treatment",
       pch = "Treatment") +
  theme_classic() +
  facet_wrap(~period, scales = "free_x") +
  theme(plot.tag = element_text(face= "bold"))
mdt_plot

# Figure 1 in main text = maximum daily temperature only
png("./figures/Fig1.png", res = 700, width = 7, height = 5, units = "in")
mdt_plot
dev.off()

# generate dataframe for average mean daily temperature plot
av_meandt <- mean_daily %>% 
  group_by(treatment, period, block) %>% 
  summarize(av_mdt = mean(mean_dt, na.rm = T), se_mdt = std.error(mean_dt),
            sd_mdt = sd(mean_dt)) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "Summer 2019", "Summer 2020"))

# mean daily temperature plot
meant_plot <- ggplot(aes(x = treatment, y = av_mdt, col = treatment, pch = treatment),
                   data = av_meandt) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(width = 0.25, size = 0.8) +
  theme(strip.text = element_text(size = 14)) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  scale_linetype_manual(values = lty.trt) +
  labs(y = "Mean temperature (ºC)", 
       x = "Treatment",
       col = "Treatment",
       pch = "Treatment") +
  theme_classic() +
  facet_wrap(~period, scales = "free_x") +
  theme(
        plot.tag = element_text(face= "bold"))
meant_plot

# mean temperature plot for supplement (Figure A4)
png("./figures/FigureA4.png", res = 700, width = 7, height = 5, units = "in")
meant_plot
dev.off()


# Appendix 1 Fig A2: Example trace showing treatment effect over one summer low tide series
# Isolate one fragment of temperature data
sample_data <- ibutton_clean %>% 
  filter(record_date %in% c("2019-07-28","2019-07-29",
                            "2019-07-30","2019-07-31","2019-08-01", "2019-08-02", "2019-08-03")) %>% 
  mutate(current_trt = if_else(date_time < "2020-04-03 16:00:00", 
                               trt_y1, if_else(trt_y2 == "Rock", trt_y2,substr(trt_y2,2,2)))) %>% 
  group_by(current_trt, date_time) %>% 
  summarise(mean_temp = mean(temp, na.rm = TRUE), mean_shore_level = mean(new_shore_level))

# Isolate analogous chunk of tide height data
tides_trace <- tides %>% 
  mutate(record_date = date(date_time)) %>% 
  filter(record_date %in%  c("2019-07-28","2019-07-29",
                             "2019-07-30","2019-07-31","2019-08-01", "2019-08-02", "2019-08-03"))

# Join together the data, indicating when the tide is above the tiles and vice versa
sample_data2 <- sample_data %>% left_join(tides_trace) %>% 
  mutate(above_below = if_else(tide_height >= mean_shore_level, "below", "above")) %>% 
  mutate(current_trt = if_else(current_trt == "C", "Cool", if_else(
    current_trt == "W", "Warm", "Rock")))

sample_trace <- ggplot(aes(x = date_time), data = sample_data2) +
  geom_tile(aes(y = 20, height = 40, fill = above_below), alpha = 0.2) +
  geom_line(aes(y = mean_temp, colour = current_trt), lwd = 0.7) +
  theme_classic() +
  labs(x = "Date", y = "Average hourly temperature (˚C)", colour = "Treatment") +
  scale_colour_manual(values = c("cadetblue","grey50","firebrick4")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  geom_line(aes(y = tide_height*10), lwd = 0.5, lty = "dashed", col = "darkblue") +
  scale_y_continuous(sec.axis = sec_axis(~.*0.1, name="Tide height (m)\n")) +
  scale_fill_manual(values = c("lemonchiffon","lightblue")) +
  guides(fill = "none")
sample_trace


png("./figures/FigA2.png", width = 12, height = 5, units = "in",
    res = 700)
sample_trace
dev.off()


##########################
# Models of mean & maximum daily temperature by treatment

# Model of maximum daily temperature in y1

# Isolate appropriate data to model (collected in Year 1 of the experiment)
trt.y1 <- max_daily %>% filter(period == "Year 1") %>% mutate(treatment = factor(treatment))

# Include main effects of treatment, cloud cover, and shore level of low tide
mod.amdt1 <- glmmTMB(mdt ~ treatment + cloud_total +  min_tide +
                       # nested random effect for individual tile within block
                       # and disperson formula to account for unequal variance between treatments
                    (1|block/number), dispformula = ~treatment,
                    data = trt.y1)

# generate residuals and visualize
res.mdt.1 <- simulateResiduals(mod.amdt1)
plot(res.mdt.1) # violation of normality in residuals

# model summary and Anova (Type II)
summary(mod.amdt1)
Anova(mod.amdt1)

# repeat for model in second year
trt.y2 <- max_daily2 %>% filter(period == "Year 2") %>% mutate(treatment = factor(treatment))

# model does not require additional covariates
mod.amdt2 <- glmmTMB(mdt ~ treatment +
                    (1|block/number), dispformula = ~treatment,
                  data = trt.y2)

res.mdt.2 <- simulateResiduals(mod.amdt2)
plot(res.mdt.2) # assumption of normality again violated
summary(mod.amdt2)
Anova(mod.amdt2)

# look at significance of contrasts between treatment groups using emmeans
emm1 <- emmeans(mod.amdt1, specs = pairwise ~ treatment)
emm1$emmeans
emm1$contrasts

emm2 <- emmeans(mod.amdt2, specs = pairwise ~ treatment)
emm2$emmeans
emm2$contrasts

# look at overall mean ± SE for treatments and bedrock
overall_mean <-  daily_temps %>% 
  group_by(treatment, date, period, block, number) %>% 
  summarize(mdt = max(temp, na.rm = T)) %>% ungroup() %>% 
  mutate(treatment = as.factor(treatment)) %>% 
  group_by(treatment, period) %>% 
  summarize(mmdt = mean(mdt), se_mdt = std.error(mdt)) %>% ungroup()


# same modeling scheme, but for mean temperatures (Appendix 2)
trt.y1 <- mean_daily %>% filter(period == "Year 1") %>% mutate(treatment = factor(treatment))

mod.meant1 <- glmmTMB(mean_dt ~ treatment +  cloud_total +  min_tide +
                        (1 | block/number), dispformula = ~treatment, 
                      data = trt.y1)

res.meant.1 <- simulateResiduals(mod.meant1)
plot(res.meant.1)
summary(mod.meant1)
Anova(mod.meant1)

trt.y2 <- mean_daily %>% filter(period == "Year 2") %>% mutate(treatment = factor(treatment),
                                                               julian_day = yday(date))

mod.meant2 <- glmmTMB(mean_dt ~ treatment +  
                        (1|block/number), dispformula = ~treatment, 
                      data = trt.y2)
res.meant.2 <- simulateResiduals(mod.meant2)
plot(res.meant.2)
summary(mod.meant2)
Anova(mod.meant2)


emm1 <- emmeans(mod.meant1, specs = pairwise ~ treatment)
emm1$emmeans
emm1$contrasts

emm2 <- emmeans(mod.meant2, specs = pairwise ~ treatment)
emm2$emmeans
emm2$contrasts


#################### 
# Mean temperature plot (Appendix 1: Fig. A5)

# generate dataframe for mean daily maximum temperature plot
av_meandt <- mean_daily %>% 
  group_by(treatment, period, block) %>% 
  summarize(av_mdt = mean(mean_dt, na.rm = T), se_mdt = std.error(mean_dt),
            sd_mdt = sd(mean_dt)) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "Summer 2019", "Summer 2020"))

# define palettes for colour and shape for plots
pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B", "grey50")
pch.trt <- c(1,1,16,16,16,16,16)

# maximum daily temperature plot
meandt_plot <- ggplot(aes(x = treatment, y = av_mdt, col = treatment, pch = treatment),
                   data = av_meandt) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(width = 0.25, size = 0.8) +
  theme(strip.text = element_text(size = 14)) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  scale_linetype_manual(values = lty.trt) +
  labs(y = "Mean daily maximum temperature (ºC)", 
       x = "Treatment",
       col = "Treatment",
       pch = "Treatment") +
  theme_classic() +
  facet_wrap(~period, scales = "free_x") +
  theme(plot.tag = element_text(face= "bold"))
meandt_plot

# Figure 1 in main text = maximum daily temperature only
png("./figures/FigA5.png", res = 700, width = 7, height = 5, units = "in")
meandt_plot
dev.off()


