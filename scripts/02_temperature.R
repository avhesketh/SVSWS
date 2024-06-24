# Amelia Hesketh January 2023
########################### 
# Reading in, cleaning, and visualizing iButton temperature data

# Load packages
pkgs <- c("tidyverse","lubridate","lme4", "lmerTest", "car", "stringi", 
          "DHARMa", "patchwork", "glmmTMB", "plotrix", "emmeans",
          "ggnewscale")
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
  mutate(date = date(date_time))

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

# We approximated tile heights based on the tide height when temperatures drop 
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
daily_temps <- ibutton_clean %>% 
  # if date is before treatment swap date, then the "active" treatment is that of y1
  # else it is that of y2
  mutate(treatment = if_else(date_time < swap_date, trt_y1, trt_y2)) %>% 
  mutate(period = if_else(date_time < swap_date, "Year 1", "Year 2")) %>% 
  # only retain those temperatures collected during daytime low tides
  right_join(daylight_tides) %>% 
  left_join(tile_heights) %>% 
  # and only retain temperatures collected when water levels falls below level of the tile
  filter(date > "2019-06-06" & tide_height < tile_level) %>% 
# retain only data collected between 15 June and 31 August of each year
  filter(is.na(temp) == F)  %>% 
  select(-tide_height) %>% 
  filter(date >= "2019-06-01" & date < "2019-09-01" |
           date >= "2020-06-01" & date < "2020-09-01") %>% 
# some loggers, on inspection of data and residuals plots, seem to have malfunctioned. REMOVE these data
  filter(((block == "A" & number == "5" & date >= "2019-07-11" & date <= "2019-10-18") |
         (block == "A" & number =="3" & date >= "2020-06-04" & date <= "2020-08-13") |
         (block == "B" & number == "7" & date >= "2019-07-16" & date <= "2019-08-14")) == F)

# mean daily temperatures
mean_daily <- daily_temps %>% 
  group_by(treatment, date, period, block, number) %>% 
  summarize(mean_dt = mean(temp, na.rm = T)) %>% ungroup()

# ean daily maximum (MDM) temperatures
max_daily <- daily_temps %>% 
  group_by(treatment, date, period, block, number) %>% 
  summarize(mdt = max(temp, na.rm = T)) %>% ungroup() 

# residual MDM temperatures
deviance <- max_daily %>% 
  group_by(date, period) %>% 
  mutate(grand_mean = mean(mdt)) %>% ungroup() %>% 
  group_by(treatment, date, period, block, number) %>% 
  summarize(residual = mdt - grand_mean)

# absolute maximum temperatures
absolute_max <- daily_temps %>% 
  group_by(treatment, period, block, number) %>% 
  summarize(max = max(temp, na.rm = T))


##########################
# Modeling residual MDM temperature


# Model Year 1 and Year 2 temperatures separately due to changes in design from
# 2 (C and W) to 4 (CC, CW, WC, WW) treatments (+ bedrock)

# Filter for the first year of temperature data only and model
trt.y1 <- deviance %>% filter(period == "Year 1") %>% mutate(treatment = factor(treatment))
trt.y1 %>% group_by(treatment) %>% summarize(mean = mean(residual),
                                             se = std.error(residual))

# model as a function of treatment and include random effect for block and tile id
mod.dev1 <- lmer(residual ~ treatment + (1|block/number), data = trt.y1)
plot(mod.dev1)
summary(mod.dev1)
anova(mod.dev1, type = 2)

# Use emmeans for pairwise comparisons of treatments
emm1 <- emmeans(mod.dev1, "treatment")
contrast(emm1, method = "pairwise", adjust = "tukey")

# Filter for the second year of temperature data only and model
trt.y2 <- deviance %>% filter(period == "Year 2") %>% mutate(treatment = factor(treatment))
trt.y2 %>% group_by(treatment) %>% summarize(mean = mean(residual),
                                             se = std.error(residual))
mod.dev2 <- lmer(residual ~ treatment + (1|block/number), data = trt.y2)
plot(mod.dev2)
summary(mod.dev2)
anova(mod.dev2, type = 2)

emm2 <- emmeans(mod.dev2, "treatment")
contrast(emm2, method = "pairwise", adjust = "tukey")

# here, generating labels for eventual plot of MDM temperature based on 
# significant differences by emmeans comparisons
labels.res <- as.data.frame(cbind(
  c("C","W","Rock", "CC","CW","WC","WW","Rock"),
  c("Summer 2019","Summer 2019","Summer 2019",
    "Summer 2020", "Summer 2020", "Summer 2020", 
    "Summer 2020", "Summer 2020"),
  c(1.9,4.5,4.5,1.9,4.3,1.9,4.3,4.3),
  c("a","b","b","c","d","c","d","cd"))
)
colnames(labels.res) <- c("treatment","period", "mean_res","label")
labels.res <- labels.res %>% mutate(mean_res = as.numeric(mean_res))

# examining at overall mean ± SD for treatments and bedrock for reporting
overall_mean <-  daily_temps %>% 
  group_by(treatment, date, period, block, number) %>% 
  summarize(mdt = max(temp, na.rm = T)) %>% ungroup() %>% 
  mutate(treatment = as.factor(treatment)) %>% 
  group_by(treatment, period) %>% 
  summarize(mmdt = mean(mdt), sd_mdt = sd(mdt)) %>% ungroup()


###########################
# Visualizing residual MDM temperatures in each treatment

# generate dataframe for mdm deviance plot
av_deviance <- deviance %>% 
  group_by(treatment, period, block, number) %>% 
  summarize(mean_res = mean(residual, na.rm = T), se_res = std.error(residual),
            sd_res = sd(residual)) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "Summer 2019", "Summer 2020"))

# define palettes for colour and shape for plots
pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B", "grey50")
pch.block <- c(1,2,3,4,5,6)

# maximum daily temperature plot
res_plot <- ggplot(data = av_deviance, aes(x = treatment, y = mean_res, col = treatment)) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(aes(pch = block), height = 0, width = 0.25, alpha = 0.8, size = 0.8,
              show.legend = FALSE) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.block) +
  labs(y = "Residual MDM temperature (ºC)", 
       x = "Treatment",
       col = "Treatment",
       pch = "Block") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~period, scales = "free_x") +
  theme(plot.tag = element_text(face = "bold"),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  geom_text(data = labels.res, aes(label = label), size = 3, col = "black", fontface = "bold")
res_plot

ggsave(res_plot, filename = "./figures/Fig1.pdf", device = cairo_pdf, 
       width = 8.5, height = 9, units = "cm", dpi = 500)


##########################
# ADDITIONAL MODELS AND VISUALIZATIONS OF TEMPERATURE DATA (SEE APPENDIX 2)

# MDM temperatures (not residuals) in each treatment

# Model of maximum daily temperature in y1

# Isolate appropriate data to model (collected in Year 1 of the experiment)
trt.y1 <- max_daily %>% filter(period == "Year 1") %>% mutate(treatment = factor(treatment))

trt.y1 %>% group_by(treatment) %>% summarize(mean = mean(mdt),
                                             sd = sd(mdt))

# Include main effects of treatment
mod.amdt1 <- lmer(mdt ~ treatment + 
                    # nested random effect for individual tile within block
                    # and for date to account for repeated measures
                    (1|block/number) + (1|date),
                  data = trt.y1)

# visualize residuals
plot(mod.amdt1) # violation of normality in residuals

# model summary and anova (Type I)
summary(mod.amdt1)
anova(mod.amdt1, type = 1)

# repeat for model in second year
trt.y2 <- max_daily %>% filter(period == "Year 2") %>% mutate(treatment = factor(treatment))

trt.y2 %>% group_by(treatment) %>% summarize(mean = mean(mdt),
                                          sd = sd(mdt))

# model does not require additional covariates
mod.amdt2 <- lmer(mdt ~ treatment +
                    (1|block/number) + (1|date),
                  data = trt.y2)

plotResiduals(mod.amdt2)
plot(mod.amdt2)
summary(mod.amdt2)
anova(mod.amdt2, type = 1)

# look at significance of contrasts between treatment groups using emmeans=
emm1 <- emmeans(mod.amdt1, "treatment")
contrast(emm1, method = "pairwise", adjust = "tukey")

emm2 <- emmeans(mod.amdt2, "treatment")
contrast(emm2, method = "pairwise", adjust = "tukey")

# Mean temperature modeling (Appendix 2)

trt.y1 <- mean_daily %>% filter(period == "Year 1") %>% mutate(treatment = factor(treatment))

trt.y1 %>% group_by(treatment) %>% summarize(mean = mean(mean_dt),
                                             sd = sd(mean_dt))

mod.meant1 <- lmer(mean_dt ~ treatment +
                        (1|block/number) + (1|date),
                      data = trt.y1)

plot(mod.meant1)
summary(mod.meant1)
anova(mod.meant1, type = 1)

trt.y2 <- mean_daily %>% filter(period == "Year 2") %>% mutate(treatment = factor(treatment),
                                                               julian_day = yday(date))

trt.y2 %>% group_by(treatment) %>% summarize(mean = mean(mean_dt),
                                             sd = sd(mean_dt))

mod.meant2 <- lmer(mean_dt ~ treatment +  
                        (1|block/number) + (1|date), 
                      data = trt.y2)

plot(mod.meant2)
summary(mod.meant2)
anova(mod.meant2, type = 1)
emm_options(opt.digits = FALSE)
emm1 <- emmeans(mod.meant1, "treatment")
contrast(emm1, method = "pairwise", adjust = "tukey")

emm2 <- emmeans(mod.meant2, "treatment")
contrast(emm2, method = "pairwise", adjust = "tukey")

# Absolute maximum temperature modeling

trt.y1 <- absolute_max %>% filter(period == "Year 1") %>% mutate(treatment = factor(treatment))

trt.y1 %>% group_by(treatment) %>% summarize(mean_max = mean(max),
                                             sd_max = sd(max))

mod.amt1 <- lmer(max ~ treatment +
                     (1|block),
                   data = trt.y1)

plot(mod.amt1)
summary(mod.amt1)
anova(mod.amt1, type = 1)

trt.y2 <- absolute_max %>% filter(period == "Year 2") %>% mutate(treatment = factor(treatment))

trt.y2 %>% group_by(treatment) %>% summarize(mean_max = mean(max),
                                                 sd_max = sd(max))

mod.amt2 <- lmer(max ~ treatment +  
                     (1|block), 
                   data = trt.y2)

plot(mod.amt2)
summary(mod.amt2)
anova(mod.amt2, type = 1)


emm1 <- emmeans(mod.amt1, "treatment")
contrast(emm1, method = "pairwise", adjust = "tukey")

emm2 <- emmeans(mod.amt2, "treatment")
contrast(emm2, method = "pairwise", adjust = "tukey")


# Fig A4

ibutton_clean <- read_csv("./clean_data/SVSWS_temp_clean.csv")

# Here, we visualize temperature traces for the whole experiment

# Generate a summary dataframe for plotting with mean and maximum daily temperatures
temp_summary <- ibutton_clean %>% 
  mutate(current_trt = if_else(date_time < "2020-04-03", trt_y1, trt_y2)) %>% 
  group_by(current_trt, date) %>% 
  summarise(mean_temp = mean(temp, na.rm = TRUE), max_temp = max(temp, na.rm = T)) %>% 
  mutate(overlay = factor(if_else(current_trt %in% c("W","CW","WW"),"Warm",
                                  if_else(current_trt %in% c("C","WC","CC"), "Cool", "Rock")),
                          levels = c("Cool", "Warm", "Rock")),
         current_trt = factor(current_trt, levels = c("C","CC","CW","W","WW","WC","Rock")))

# define color palette
pal.traces <- c("#014779", "#014779","#7985CB", "#EE4B2B","#EE4B2B","#9C0098","grey50")

# read in survey times to show those on the plot
surveys <- read_csv("./raw_data/design/SVSWS_survey_times.csv") %>% 
  group_by(survey_no) %>% 
  summarize(date = max(date)) %>% 
  # this variable is just showing post-summer and winter timepoints when data were
  # used in biological analyses 
  mutate(important = if_else(survey_no %in% c(8,10,14,16), "y","n"))

# First panel: Average maximum temperature
FigS4a <- ggplot(aes(y = max_temp, x = date, colour = current_trt), 
                        data = temp_summary) +
  geom_line(lwd = 0.5, alpha = 0.7) +
  facet_grid(rows = vars(overlay)) +
  theme_classic() +
  scale_color_manual(values = pal.traces) +
  labs(x = "Time", y = "Average maximum temperature (˚C)", colour = "Treatment") +
  new_scale_colour() +
  geom_vline(aes(xintercept = date, color = important), lwd = 0.4, lty = "dashed", data = surveys, show.legend = F) +
  scale_color_manual(values = c("grey80","grey20")) +
  theme(panel.grid.major.y = element_line(linewidth= 0.3),
        panel.grid.minor.y = element_line(linewidth = 0.2),
        plot.tag = element_text(face = "bold")) +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS4a

# Second panel: Mean daily temperature
FigS4b <- ggplot(aes(y = mean_temp, x = date, colour = current_trt), 
                         data = temp_summary) +
  geom_line(lwd = 0.5, alpha = 0.7) +
  facet_grid(rows = vars(overlay)) +
  theme_classic() +
  scale_color_manual(values = pal.traces) +
  labs(x = "Time", y = "Average mean temperature (˚C)", colour = "Treatment") +
  new_scale_colour() +
  geom_vline(aes(xintercept = date, color = important), lwd = 0.4, lty = "dashed", data = surveys, show.legend = F) +
  scale_color_manual(values = c("grey70","grey20")) +
  theme(panel.grid.major.y = element_line(linewidth= 0.3),
        panel.grid.minor.y = element_line(linewidth= 0.2),
        plot.tag = element_text(face = "bold")) +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS4b

# A3C: Example trace showing treatment effect over one summer low tide series
# Isolate one fragment of temperature data
sample_data <- ibutton_clean %>% 
  filter(date %in% c("2019-07-28","2019-07-29",
                            "2019-07-30","2019-07-31","2019-08-01", "2019-08-02", "2019-08-03")) %>% 
  mutate(current_trt = if_else(date_time < "2020-04-03 16:00:00", 
                               trt_y1, if_else(trt_y2 == "Rock", trt_y2,substr(trt_y2,2,2)))) %>% 
  group_by(current_trt, date_time) %>% 
  summarise(mean_temp = mean(temp, na.rm = TRUE), mean_shore_level = mean(new_shore_level))

# Isolate analogous chunk of tide height data
tides_trace <- tides %>% 
  mutate(date = date(date_time)) %>% 
  filter(date %in%  c("2019-07-28","2019-07-29",
                             "2019-07-30","2019-07-31","2019-08-01", "2019-08-02", "2019-08-03"))

# Join together the data, indicating when the tide is above the tiles and vice versa
sample_data2 <- sample_data %>% left_join(tides_trace) %>% 
  mutate(above_below = if_else(tide_height >= mean_shore_level, "below", "above")) %>% 
  mutate(current_trt = factor(current_trt, levels = c("C","W","Rock")))

# Third panel
FigS4c <- ggplot(aes(x = date_time), data = sample_data2) +
  geom_tile(aes(y = 20, height = 40, fill = above_below), alpha = 0.2) +
  geom_line(aes(y = mean_temp, colour = current_trt), lwd = 0.5, show.legend = F) +
  theme_classic() +
  labs(x = "Date", y = "Average hourly temperature (˚C)", colour = "Treatment") +
  scale_colour_manual(values = c("#014779","#EE4B2B","grey50")) +
  geom_line(aes(y = tide_height*10), lwd = 0.5, lty = "dashed", col = "darkblue") +
  scale_y_continuous(sec.axis = sec_axis(~.*0.1, name="Tide height (m)\n")) +
  scale_fill_manual(values = c("lemonchiffon","lightblue")) +
  guides(fill = "none") +
  theme(plot.tag = element_text(face = "bold"))
FigS4c

# Assemble multipanel figure using patchwork
FigS4 <- (FigS4a / FigS4b / FigS4c) + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
FigS4

# Save file
png("./figures/FigS4.png", width = 10, height = 10, units = "in",
    res = 700)
FigS4
dev.off()

# Note that some data are missing for rock temperatures between 19 July to 19 August 2020.
# iButtons were pre-programmed incorrectly, and started recording too early.

# Max & mean temperature plots (Appendix 1: Fig. A5)

# generate dataframe for mean daily maximum temperature plot
av_mdt <- max_daily %>% 
  group_by(treatment, period, block, number) %>% 
  summarize(mean_mdt = mean(mdt, na.rm = T), se_mdt = std.error(mdt),
            sd_mdt = sd(mdt)) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "Summer 2019", "Summer 2020"))

# generate dataframe for mean daily temperature plot
av_meandt <- mean_daily %>% 
  group_by(treatment, period, block, number) %>% 
  summarize(av_mdt = mean(mean_dt, na.rm = T), se_mdt = std.error(mean_dt),
            sd_mdt = sd(mean_dt)) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "Summer 2019", "Summer 2020"))

# generate dataframe for absolute maximum temperature plot
max_plot <- absolute_max %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "Summer 2019", "Summer 2020"))

# define palettes for colour and shape for plots
pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B", "grey50")
pch.block <- c(1,2,3,4,5,6)

# Maximum daily temperature (not residuals) plot


# Raw (not residual) mean daily maximum temperature plot
FigS5a <- ggplot(data = av_mdt, aes(x = treatment, y = mean_mdt, col = treatment)) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(aes(pch = block), height = 0, width = 0.25, alpha = 0.8, size = 0.8,
              show.legend = FALSE) +
  theme(strip.text = element_text(size = 14)) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.block) +
  labs(y = "Mean daily max temperature (ºC)", 
       x = "Treatment",
       col = "Treatment",
       pch = "Block") +
  theme_classic() +
  facet_wrap(~period, scales = "free_x") +
  theme(plot.tag = element_text(face= "bold"))
FigS5a

# mean daily tempearture plot
FigS5b <- ggplot(aes(x = treatment, y = av_mdt, col = treatment),
                   data = av_meandt) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(aes(pch = block), alpha = 0.8, height = 0, size = 0.8, show.legend = F) +
  theme(strip.text = element_text(size = 14)) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.block) +
  labs(y = "Mean temperature (ºC)", 
       x = "Treatment",
       col = "Treatment",
       pch = "Treatment") +
  theme_classic() +
  facet_wrap(~period, scales = "free_x") +
  theme(plot.tag = element_text(face= "bold"))
FigS5b

# absolute maximum temperature
FigS5c <- ggplot(aes(x = treatment, y = max, col = treatment),
                 data = max_plot) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(aes(pch = block), alpha = 0.8, height = 0, size = 0.8, show.legend = F) +
  theme(strip.text = element_text(size = 14)) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.block) +
  labs(y = "Maximum temperature (ºC)", 
       x = "Treatment",
       col = "Treatment",
       pch = "Treatment") +
  theme_classic() +
  facet_wrap(~period, scales = "free_x") +
  theme(plot.tag = element_text(face= "bold"))
FigS5c

# assemble multipanel plot
FigS5 <- (FigS5a / FigS5b/ FigS5c) + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect")

# Save figure
png("./figures/FigS5.png", res = 700, width = 6, height = 9, units = "in")
FigS5
dev.off()


