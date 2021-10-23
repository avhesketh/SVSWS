# script to read in and clean ibutton data

pkgs <- c("tidyverse","lubridate","lmer", "car")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# now figure out where the files are and their names

setwd("./raw_data/temperature")

file_names <- list.files()

ibuttons <- data.frame()

for (i in 1:length(file_names)){
  temp_data <- read_csv(file_names[i], skip = 14)
  temp_data$info <- file_names[i]
  ibuttons <- rbind(ibuttons, temp_data)
  }

## split date and time column

ibuttons$`Date/Time` <- as.character(ibuttons$`Date/Time`)

ibutton_separate <- ibuttons %>%
  mutate(date_time = dmy_hms(`Date/Time`)) %>% 
  select(-Unit, -`Date/Time`) %>% 
  # separate out file information into useful columns
  separate(info, into = c("project","block","number","X", "collect"), sep = "_") %>% 
  select(-project, -X) %>% 
  rename(temp = Value) %>% 
  mutate(number = as.numeric(number),
         # clean up to create a collection date column
         collect_date = str_remove_all(collect, ".csv")) %>% 
  select(-collect) %>% 
  mutate(collect_date = ymd(collect_date)) %>% 
  separate(date_time, into = c("record_date","record_time"), sep = " ", remove = F)
setwd("../../")

# this file has info about where tiles were during settlement vs. during the experiment
block_design <- read_csv("./raw_data/design/SVSHW_tilesetup.csv")

# ibuttons were sometimes set up or data were collected before/after field work 
# due to laptop not being around or time/personnel constraints
# need to filter these false data out

ibutton_filter <- ibutton_separate %>% 
  mutate(remove = ifelse((collect_date == "2019-10-18" & date_time > "2019-10-18 3:00:00")|
                           (collect_date == "2020-03-14" & date_time < "2020-01-22 20:00:00")|
                           (collect_date == "2019-09-14" & date_time > "2019-09-14 12:00:00")|
                           (date_time >= "2021-02-24 19:00:00"),
                         TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  select(-remove)

# and now divide data into before tile move and after
# to properly assign tile numbers to ibuttons for tiles that were renamed

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
  mutate(remove = ifelse(is.na(date_time), TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  select(-remove) %>% 
  unite(trt, c(colour_y1, colour_y2), sep = "_", remove = FALSE) %>% 
  mutate(trt = ifelse(number == 0, "rock", trt))
  
ibutton_y1 <- ibutton_position %>% 
  filter(date_time < "2020-04-03 17:00:00") %>% 
  mutate(current_trt = ifelse(trt == "rock", trt, colour_y1))

ibutton_clean <- ibutton_position %>% 
  filter(date_time >= "2020-04-03 17:00:00") %>% 
  mutate(current_trt = ifelse(trt == "rock", trt, colour_y2)) %>% 
  full_join(ibutton_y1) %>% 
  mutate(hour = hour(date_time)) %>% 
  mutate(date_time = ymd_h(paste(record_date, hour, sep = " ")))

cols_to_remove <- c("record_date", "record_time", "original_angle","new_angle",
                    "new_angle","new_compass","hour","trt", "colour_y1",
                    "colour_y2", "survived","original_herb_trt")

ibutton_clean <- ibutton_clean %>% select(-cols_to_remove) %>% filter(is.na(temp) == F)

#write_csv(ibutton_clean, "./clean_data/SVSHW_temp_clean.csv")

#looking at the differences in tile temperatures

temp_summary <- ibutton_clean %>% 
  group_by(current_trt, date_time) %>% 
  summarise(mean_temp = mean(temp, na.rm = TRUE)) %>% 
  mutate(current_trt = str_replace_all(current_trt, c("black" = "warm", "white"="cool")))


trt_temps <- ggplot(aes(y = mean_temp, x = date_time, colour = current_trt), 
                    data = temp_summary %>% filter(current_trt %in% c("warm", "cool"))) +
  geom_line() +
  facet_wrap(~current_trt) +
  theme_bw() +
  labs(x = "Time", y = "Average hourly temperature (˚C)", colour = "Treatment") +
  scale_colour_manual(values = c("cadetblue","firebrick4")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ylim(c(-5,40))

View(temp_summary)
# need to isolate when tiles are out of the water

tides <- read_csv("./raw_data/design/SVSHW_tides.csv", col_names = FALSE) %>% 
  separate(X1, c("date","blank", "time", "time_zone", "tide_height"), sep = "([\\  \\ ])") %>% 
  select(-blank) %>% 
  unite(c("date", "time"), col = date_time, sep = " ") %>% 
  mutate(tide_height = as.numeric(tide_height), date_time = ymd_hm(date_time))

ibutton_rocks <- ibutton_clean %>% 
  filter(current_trt == "rock")

above_water <- ibutton_clean %>% 
  left_join(tides) %>% 
  filter(ifelse(collect_date <= "2019-06-06", tide_height < original_shore_level,
                tide_height < new_shore_level)) %>% 
  separate(date_time, into = c("date","time"), sep = " ")

# summary stats

max_daily <- above_water %>% 
  group_by(current_trt, date, block, number) %>% 
  summarize(mdt = max(temp, na.rm = T)) %>% 
  mutate(date = as.factor(date))

mod.amdt <- lmer(mdt ~ current_trt + (1|block) + (1|date), 
                 data = max_daily)

plot(mod.amdt)
summary(mod.amdt)
Anova(mod.amdt)

average_max_daily <- max_daily %>% 
  group_by(current_trt, date) %>% 
  summarize(amdt = mean(mdt, na.rm=T), se_amdt = sd(mdt,na.rm = T)/sqrt(length(mdt))) %>% 
  #mutate(date = as.Date(date)) %>% 
  mutate(current_trt = str_replace_all(current_trt, c("black" = "warm", "white"="cool")),
         current_trt = factor(current_trt, levels = c("rock","cool","warm")),
         date = as.Date(date))

amd_plot <- ggplot(aes(x = date, y = amdt, colour = current_trt), 
                   data = average_max_daily %>% filter(current_trt %in% c("warm", "cool"))) +
  geom_line() +
  labs(x = "Date", y = "Average maximum aerial temperature (˚C)", colour = "Treatment") +
  theme_bw() + 
  scale_colour_manual(values = c("cadetblue","firebrick4")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  facet_wrap(~current_trt) +
  geom_rect(aes(xmin = as.Date("2019-10-15"), xmax = as.Date("2020-03-15"), ymin = 0, ymax = 45), colour = "grey90", fill = "grey90", alpha = 0.02) +
  geom_rect(aes(xmin = as.Date("2020-10-15"), xmax = as.Date("2021-02-24"), ymin = 0, ymax = 45), colour = "grey90", fill = "grey90", alpha = 0.02) +
  ylim(c(0,45))

amd_plot

?geom_
amd_plot <- ggplot(aes(x = Treatment, y = amdt, colour = Treatment), data = average_max_daily) +
  geom_point() +
  geom_errorbar(aes(ymax = amdt + se_amdt, ymin = amdt - se_amdt)) +
  theme_bw()
amd_plot

# only when tides are daytime...

ibutton_day <- above_water %>% 
  filter((date < "2019-09-15" | date >= "2020-04-15" & date <= "2020-09-15") == TRUE) %>% 
  group_by(current_trt, date, block, number) %>% 
  summarize(mdt = max(temp, na.rm = T))

av_mdt <- ibutton_day %>% group_by(current_trt) %>% summarize(av_mdt = mean(mdt),
                                                                             se_mdt = sd(mdt)/sqrt(length(mdt)))

mdt_plot <- ggplot(aes(x = current_trt, y = av_mdt, col = current_trt),
                   data = av_mdt) +
  geom_point() +
  geom_errorbar(aes(ymax = av_mdt + se_mdt, ymin = av_mdt - se_mdt)) +
  theme_bw()
mdt_plot

model.temp <- lmer(mdt ~ current_trt + (1|block) + 
                    (1|date), data = only_2020)
summary(model.temp)
Anova(model.temp)


lineplot_mdt <- ggplot(aes(x = as.Date(date), color = current_trt, y = mdt), data = ibutton_day) +
  geom_line() + 
  facet_wrap(~current_trt)
lineplot_mdt

only_2020 <- ibutton_day %>% filter(date > "2020-01-01")


