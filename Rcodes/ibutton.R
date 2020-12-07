# script to read in and clean ibutton data

pkgs <- c("tidyverse","lubridate")
lapply(pkgs, library, character.only = TRUE)

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
  mutate(date_time = round(strptime(`Date/Time`, format = "%d/%m/%y %I:%M:%S %p"), units = "hours")) %>% 
  select(-Unit, -`Date/Time`) %>% 
  # separate out file information into useful columns
  separate(info, into = c("project","block","number","X", "collect"), sep = "_") %>% 
  select(-project, -X) %>% 
  rename(temp = Value) %>% 
  mutate(number = as.numeric(number),
         # clean up to create a collection date column
         collect_date = str_remove_all(collect, ".csv")) %>% 
  select(-collect) %>% 
  mutate(collect_date = as.Date(collect_date, format = "%Y%m%d"))

setwd("../../")

# this file has info about where tiles were during settlement vs. during the experiment
block_design <- read_csv("./raw_data/design/SVSHW_tilesetup.csv")

# ibuttons have incorrect data for the October 2019 collection after 3 AM, 
# and the January 2020 collection after 11 PM - need to filter these outs

ibutton_filter <- ibutton_separate %>% 
  mutate(remove = ifelse((collect_date == "2019-10-18" & date_time > "2019-10-18 3:00:00")|
                           (collect_date == "2020-03-14" & date_time < "2020-01-22 20:00:00"),
                         TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  select(-remove)

# and now divide data into before tile move and after to 
# properly assign tile numbers to ibuttons for tiles that were renamed

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
  full_join(ibutton_y1)

#write_csv(ibutton_clean, "./clean_data/SVSHW_temp_clean.csv")

#looking at the differences in tile temperatures

temp_summary <- ibutton_clean %>% 
  group_by(current_trt, date_time) %>% 
  summarise(mean_temp = mean(temp, na.rm = TRUE)) %>% 
  rename(Treatment = current_trt) %>% 
  mutate(Treatment = str_replace_all(Treatment, c("black" = "warm", "white"="ambient")))
  
trt_temps <- ggplot(aes(y = mean_temp, x = date_time, colour = Treatment), data = temp_summary) +
  geom_line() +
  facet_wrap(~Treatment) +
  theme_classic() +
  labs(x = "Time", y = "Average hourly temperature (˚C)") +
  scale_colour_manual(values = c("steelblue","grey30","indianred3")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14))
trt_temps

# need to isolate when tiles are out of the water

tides <- read_csv("./raw_data/SVSHW_tides.csv", col_names = FALSE) %>% 
  separate(X1, into = c("date","X","time", "tz", "tide_height"), sep = " ") %>% 
  select(-X, -tz) %>% 
  unite(date_time, c(date, time), sep = " ") %>% 
  mutate(date_time = as.POSIXct(date_time))

ibutton_rocks <- ibutton_clean %>% 
  filter(current_trt == "rock")

above_water <- ibutton_clean %>% 
  full_join(tides) %>% 
  filter(ifelse(collect_date <= "2019-06-06", tide_height < original_shore_level,
                tide_height < new_shore_level)) %>% 
  full_join(ibutton_rocks) %>% 
  mutate(date = strptime(date_time, format = "%Y-%m-%d")) 

# summary stats

max_daily <- above_water %>% 
  group_by(current_trt, date, block, number) %>% 
  summarize(mdt = max(temp)) %>% 
  mutate(date = as.factor(date))

mod.amdt <- lmer(mdt ~ current_trt + (1|date) + (1|block), 
                 data = max_daily)

plot(mod.amdt)
summary(mod.amdt)
Anova(mod.amdt)

average_max_daily <- max_daily %>% 
  group_by(current_trt, date, block) %>% 
  summarize(amdt = mean(mdt), se_amdt = sd(mdt)/sqrt(length(mdt))) %>% 
  mutate(date = as.Date(date)) %>% 
  rename("Treatment" = current_trt) %>% 
  mutate(Treatment = str_replace_all(Treatment, c("white" = "ambient",
                                                  "black" = "warm")),
         Treatment = factor(Treatment, levels = c("rock","ambient","warm")))

amd_plot <- ggplot(aes(x = date, y = amdt, colour = Treatment), data = average_max_daily) +
  geom_line() +
  labs(x = "Date", y = "Average maximum aerial temperature (˚C)") +
  theme_bw() + 
  scale_colour_manual(values = c("grey30","steelblue","indianred3")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  facet_wrap(~Treatment)

amd_plot
