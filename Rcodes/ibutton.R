### May 2019 
# script to read in ibutton data

library(tidyverse)
library(lubridate)

# now figure out where the files are and their names

setwd("./ibuttons")

file_names <- list.files()

ibuttons <- data.frame()

for (i in 1:length(file_names)){
  temp_data <- read_csv(file_names[i], skip = 14)
  temp_data$info <- file_names[i]
  ibuttons <- rbind(ibuttons, temp_data)
  }

## split date and time column

ibuttons$`Date/Time` <- as.character(ibuttons$`Date/Time`)

ibutton_dates <- ibuttons %>%
  separate("Date/Time", into = c("date", "time", "am"), sep = " ") %>% 
  unite(time_12hr, time, am, sep = " ")

# get dates and times in useful formats and then rejoin

ibutton_dates$date <- as.Date(ibutton_dates$date, format = "%d/%m/%y")

ibutton_dates$time_24hr <- format(strptime(ibutton_dates$time_12hr, "%I:%M:%S %p"),  format="%H:%M:%S")

# extract block and tile info from file name
ibutton_data <- ibutton_dates %>% 
  select(-2, -3) %>% 
  unite(date_time, date, time_24hr, sep = " ") %>% 
  separate(info, into = c("tileno","timestamp"), sep = "_") 

# and now split up the block and tile number
ibutton_data <- ibutton_data %>% 
  mutate(block = substr(ibutton_data$tileno, 1, 1),
         number = substr(ibutton_data$tileno, 2, length(ibutton_data$tileno))) %>% 
  rename(temp = Value)

ibutton_data <- ibutton_data %>% 
  dplyr::select(-tileno, -timestamp)

# convert the double-digit integers (e.g. 01) to regular numbers
# for rock temperature, preserve that identifier
for (i in 1:length(ibutton_data$number)){
  if (is.na(ibutton_data$number[i])) {
    ibutton_data$number[i] = "rock"
  }
}

# now add in the colour of the tile

# this file has info about when tiles got moved and tile colours
block_design <- read_csv("../moving_guide.csv")

# remove the chunk of time when ibuttons were all being replaced
ibutton_data_remove <- ibutton_data %>% 
  filter(date_time >= "2019-06-05 12:00:00" &  
           date_time <= "2019-06-06 18:00:00")

ibutton_useful <- ibutton_data %>% 
  anti_join(ibutton_data_remove)

# separate out rock and tile temperature records
ibutton_tiles <- ibutton_useful %>% 
  filter(number != "rock" & number != "Rock")

ibutton_rocks <- ibutton_useful %>%
  filter(number == "rock" | number == "Rock") %>% 
  mutate(colour = "rock", number = "rock")

write.csv(ibutton_tiles, "ibuttons_raw.csv")
write.csv(ibutton_rocks, "rock_ibuttons_raw.csv")

# and now divide data into before tile move and after to 
# properly assign tile numbers to ibuttons for tiles that were renamed
ibutton_premove <- ibutton_tiles %>% 
  filter(date_time < "2019-06-05 12:00:00")

ibutton_postmove <- ibutton_tiles %>% 
  filter(date_time > "2019-06-06 18:00:00")

# need to get block and tile number into one column to allow comparison
# to block design guide
ibutton_joint <- ibutton_premove %>% 
  unite(b_no, c(block, number), sep = "_", remove = TRUE)

joint_legend <- block_design %>% 
  unite(b_no_original, c(original_block, original_no), sep = "_", remove = TRUE) %>% 
  unite(b_no_new, c(new_block, new_no), sep = "_", remove = TRUE)

# now compare each entry against block design guide. if a tile
# was moved, then replace the old tile number/block with the 
# one after the move
for (i in 1:length(ibutton_joint$b_no)){
  for (j in 1:length(joint_legend$b_no_original)){
    if (ibutton_joint$b_no[i] == joint_legend$b_no_original[j]){
      ibutton_joint$b_no[i] = joint_legend$b_no_new[j]
    }
  }
}

# and now get columns separate again
ibutton_pre <- ibutton_joint %>% 
  separate(b_no, into = c("block", "number"), sep = "_", remove = TRUE)

block_design_new <- block_design %>% 
  select(-original_block, -original_no) %>% 
  rename(block = new_block, number = new_no)

ibutton_pre$number <- as.numeric(ibutton_pre$number)
ibutton_pre <- na.omit(ibutton_pre)

ibutton_premove <- ibutton_pre %>% 
  full_join(block_design_new)

ibutton_postmove$number <- as.numeric(ibutton_postmove$number)

# join data together!
ibutton_postmove <- ibutton_postmove %>% 
  full_join(block_design_new)

ibuttons_all <- ibutton_premove %>% 
  rbind(ibutton_postmove) %>% 
  rbind(ibutton_rocks) %>% 
  na.omit()
write.csv(ibuttons_all, "../clean_data.csv")

# now let's plot it
library(tidyverse)
cd_ibutton <- read_csv("clean_data.csv")

glimpse(cd_ibutton)

rocks <- cd_ibutton %>% 
  filter(colour == "rock") %>% 
  select(-1, -5, -6)
glimpse(rocks)

write_csv(rocks, "rock_temp_summer19.csv")

# create summary dataframe of maximum daily temperature
# and mean daily temperature for each treatment over time

# mean daily
ibuttons_mean <- ibuttons_all %>% 
  separate(date_time, c("date", "time"), sep = " ") %>% 
  group_by(colour, date, block) %>% 
  summarize(mean_temp = mean(temp), se_temp = sd(temp)/sqrt(length(temp)))

# max daily
ibuttons_max_tiles <- ibuttons_all %>% 
  separate(date_time, c("date", "time"), sep = " ") %>% 
  group_by(colour, date, block, number) %>% 
  summarize(max_temp = max(temp))

ibuttons_max_trts <- ibuttons_max_tiles %>% 
  group_by(colour, date, block) %>% 
  summarize(mdm = mean(max_temp), se_mdm = sd(max_temp)/sqrt(length(max_temp)))

# plots
# first, mean

ibuttons_mean$date <- as.Date(ibuttons_mean$date)

temp_time <- ggplot(aes(x = date, y = mean_temp, colour = colour), data = ibuttons_mean) +
  geom_point(size = 0.5) +
  geom_line()
temp_time

# calculating means - strangely, black tiles = white = rock
# could it be that it's only relevant for low tides? if so, can I isolate
# low tide temps?

grand_mean <- ibuttons_mean %>% 
  group_by(colour) %>% 
  summarize(overall_mean = mean(mean_temp))

grand_mdm <- ibuttons_max_trts %>% 
  group_by(colour) %>% 
  summarize(overall_mean = mean(mdm))

# ibuttons summer data only

ibuttons_summer <- ibuttons_all %>% 
  separate(date_time, c("date", "time"), sep = " ") %>% 
  filter(date < "2019-09-01" & date > "2019-06-01")

mean_summer <- ibuttons_summer %>% 
  group_by(colour) %>% 
  summarize(overall_mean = mean(temp))

max_summer <- ibuttons_summer %>% 
  group_by(colour, block, number, date) %>% 
  summarize(max_temp = max(temp)) %>% 
  group_by(colour) %>% 
  summarize(mdm = mean(max_temp))

# ok, so the difference does hold during the summer, but black tiles apparently
# more similar to rock temp?

mean_summ <- ibuttons_summer %>% 
  group_by(colour, date, block) %>% 
  summarize(mean_temp = mean(temp))

ibuttons_summer$date <- as.Date(ibuttons_summer$date)

temp_time <- ggplot(aes(x = date, y = mean_temp, colour = colour), data = mean_summ) +
  geom_point(size = 0.5) +
  geom_line()
temp_time

mdm_summ <- ibuttons_summer %>% 
  group_by(colour, date, block, number) %>% 
  summarize(max_daily = max(temp))

mdm_summ <- mdm_summ %>% 
  group_by(colour, date, block) %>% 
  summarize(mdm = mean(max_daily),
            se_mdm = sd(max_daily)/sqrt(length(mdm_summ$max_daily))) #%>% 
 # filter(colour != "rock")

mdm_time <- ggplot(aes(x=date,y=mdm,colour=colour), data = mdm_summ)+
  geom_line(se = TRUE)
mdm_time



