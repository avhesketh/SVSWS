## exploring tide data
library(tidyverse)
library(lubridate)

setwd("~/Documents/UBC/R_projects/natgeo")

# read in data from DFO hourly tides at Fulford

tides <- read_csv("tide_data_fulford.csv", col_names = FALSE) %>% 
  separate(X1, c("date","blank", "time", "time_zone", "tide_level"), sep = "([\\  \\ ])") %>% 
  select(-blank) #%>% 
  #unite(date_time, c(date, time), sep = " ") 

tides$tide_level <- as.numeric(tides$tide_level)
tides$date <- as.Date(tides$date)

# before plotting convert to POSIXct format to preserve date and time info as a continuous variable
#tides$date_time <- as.POSIXct(tides$date_time)

tide_plot <- ggplot(aes(x = date_time, y = tide_level), data = tides) +
  geom_line()
tide_plot

# try to filter all the data for individuals blocks for when tiles are actually out of water

# break up tide level data into before and after move chunks
tides_pre <- tides %>% 
  filter(date < "2019-06-05") %>% 
  unite(date_time, c(date, time), sep = " ")
tides_post <- tides %>% 
  filter(date >= "2019-06-06") %>% 
  unite(date_time, c(date, time), sep = " ")

# create list of dataframes with irrelevant tide data filtered out
tile_height <- read_csv("tile_shorelevel.csv")

tiles_aw <- data.frame(matrix(ncol = 4, nrow = 0))
colnames <- c("date_time", "time_zone", "tide_level", "location")
colnames(tiles_aw) <- colnames 
tiles_aw$date_time <- as.character(tiles_aw$date_time)
tiles_aw$time_zone <- as.character(tiles_aw$time_zone)
tiles_aw$tide_level <- as.numeric(tiles_aw$tide_level)
tiles_aw$location <- as.character(tiles_aw$location)

for (i in 1:length(tile_height$block)) {
  if (tile_height$old_new[i] == "new") {
    AW <- tides_post %>% 
      filter(tide_level < tile_height$tide_height[i]) %>% 
      mutate(location = paste(tile_height$block[i], tile_height$old_new[i], sep = "_"))
  }
  if (tile_height$old_new[i] == "old") {
    AW <- tides_pre %>% 
      filter(tide_level < tile_height$tide_height[i]) %>% 
      mutate(location = paste(tile_height$block[i], tile_height$old_new[i], sep = "_"))
  }
  tiles_aw <- tiles_aw %>% 
    full_join(AW)
}

# now apply to the ibutton data and get temperature only at times where tiles are out of water

ibuttons_tiles <- read_csv("ibuttons_raw.csv") %>% 
  select(-X1)
ibuttons_rocks <- read_csv("rock_ibuttons_raw.csv") %>% 
  select(-X1)

ibuttons_all <- ibuttons_tiles %>% 
  full_join(ibuttons_rocks)

ibuttons_pre <- ibuttons_all %>% 
  filter(date_time < "2019-06-05 12:00:00") %>% 
  mutate(old_new = "old")
ibuttons_all <- ibuttons_all %>% 
  filter(date_time > "2019-06-06 18:00:00") %>% 
  mutate(old_new = "new") %>% 
  full_join(ibuttons_pre) %>% 
  mutate(location = paste(block, old_new, sep = "_"))

# get everything into a valid date format
tiles_aw$date_time <- as.POSIXlt(tiles_aw$date_time, format = "%Y-%m-%d %H:%M")
ibuttons_all$date_time <- as.POSIXct(ibuttons_all$date_time, format = "%Y-%m-%d %H:%M")

# easiest to just round to the nearest hour and then join df
tiles_aw$date_time <- format(round(tiles_aw$date_time, units = "hours"), format = "%Y-%m-%d %H:%M")
ibuttons_all$date_time <- format(round(ibuttons_all$date_time, units = "hours"), format = "%Y-%m-%d %H:%M")

low_tide <- ibuttons_all %>% 
  inner_join(tiles_aw) %>% 
  select(-old_new, -location, -colour, -time_zone)

# now have a df of all ibuttons measurements from when the tide was low that I can feed into previous code
######

ibutton_lt_pre <- low_tide %>% 
  filter(date_time < "2019-06-05 12:00:00")

ibutton_lt_post <- low_tide %>% 
  filter(date_time > "2019-06-06 18:00:00")

# need to get block and tile number into one column to allow comparison
# to block design guide
ibutton_joint <- ibutton_lt_pre %>%
  unite("b_no", c(block, number), sep = "_", remove = FALSE)

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

ibutton_rock <- ibutton_joint %>% 
  separate(b_no, into = c("block", "number"), sep = "_", remove = TRUE) %>% 
  filter(number == "rock") %>% 
  unique()

ibutton_rock_full <- ibutton_lt_post %>% 
  filter(number == "rock") %>%
  unique() %>% 
  full_join(ibutton_rock) %>% 
  mutate(colour = "rock")

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
ibutton_lt_post$number <- as.numeric(ibutton_lt_post$number)
ibutton_lt_post <- na.omit(ibutton_lt_post)

ibutton_postmove <- ibutton_lt_post %>% 
  full_join(block_design_new)

ibuttons_all_lowtide <- ibutton_premove %>% 
  rbind(ibutton_postmove) %>% 
  rbind(ibutton_rock_full) %>% 
  na.omit() %>% 
  select(-tide_level)

###
ibuttons_max_lt <- ibuttons_all_lowtide %>% 
  separate(date_time, c("date", "time"), sep = " ") %>% 
  group_by(colour, date, block, number) %>% 
  summarize(max_temp = max(temp))

ibuttons_max_lt_trts <- ibuttons_max_lt %>% 
  group_by(colour, date, block) %>% 
  summarize(mdm = mean(max_temp), se_mdm = sd(max_temp)/sqrt(length(max_temp))) %>% 
  filter(colour != "rock")

# for only times when tiles are out of water, black = rock max daily temp, white = lower ... hmmmmm
ibuttons_max_lt_trts %>% group_by(colour) %>% summarize(mdm = mean(mdm))

ibuttons_max_lt_trts$date <- as.Date(ibuttons_max_lt_trts$date)
lowtide_mdt <- ggplot(aes(x = date, y = mdm, colour = colour), data = ibuttons_max_lt_trts) +
  geom_line() + 
  facet_wrap(~block)

lowtide_mdt
