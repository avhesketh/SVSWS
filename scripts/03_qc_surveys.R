# Cleaning raw biological data
# Amelia Hesketh January 2022

###################
# Cleaning tile survey data########

# load necessary packages
library(assertr)
library(tidyverse)

# read in raw data
survey_data <- read_csv("./raw_data/tile_surveys/SVSWS_surveys.csv",
                        col_types = c("D","f","f","c","n","n","c","?")) %>% 
  mutate(species = str_replace_all(species, c("Lottia_sp" = "Lottia_sp_recruits", 
                                              "scunge" = "Scunge"))) %>% 
  select(-8)

# isolate algae and invertebrates separately

sp_list <- survey_data$species %>% unique %>% 
  as.data.frame() %>% rename(species = 1)

alga_list <- as.data.frame(c("Ulothrix_sp", "Ulva_sp", "Savoiea_robusta", "Pyropia_sp",
           "Leathesia_marina", "Fucus_distichus","Endocladia_muricata",
           "Petalonia_fascia", "Scunge","Mastocarpus_sp_crust")) %>% 
  rename(species = 1)

invert_list <- sp_list %>% anti_join(alga_list)

spp_list <- alga_list %>% mutate(type = "alga") %>% 
  rbind(invert_list %>% mutate(type = "invert")) %>% 
  write_csv("./clean_data/SVSWS_species_list.csv")

# check that all algae have percent cover recorded and that all invertebrates have count data
survey_data %>% filter(species %in% invert_list$species) %>% verify(is.na(percent_cover)) %>% verify(count >= 0)

# remove all NA count data
survey_data_qc <- survey_data %>% mutate(remove = if_else(species %in% invert_list$species & is.na(count),
                                                         T, F)) %>% filter(remove == F) %>% select(-remove)
# verify that all remaining count data are > 0
survey_data_qc %>% filter(species %in% invert_list$species) %>% verify(count >= 0)

survey_data_alga <- survey_data_qc %>% filter(species %in% alga_list$species) 

survey_data_alga %>% verify(percent_cover >= 0)

# cleaning the data. put percent cover and count data into separate columns
survey_data_qc2 <- survey_data_qc %>% mutate(percent_cover = if_else(species %in% alga_list$species & (is.na(count) == F), 
                                                                    count, percent_cover)) %>%
  # make all count column entries NA for algal species
  mutate(count = case_when(species %in% alga_list$species & (is.na(count) == F) ~ NA_real_,
                            TRUE ~ count)) %>% 
  # remove all empty percent cover data
  mutate(remove = if_else(species %in% alga_list$species & is.na(percent_cover), T, F)) %>% 
  filter(remove == F) %>% select(-remove)

survey_data_alga <- survey_data_qc2 %>% filter(species %in% alga_list$species) 

survey_data_alga %>% verify(percent_cover >= 0)
# now it all looks good

## Need to also account for the fact that tiles were moved around 
## in June 2019 due to changes in the experimental design.

# links overall unique tile ID to block and number within block
block_design <- read_csv("./raw_data/design/SVSWS_tilesetup.csv")

# filter out surveys after tiles were moved around
# in this case, the block and number recorded are the "new" block and number
survey_postmove <- survey_data_qc2 %>% 
  filter(date > "2019-06-06") %>% 
  # relabel as new number/block
  rename("new_no" = number,
         "new_block" = block) %>%
  # this will map to the tile design file
  left_join(block_design) %>% 
  rename(angle = new_angle,
         shore_level = new_shore_level,
         block = new_block) %>% 
  select(-new_no, -survived, -original_no, -original_shore_level,
         -original_angle, -colour_y1, -colour_y2, -original_block)

# filter out the surveys before tiles were moved around
# in this case, the block and number are the "old" block and number
survey_premove <- survey_data_qc2 %>% 
  filter(date <= "2019-06-06") %>% 
  rename("original_no" = number,
         "original_block" = block) %>% 
  left_join(block_design) %>% 
  rename(angle = original_angle,
         shore_level = original_shore_level,
         block = original_block) %>% 
  select(-original_no, -survived, -new_no, -new_shore_level,
         -new_angle, -colour_y1, -colour_y2,  -new_block)

# join all these data together such that individual tiles are now tracked by their ID,
# not by block/number combo since that changed through time.
survey_all <- survey_premove %>% full_join(survey_postmove)

# decompose treatment into overall treatment in year one and year two treatments
surveys_clean <- survey_all %>%
  separate(treatment, into = c("trt_y1", "trt_y2"), sep = 1, remove = FALSE)

# also need a vector associating sampling dates with a survey id number
survey_list <- surveys_clean %>% select(date) %>% unique() %>% arrange(date)

# note that some surveys occurred over multiple dates if tides were not favourable
survey_no <- c(1,1,2,2,3,4,4,5,5,6,6,7,8,8,9,10,11,12,13,14,15,16)
survey_join <- cbind(survey_list, survey_no)

# integrate survey id into cleaned data
surveys_clean2 <- surveys_clean %>% left_join(survey_join)

tile_info <- surveys_clean2 %>% select(tile_id, survey_no, block, 
                                       first_herb_trt,
                                      second_herb_trt, angle, 
                                       new_compass,shore_level,treatment,
                                       trt_y1, trt_y2, notes) %>% 
  unique() %>% write_csv("./clean_data/SVSWS_tile_treatments.csv")

tile_data <- surveys_clean2 %>% select(tile_id, survey_no, species, count, percent_cover)

# fill in the implicit zeroes with EXPLICIT zero data

combos <- tile_data %>% 
  select(tile_id, survey_no) %>% unique()

# for the count data column
count_data <- tile_data %>% 
  select(-percent_cover) %>% 
  na.omit %>% 
  right_join(combos) %>% 
  complete(nesting(tile_id, survey_no), species, fill = list(count = 0)) %>% 
  na.omit()

# for the cover data column
cover_data <- tile_data %>% 
  select(-count) %>% 
  na.omit() %>% 
  right_join(combos) %>% 
  complete(nesting(tile_id, survey_no), species, fill = list(percent_cover = 0)) %>% 
  na.omit()

# re-join these data
all_data <- count_data %>% full_join(cover_data) %>% left_join(tile_info)

#  save the clean data
write_csv(all_data, "./clean_data/SVSWS_survey_clean.csv")

###################
# Cleaning epifaunal data (destructive sampling)

# load in the two dataframes in question
sept20_samples <- read_csv("./raw_data/epifauna/SVSWS_20200914_epifauna.csv") %>% mutate(date = "2020-09-14")
feb21_samples <- read_csv("./raw_data/epifauna/SVSWS_20210224_epifauna.csv") %>% mutate(date = "2021-02-24")

# need to clean up the species names to be consistently formatted

epifauna_original <- as.data.frame(unique(c(sept20_samples$taxon, feb21_samples$taxon))) %>% 
  rename(taxon = 1)

# load name repair file 

repair <- read_csv("./raw_data/epifauna/taxonomic_codes.csv")

# drop messy names, sum abundances of each taxon
all_epifauna <- sept20_samples %>% full_join(feb21_samples) %>% 
  left_join(repair) %>% select(-taxon) %>% 
  rename(taxon = taxon_repaired) %>% select(-notes) %>% 
  group_by(date, block, number, taxon) %>% 
  transmute(total_abund = sum(abund)) %>% ungroup()

block_design_useful <- block_design %>% select(new_no, new_block, treatment, tile_id)

# add in useful information about tile id and block for later analysis
all_epifauna_trts <- all_epifauna %>% rename(new_no = number, new_block = block) %>% left_join(block_design_useful) %>% 
  select(-new_no) %>% rename(block = new_block) %>% unique()

# save these data
write_csv(all_epifauna_trts, "./clean_data/SVSHS_epifauna_clean.csv")


##################
# Appendix plot Fig A3 of number of tiles in each treatment over time

# read in survey dates
survey_dates <- read_csv("./raw_data/design/SVSWS_survey_times.csv") %>% 
  # sometimes surveys occurred over multiple tides (if tide window was short)
  # here, we just take the date of survey completion as the survey date
  group_by(survey_no) %>% 
  summarize(date = max(date))

# read in the tile data
tile_numbers <- read_csv("./clean_data/SVSWS_survey_clean.csv") %>% 
  mutate(treatment = if_else(survey_no <= 10, trt_y1, treatment)) %>% 
  filter((survey_no > 10 & treatment %in% c("C","W")) == F) %>% 
  group_by(survey_no, treatment) %>% 
  # calculate the number of tiles present at each survey timepoint
  summarize(number_tiles = length(unique(tile_id))) %>% 
  left_join(survey_dates)

# for site visits where the number of tiles quickly changed (April 2020 - treatments changed for Year 2;
# September 2020 - half of tiles destructively sampled), add in these data manually
tiles_extrapts <- data.frame(treatment = c("C","W","CC","CW","WC","WW",
                                           "CC","CW","WC","WW"),
         number_tiles = c(41,46,22,19,20,25,12,10,11,16), 
         date = c(rep(ymd("2020-04-03"), times = 2, each =1),
                  rep(ymd("2020-04-03"), times = 4, each = 1),
                  rep(ymd("2020-09-14"), times = 4, each = 1)))

# join together the data
tile_numbers <- tile_numbers %>% full_join(tiles_extrapts) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

# create plot of changes in treatment sample sizes (tile numbers) over time
tile.time <- ggplot(tile_numbers, 
                    aes(x=date, y = number_tiles, col = treatment)) +
  geom_line(lwd = 0.8) +
  theme_classic() +
  scale_color_manual(values = pal.trt) +
  geom_vline(aes(xintercept = as.Date("2020-04-03")), col = "grey50", lwd = 0.8) +
  geom_vline(aes(xintercept = ymd("2019-08-27")), col = "grey60", lty = "dashed", lwd = 0.8) +
  geom_vline(aes(xintercept = ymd("2019-06-06")), col = "grey70", lty = "dotdash", lwd = 0.8) +
  labs(x = "Date", y = "Sample size", col = "Treatment", lty = "Treatment") +
  annotate("segment", x=as.Date("2020-09-14"), lwd = 0.8, xend = as.Date("2020-09-14"),
           y = 33, yend = 28, arrow = arrow(length = unit(0.08, "inches"))) +
  ylim(c(0,50))
tile.time

# save this figure
png("./figures/FigS3.png", res = 700, width = 6, height = 3.5, units = "in")
tile.time
dev.off()
