## QC the raw data

survey_data <- read_csv("./raw_data/tile_surveys/SVSHW_surveys.csv",
                        col_types = c("D","f","f","c","n","n","c","?")) %>% 
  mutate(species = str_replace_all(species, c("Lottia_sp" = "Lottia_sp_recruits", "scunge" = "Scunge"))) %>% 
  select(-8)

library(assertr)

id <- c("i","i","a","i", "i", "i", "i", "i", "i", "i", "a", 
        "a","i","i","i","i","i","i", "a","a","i","a","a","a","a", "a")

classify_sp <- cbind(sp_list, id)

invert_list <- classify_sp %>% filter(id == "i")
alga_list <- classify_sp %>% filter(id == "a")

str(survey_data)

survey_data %>% filter(species %in% invert_list$species) %>% verify(is.na(percent_cover)) %>% verify(count >= 0)

survey_data_qc <- survey_data %>% mutate(remove = if_else(species %in% invert_list$species & is.na(count),
                                                         T, F)) %>% filter(remove == F) %>% select(-remove)

survey_data_qc %>% filter(species %in% invert_list$species) %>% verify(count >= 0)

survey_data_alga <- survey_data_qc %>% filter(species %in% alga_list$species) 

survey_data_alga %>% verify(percent_cover >= 0)

survey_data_qc2 <- survey_data_qc %>% mutate(percent_cover = if_else(species %in% alga_list$species & (is.na(count) == F), 
                                                                    count, percent_cover)) %>% 
  mutate(remove = if_else(species %in% alga_list$species & is.na(percent_cover), T, F)) %>% filter(remove == F) %>% 
           select(-remove)

survey_data_alga <- survey_data_qc2 %>% filter(species %in% alga_list$species) 

survey_data_alga %>% verify(percent_cover >= 0)

## change to final block...

