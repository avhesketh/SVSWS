## Community trajectories
## Amelia Hesketh: October 2021

# Step 1: generate matrices for ALL visual survey data(!)
# Will need to store these in a list where each entry is associated with a block and a survey date.

survey_data <- read_csv("./raw_data/tile_surveys/SVSHW_surveys.csv") %>%   
  mutate(species = str_replace_all(species, c("Lottia_sp" = "Lottia_sp_recruits", "scunge" = "Scunge")))

# First, generate a list of unique species (26 in this case)

sp_list <- survey_data %>% select(species) %>% unique()

# this species list will need to be joined to any matrices so that zeroes are 
# matrices are of the same dimensions in the final list.

# also need a vector associating sampling dates with a survey number

survey_list <- survey_data %>% select(date) %>% unique()

survey_no <- c(1,1,2,2,3,4,4,5,5,6,6,7,8,8,9,10,11,12,13,14,15,16)

survey_join <- cbind(survey_list, survey_no)

survey_data_numeric <- survey_data %>% full_join(survey_join)

# and one more df for max cover/abund for each species to allow standardizing

max_abund <- survey_data %>% group_by(species) %>% mutate(max_abund = if_else(is.na(count), max(percent_cover),
                                                                              max(count))) %>% 
  select(species, max_abund) %>% unique() %>% na.omit()
View(survey_data)
# Now, create a loop to generate this beautiful list.

for (survey in 1:16){
  
  # first, subset data by survey number
  survey = 1
  survey_subset <- survey_data_numeric %>% filter(survey_no == survey) %>% 
    select(block, species, count, percent_cover)
  
  # then, make a separate df for algae and inverts, and rename those columns to allow a join
  # might as well standardize the values at the same time to at least have on a similar scale...
  
  algae <- survey_subset %>% select(block, percent_cover, species) 
  
  
  %>% na.omit() %>%
    filter(percent_cover > 0) %>% 
    group_by(species) %>% mutate(max_cover = max(percent_cover)) %>% ungroup() %>% 
    mutate(st_abund = percent_cover/max_cover) %>% select(-percent_cover)
    
  inverts <- survey_subset %>% select(block, count, species) %>% na.omit() %>% 
    group_by(species) %>% mutate(max_count = max(count)) %>% ungroup() %>% 
    mutate(st_abund = count/max_count) %>% select(-count)
  
  # put these back together
  
  all_spp <- full_join(algae, inverts) %>% full_join(species_list)
  
}
?replace_na
        