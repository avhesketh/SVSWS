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

## cleaned up! now need to get everything on the same page with tile numbers

# links overall unique tile ID to block and number within block
block_design <- read_csv("./raw_data/design/SVSHW_tilesetup.csv")

# filter out surveys after tiles were moved around
survey_postmove <- survey_data_qc2 %>% 
  filter(date > "2019-06-06") %>% 
  rename("new_no" = number,
         "new_block" = block) %>% 
  left_join(block_design) %>% 
  rename(angle = new_angle,
         shore_level = new_shore_level,
         block = new_block) %>% 
  select(-new_no, -survived, -original_no, -original_shore_level,
         -original_angle, -colour_y1, -colour_y2, -notes, -original_block)

survey_premove <- survey_data_qc2 %>% 
  filter(date <= "2019-06-06") %>% 
  rename("original_no" = number,
         "original_block" = block) %>% 
  left_join(block_design) %>% 
  rename(angle = original_angle,
         shore_level = original_shore_level,
         block = original_block) %>% 
  select(-original_no, -survived, -new_no, -new_shore_level,
         -new_angle, -colour_y1, -colour_y2, -notes, -new_block)

survey_all <- survey_postmove %>% full_join(survey_premove)

surveys_clean <- survey_all %>%
  separate(treatment, into = c("trt_y1", "trt_y2"), sep = 1, remove = FALSE)

# also need a vector associating sampling dates with a survey number

survey_list <- surveys_clean %>% select(date) %>% unique() %>% arrange(date)

survey_no <- c(1,1,2,2,3,4,4,5,5,6,6,7,8,8,9,10,11,12,13,14,15,16)

survey_join <- cbind(survey_list, survey_no)

surveys_clean2 <- surveys_clean %>% left_join(survey_join)

write_csv(surveys_clean2, "./clean_data/SVSHW_survey_clean.csv")

## clean infaunal data

sept20_samples <- read_csv("./raw_data/final_samples/SVSHW_20200914_infauna.csv") %>% mutate(date = "2020-09-14")
feb21_samples <- read_csv("./raw_data/final_samples/SVSHW_20210224_infauna.csv") %>% mutate(date = "2021-02-24")

unique(sept20_samples$taxon)

infauna_original <- as.data.frame(unique(c(sept20_samples$taxon, feb21_samples$taxon))) %>% 
  rename(taxon = 1)

infauna_original

infauna_repair <- c("Amphipoda", "Insecta","Littorina_scutulata","Littorina_sitkana",
                    "Lottia_paradigitalis","Polychaeta","Copepoda","Cirripedia_larva","Amphipoda",
                    "Mytilus_trossulus","Nemertean","Isopoda","Amphipoda_larva",
                    "Lottia_pelta","Anthopleura_elegantissima","Polychaeta","Isopoda","Lasaea_rubra",
                    "Lottia_digitalis","Onchidoris_bilamellata","Nemertean","Pagurus_hirsutiusculus",
                    "Lottia_spp_recruit","Polychaeta","Oedoparena_sp_larva","Insecta","Amphipoda","Lottia_spp_recruits",
                    "All","Cirripedia_larva","Polychaeta","Isopoda","Lottia_scutum","Neostylidium_eschrichtii",
                    "Platyhelminthes","Oedoparena_sp_larva","Oedoparena_sp_larva","Annelida","Insecta")

infauna_join <- cbind(infauna_original, infauna_repair)

all_infauna <- sept20_samples %>% full_join(feb21_samples) %>% left_join(infauna_join) %>% select(-taxon) %>% 
  rename(taxon = infauna_repair) %>% select(-notes) %>% group_by(date, block, number, taxon) %>% 
  transmute(total_abund = sum(abund)) %>% ungroup()

block_design_useful <- block_design %>% select(new_no, new_block, treatment, tile_id)

all_infauna_trts <- all_infauna %>% rename(new_no = number, new_block = block) %>% left_join(block_design_useful) %>% 
  select(-new_no) %>% rename(block = new_block) %>% unique()

write_csv(all_infauna_trts, "./clean_data/SVSHW_infauna_clean.csv")
