#abundance data stuff
pkgs <- c("tidyverse","glmmTMB")
lapply(pkgs, library, character.only = TRUE)

# read in dataframes

surveys <- read_csv("./raw_data/tile_surveys/SVSHW_surveys.csv") %>% 
  select(-X8) 

block_design <- read_csv("./raw_data/design/SVSHW_tilesetup.csv")

sum(pansy)

# need to first assign treatments to each tile, but need to take into account when they got moved

survey_postmove <- surveys %>% 
  filter(date > "2019-06-06") %>% 
  rename("new_no" = number,
         "new_block" = block) %>% 
  full_join(block_design) %>% 
  rename(number = new_no,
         block = new_block,
         angle = new_angle,
         shore_level = new_shore_level)

surveys_tiles <- surveys %>% 
  filter(date <= "2019-06-06") %>% 
  rename("original_no" = number,
         "original_block" = block) %>% 
  full_join(block_design) %>% 
  rename(number = original_no,
         block = original_block,
         angle = original_angle,
         shore_level = original_shore_level) %>% 
  full_join(survey_postmove)

surveys_clean <- surveys_tiles %>% 
  mutate(date = str_replace_all(date , c("2019-05-08" = "2019-05-09",
                                         "2019-06-04" = "2019-06-05",
                                         "2019-07-03" = "2019-07-04",
                                         "2019-07-17" = "2019-07-18", 
                                         "2019-07-30" = "2019-07-31",
                                         "2019-10-18" = "2019-10-20")),
          treatment = ifelse(colour_y1 == "white" & colour_y2 == "black", "AW",
                            ifelse(colour_y1 == "white" & colour_y2 == "white", "AA",
                                   ifelse(colour_y1 == "black" & colour_y2 == "white", "WA",
                                          ifelse(colour_y1 == "black" & colour_y2 == "black", "WW", NA))))) %>% 
  filter(treatment %in% c("WW","AW","WA","AA")) %>% 
  select(1:6,8,10,23)

#write_csv(surveys_clean, "./clean_data/SVSHW_survey_clean.csv")

# barnacle recruitment 2020

recruitment <- read_csv("./raw_data/tile_surveys/SVSHW_bncle_recruit.csv") %>% 
  mutate(treatment = ifelse(trt == "A" & consecutive == "N", "AW",
                            ifelse(trt == "A" & consecutive == "Y", "AA",
                                  ifelse(trt == "W" & consecutive == "N", "WA", "WW")))) 

# destructive sampling sept 2020

samples_2020 <- read_csv("./raw_data/final_samples/SVSHW_infaunal_sept2020.csv") %>% 
  select(-X5)

block_design_sep2020 <- read_csv("./raw_data/design/SVSHW_tilesetup.csv") %>% 
  mutate(treatment = ifelse(colour_y1 == "white" & colour_y2 == "black", "AW",
                            ifelse(colour_y1 == "white" & colour_y2 == "white", "AA",
                                    ifelse(colour_y1 == "black" & colour_y2 == "white", "WA", "WW")))) %>% 
  select(-original_no, -original_block, -original_shore_level, -original_angle,-original_herb_trt, -colour_y1, -colour_y2, -survived) %>% 
  rename(number = new_no, block = new_block, angle = new_angle, aspect = new_compass, shore_level = new_shore_level) %>% 
  na.omit()

samples_clean <- samples_2020 %>% 
  left_join(block_design_sep2020)

write_csv(samples_clean, "./clean_data/SVSHW_infauna_20200914.csv")

######## modeling stuff

balanus_y1 <- read_csv("./clean_data/SVSHW_survey_clean.csv") %>% 
  filter(species == "balanus" & date <= "2020-03-15") %>% 
  mutate(trty1 = substring(treatment, 1, 1),
         trty2 = substring(treatment, 2, 2),
         timesincestart = difftime(date, min(date), units = c("weeks")),
         remove = ifelse(count == "NA", TRUE, FALSE), 
         count = as.integer(count)) %>% 
  filter(remove == FALSE) %>% 
  select(-remove)

balanus_mar20 <- balanus_y1 %>% 
  filter(date == "2019-10-20")

M1_bal <- glmmTMB(count ~ trty1 + (1|block), 
                  data = balanus_mar20, 
                  dispformula = ~block,
                  family = nbinom1)

M1_bal_resdh <- simulateResiduals(M1_bal)
plot(M1_bal_resdh)

M2_bal <- glmmTMB(count ~ treatment + (1|block), 
                  data = balanus_mar20, 
                  family = nbinom1)

AIC(M1_bal, M2_bal) # M1 better

Anova(M1_bal)

## now for ulothrix

ulo_y1 <- read_csv("./clean_data/SVSHW_survey_clean.csv") %>% 
  filter(species == "ulothrix" & date <= "2020-03-15") %>% 
  mutate(trty1 = substring(treatment, 1, 1),
         trty2 = substring(treatment, 2, 2),
         timesincestart = difftime(date, min(date), units = c("weeks")),
         remove = ifelse(percent_cover == "NA", TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  select(-remove)

ulo_y1_oct <- ulo_y1 %>%
  filter(date == "2019-10-20") %>% 
  mutate(percent_cover = (percent_cover*(length(percent_cover) -1) + 0.5)/length(percent_cover)) %>% 
  mutate(percent_cover = as.numeric(percent_cover), trty1 = as.factor(trty1))

M1_ulo <- glmmTMB(percent_cover/100 ~ trty1 + (1|block),
                  data = ulo_y1_oct,
                  family = beta_family())

plot(simulateResiduals(M1_ulo))

Anova(M1_ulo)
