## Community trajectories & analysis
## Amelia Hesketh: October 2021

pkgs <- c("tidyverse","lubridate","vegan", "ecotraj")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

survey_data <- read_csv("./clean_data/SVSHW_survey_clean.csv")

# First, generate a vector for treatments...

trt_tile <- survey_data %>% select(tile_id, treatment) %>% unique()

# and make one more df for max cover/abund for each species to allow standardizing

max_abund <- survey_data %>% group_by(species) %>% mutate(max_abund = if_else(is.na(count), max(percent_cover),
                                                                              max(count))) %>% 
  select(species, max_abund) %>% unique() %>% na.omit()

# Now, create a loop to generate this beautiful list.

survey_list <- list()

survey_data_st <- data.frame()

for (survey in 1:16){
  
  # first, subset data by survey number
  survey_subset <- survey_data %>% filter(survey_no == survey) %>% 
    select(block, species, count, percent_cover, tile_id, survey_no)
  
  # then, make a separate df for algae and inverts, and rename those columns to allow a join
  # might as well standardize the values at the same time to at least have on a similar scale...
  
  algae <- survey_subset %>% select(block, tile_id, percent_cover, species) %>%
    left_join(max_abund) %>% 
    group_by(species) %>%
    mutate(st_abund = percent_cover/max_abund) %>% select(-percent_cover, -max_abund) %>% na.omit
    
  inverts <- survey_subset %>% select(block, tile_id, count, species) %>%
    left_join(max_abund) %>% 
    group_by(species) %>% 
    mutate(st_abund = count/max_abund) %>% select(-count, -max_abund) %>% na.omit
  
  # put these back together
  
  all_spp <- full_join(algae, inverts) %>% mutate(survey_no = survey) %>% left_join(trt_tile)
  
  survey_data_st <- rbind(survey_data_st, all_spp)
  
}


#want to try a plot that shows equivalence of CC/CW and WC/WW for the first year.

# Part 1: isolate year 1 & assign to the two treatments before averaging community composition at each timepoint

survey_data_spread2 <- survey_data_st %>% relocate(treatment, survey_no, .before = st_abund) %>%
  filter(survey_no %in% 1:10) %>% mutate(treatment = if_else(treatment %in% c("WW","WC"), "W",
                                                             if_else(treatment %in% c("CC","CW"), "C", treatment))) %>% 
  unique() %>% select(-block) %>% 
  group_by(treatment, survey_no, species) %>% summarize(st_abund = mean(st_abund)) %>% ungroup() %>% 
  pivot_wider(names_from = species, values_from = st_abund, values_fill = 0,
              id_cols = c("treatment","survey_no"))

# Part 2: isolate year 2 & average for each treatment
survey_data_spread3 <- survey_data_st %>% relocate(treatment, survey_no, .before = st_abund) %>%
  filter(survey_no %in% 10:16 & treatment %in% c("CC","WC","CW","WW")) %>% 
  unique() %>% select(-block) %>% 
  group_by(treatment, survey_no, species) %>% summarize(st_abund = mean(st_abund)) %>% ungroup() %>% 
  pivot_wider(names_from = species, values_from = st_abund, values_fill = 0,
              id_cols = c("treatment","survey_no"))

# Part 3: join together, replace NA with zeroes
survey_data_spread_smush <- full_join(survey_data_spread2, survey_data_spread3) %>% mutate_all(~replace_na(.,0))

# Create site vector
sites <- survey_data_spread_smush$treatment

sites <- gsub("CC", 3, sites, fixed = T)
sites <- gsub("CW", 4, sites, fixed = T)
sites <- gsub("WC", 5, sites, fixed = T)
sites <- gsub("WW", 6, sites, fixed = T)
sites <- gsub("C", 1, sites, fixed = T)
sites <- gsub("W", 2, sites, fixed = T)

# Create survey vector
surveys <- survey_data_spread_smush$survey_no

# Create matrix for metaMDS
cta_matrix <- survey_data_spread_smush %>% select( -treatment, -survey_no) %>% as.matrix()

# Get distances in nMDS
test <- metaMDS(cta_matrix, distance = "bray", k = 2, try = 9999)

# Generate plot!
oct_ord <- ordiplot(test, type = "none")
trajectoryPlot(test$points, sites, surveys, traj.colors = c("cornflowerblue","tomato3","blue","purple","orange","firebrick4"), lwd = 2)
legend("right", col = c("cornflowerblue","tomato3","blue","purple","orange","darkred"), legend = c("C","W","CC", "CW", "WC", "WW"), 
       bty = "n", lty = 1, lwd = 2)
points(x = c(0.923626972905131, 1.73887140611446), y = c(0.820536487858217, 0.911956147070956), pch = 4, lwd = 2)
points(x = c(-1.15486748997812, -0.576866176698416, -0.441740711707467, -0.0165466292382169), 
       y = c(0.320257210181677, 0.213559718152996, 0.284303427786167, -0.0123149647880705), pch = 8, lwd = 2)
points(x = c(-0.532798093768667, -0.877810259215775, -0.42835246500082848, -0.258833574482709),
       y = c(0.683368742327626, 0.262911229829701, 0.502851592978135, -1.32146970349316), pch = 0, lwd = 2)
orditorp(oct_ord,display="species",col="grey30", cex = 0.8)


### Delving into infaunal data

infauna <- read_csv("./clean_data/SVSHW_infauna_clean.csv") %>% filter(taxon != "All")

infauna_spread <- infauna %>% pivot_wider(id_cols = c(block, tile_id, date, treatment), 
                                          names_from = taxon, values_from = total_abund, values_fill = 0)

infauna_fall <- infauna_spread %>% filter(date == "2020-09-14") %>% select(-date, -tile_id)
infauna_spring <- infauna_spread %>% filter(date == "2021-02-24")  %>% select(-date, -tile_id)

blocks_fall <- infauna_fall$block
blocks_spring <- infauna_spring$block

treatments_fall <- as.data.frame(infauna_fall$treatment) %>% rename(treatment = 1) %>% 
  mutate(trt_y1 = substr(treatment, 1,1), trt_y2 = substr(treatment,2,2))

treatments_spring <- as.data.frame(infauna_spring$treatment) %>% rename(treatment = 1) %>% 
  mutate(trt_y1 = substr(treatment, 1,1), trt_y2 = substr(treatment,2,2))

infauna_fall_matrix <- infauna_fall %>% select(-block, -treatment) %>% as.matrix
infauna_spring_matrix <- infauna_spring %>% select(-block, -treatment) %>% as.matrix

#1: fall 2020 nMDS

fall2020 <- metaMDS(infauna_fall_matrix, distance = "bray", k = 2, try = 9999)

sept_ord <- ordiplot(fall2020, type = "none", xlim=c(-1,1),ylim = c(-1.5,1.5))
points(sept_ord, what = "sites", col = "black")
orditorp(sept_ord, display="species",col="grey30", cex = 0.8, air=1)
ordiellipse(sept_ord, groups = treatments_fall$treatment, 
            col = c("blue", "purple","orange","darkred"))
legend("right", col = c("blue","purple","orange","darkred"), legend = c("CC", "CW", "WC", "WW"), 
       bty = "n", lty = 1, lwd = 2)

permanova.fall2020 <- adonis(infauna_fall_matrix ~ trt_y1*trt_y2, data = treatments_fall, perm = 9999)
permanova.fall2020


disp.fall2020 <- betadisper(vegdist(infauna_fall_matrix, method = "bray"), group = treatments_fall$treatment)
disp.fall2020
anova(disp.fall2020)

#2: spring 2021 nMDS

spring2021 <- metaMDS(infauna_spring_matrix, distance = "bray", k = 2, try = 9999)

feb_ord <- ordiplot(spring2021, type = "none", xlim=c(-1,1),ylim = c(-1.5,1.5))
points(feb_ord, what = "sites", col = "black")
orditorp(feb_ord, display="species",col="grey30", cex = 0.8, air=1)
ordiellipse(feb_ord, groups = treatments_spring$treatment, 
            col = c("blue", "purple","orange","darkred"))
legend("right", col = c("blue","purple","orange","darkred"), legend = c("CC", "CW", "WC", "WW"), 
       bty = "n", lty = 1, lwd = 2)


permanova.spring2021 <- adonis(infauna_spring_matrix ~ trt_y1*trt_y2, data = treatments_spring, perm = 9999)
permanova.spring2021

disp.spring2021 <- betadisper(vegdist(infauna_spring_matrix, method = "bray"), group = treatments_spring$treatment)
disp.spring2021
anova(disp.spring2021)
