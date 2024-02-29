## Community trajectories & analysis
## Amelia Hesketh: October 2021

pkgs <- c("tidyverse","lubridate","vegan", "ecotraj")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

survey_data <- read_csv("./clean_data/SVSWS_survey_clean.csv")


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

survey_data_spread2 <- survey_data %>% 
  filter(survey_no %in% 1:10) %>% mutate(treatment = if_else(treatment %in% c("WW","WC"), "W",
                                                             if_else(treatment %in% c("CC","CW"), "C", treatment))) %>% 
  unique() %>% select(-block) %>% 
  group_by(treatment, survey_no, species) %>% 
  reframe(mean_abund = if_else(is.na(count), mean(percent_cover),mean(count))) %>%
  unique() %>%  
  pivot_wider(names_from = species, values_from = mean_abund, values_fill = 0,
              id_cols = c("treatment","survey_no"))

# Part 2: isolate year 2 & average for each treatment
survey_data_spread3 <- survey_data %>%
  filter(survey_no %in% 10:16 & treatment %in% c("CC","WC","CW","WW")) %>% 
  unique() %>% select(-block) %>% 
  group_by(treatment, survey_no, species) %>% 
  reframe(mean_abund = if_else(is.na(count), mean(percent_cover),mean(count))) %>% unique() %>%  
  pivot_wider(names_from = species, values_from = mean_abund, values_fill = 0,
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

#sites <- sites %>% as.data.frame(sites) %>% rename(site_codes = 1)
#write_csv(sites, "./plotting_df/cta_sites.csv")

# Create survey vector
surveys <- survey_data_spread_smush$survey_no
treatment <- survey_data_spread_smush$treatment

#surveys <- surveys %>% as.data.frame %>% rename(survey_codes = 1)
#write_csv(surveys, "./plotting_df/cta_surveys.csv")

# Create matrix for metaMDS
cta_matrix <- survey_data_spread_smush %>% ungroup() %>% select( -treatment, -survey_no) %>% as.matrix()

#write_csv(cta_matrix, "./plotting_df/cta_matrix.csv")

# Get distances in nMDS
cta_bray <- dbrda(cta_matrix ~ treatment, distance = "bray", k = 2, try = 999, autotransform = T)

pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B", "grey50")
lty.trt <- c("twodash","twodash","solid","solid","solid","solid")

scores(cta_bray)$sites

bray.plot <- scores(cta_bray)$sites %>% cbind(survey_data_spread_smush$treatment) %>% 
  cbind(survey_data_spread_smush$survey_no) %>% as_tibble(.) %>% rename(treatment = 3, survey_no = 4) %>% 
  mutate(dbRDA1 = as.numeric(dbRDA1),dbRDA2 = as.numeric(dbRDA2)) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

for (row in 1:nrow(bray.plot)){
  if (row < nrow(bray.plot)){
    bray.plot$x.end[row] <- bray.plot$dbRDA1[row+1]
    bray.plot$y.end[row] <- bray.plot$dbRDA2[row+1]
  }
}

bray.plot.2 <- bray.plot %>% slice(c(-10,-20,-27,-34,-41,-48))

cta <- ggplot(bray.plot.2, aes(x = dbRDA1, y =dbRDA2, xend = x.end, yend = y.end, col = treatment, lty = treatment)) +
  geom_segment(arrow = arrow(type = "closed", ends = "last", length = unit(0.1,"in"))) +
  scale_linetype_manual(values = lty.trt) +
  scale_color_manual(values = pal.trt) +
  labs(lty = "Treatment", col = "Treatment") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title= element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  geom_point(data=bray.plot.2 %>% filter(survey_no ==1), col = "black", pch = 1, size = 3,
             stroke = 1) +
  annotate("point", y = -0.9, x = -1.9, pch = 1, size = 3) +
  annotate("text", y = -0.9, x = -1.4, label = "April 2019 (start y1)", size = 4) +
  geom_point(data=bray.plot.2 %>% filter(survey_no ==7), col = "black", pch = 0, size = 3,
             stroke = 1) +
  annotate("point", y = -1.05, x = -1.9, pch = 0, size = 3) +
  annotate("text", y = -1.05, x = -1.55, label = "August 2019", size = 4) +
  geom_point(data=bray.plot.2 %>% filter(survey_no ==10), col = "black", pch = 2, size = 3,
             stroke = 1) +
  annotate("point", y = -1.2, x = -1.9, pch = 2, size = 3) +
  annotate("text", y = -1.2, x = -1.4, label = "April 2020 (start y2)", size = 4) +
  geom_point(data=bray.plot.2 %>% filter(survey_no == 14), col = "black", pch = 5, size = 3,
             stroke = 1) +
  annotate("point", y = -1.35, x = -1.9, pch = 5, size = 3) +
  annotate("text", y = -1.35, x = -1.47, label = "September 2020", size = 4) +
  geom_point(data=bray.plot.2 %>% filter(survey_no == 15), aes(x = x.end, y = y.end),
             col = "black", pch = 13, size = 3,
             stroke = 1) +
  annotate("point", y = -1.5, x = -1.9, pch = 13, size = 3) +
  annotate("text", y = -1.5, x = -1.39, label = "February 2021 (end)", size = 4)
cta

ggsave(cta, filename = "./figures/FigA11.png",dpi = 1000, width = 3.5, height = 2.5, units = "in",scale = 2)

