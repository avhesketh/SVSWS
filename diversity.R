## Diversity

surveys <- read_csv("./clean_data/SVSHW_survey_clean.csv")

richness <- surveys %>% left_join(survey_join) %>% 
  group_by(block, tile_id, survey_no) %>% 
  summarize(richness = length(species))

invert_diversity <- surveys %>% left_join(survey_join) %>% 
  filter(is.na(count) == F) %>% slice(-4559) %>% select(tile_id, survey_no, block, species, count, treatment) %>%
  pivot_wider(names_from = species, values_from = count) %>% 
  mutate_all(~replace_na(.,0))

factors_diversity <- invert_diversity %>% select(tile_id, survey_no, block, treatment)

shannon <- diversity(invert_diversity)

diversity_df <- as.data.frame(matrix(ncol = 3))
colnames(diversity_df) <- c("richness", "shannon", "evenness")

for (tile in 1:100){
  for (survey in 1:16){
    
    tile_matrix <- invert_diversity %>% filter(tile_id == tile, survey_no == survey) %>% 
      select(-tile_id, -survey_no, -block, -treatment) %>% as.matrix() 
    
    if(is.numeric(tile_matrix)==T){
      
      richness <- specnumber(tile_matrix)
      shannon <- diversity(tile_matrix)
      evenness <- shannon/log(richness)
      
      new_row <- c(richness, shannon, evenness)
      diversity_df <- rbind(diversity_df, new_row)
    }
  }
}

View(survey_join)
  
diversity_omit <- diversity_df %>% slice(-which(is.na(richness)))

diversity_join <- cbind(factors_diversity, diversity_omit)

ggplot(aes(x = treatment, y = richness, col = treatment), data = diversity_join %>% filter(survey_no == 14)) +
  geom_boxplot()

library(lme4)

shannon.y2.summ <- lm(evenness ~ treatment, data = diversity_join %>% filter(survey_no == 14))
summary(shannon.y2.summ)
plot(shannon.y2.summ)
Anova(shannon.y2.summ)

shannon.y2.wint <- lm(evenness ~ treatment, data = diversity_join %>% filter(survey_no == 16))
summary(shannon.y2.wint)
Anova(shannon.y2.wint)
plot(shannon.y2.wint)



