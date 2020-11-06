library(vegan)
library(tidyverse)

# need to get community data in useful format
surveys <- read_csv("clean_survey_data.csv")

comm_october <- surveys %>%
  filter(date == "2019-10-20") %>%
  unite("abund", c(count, percent_cover), sep = "") %>% 
  select(-timediff, -date, -X1) %>% 
  unite("b_no", c(block, number), sep = "_")

comm_october$abund <- gsub("NA", "", comm_october$abund)
comm_october$abund <- as.numeric(comm_october$abund)

# get data in wide format with col = spp, row = tile
comm_spread <- spread(comm_october, key = species, value = abund)

#get vector of treatment temps
comm_trts <- comm_spread$colour

# change b_no to row names
comm_matrix <- comm_spread %>% 
  select(-colour) %>% 
  column_to_rownames("b_no")

# replace NA values with zero
comm_matrix[is.na(comm_matrix)] <- 0

# convert to matrix
comm_matrix <- as.matrix(comm_matrix)

# ready to analyze
comm_oct <- metaMDS(comm_matrix, k = 2, try = 100)

ordiplot(comm_oct, type = "p")
orditorp(comm_oct,display="species",col="red",air=0.01)
ordiellipse(comm_oct, groups = comm_trts, 
            col = c("grey", "black"), label = TRUE)


# surveys july - shannon diversity

survey_july <- read_csv("surveys.csv") %>% 
  filter(date > "20-07-01") %>% 
  select(-1, -6, -7,-8) %>% 
  mutate(count = as.numeric(count))

barn_jul_trt <- barn_jul %>% 
  select(-5,-6,-7) %>% 
  rename(number = "tile")

survey_july <- survey_july %>% 
  full_join(barn_jul_trt) %>% 
  unique()

spread_july <- spread(survey_july, key = species, value = count) %>% 
  select(-7, -9, -11) %>% 
  mutate(mytilus = replace_na(mytilus, 0)) %>% 
  na.omit()

spread_july_data <- spread_july[,5:8]

spread_july_div <- spread_july %>% 
  mutate(shannon = diversity(spread_july_data, index = "shannon"))

div_plot <- ggplot(aes(x = trt, fill = consecutive, y = shannon), data = spread_july_div) +
  geom_bar(stat = "identity", position = "dodge")
div_plot

hist(spread_july_div$shannon, breaks=100)

library(lme4)

div.glmm.1 <- lmer(shannon ~ trt*consecutive + (1|block),
                   data = spread_july_div)

plot(div.glmm.1)

summary(div.glmm.1)

Anova(div.glmm.1)
