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


