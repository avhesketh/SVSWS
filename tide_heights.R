## June 2019 It's time to think about tide heights

library(tidyverse)

## read csv file

heights <- read_csv("tide_heights.csv")

##

## now plot it

height_plot <- ggplot(aes(x=number, y = tile_cd, colour = block), data = heights) +
  geom_point() +
  geom_smooth(se = F) +
  ylab("Shore height (m)") +
  xlab("Tile number")
height_plot
