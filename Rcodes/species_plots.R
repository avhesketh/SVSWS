## plotting patterns in counts

balanus <- read_csv("balanus.csv")

balanus_plot <- balanus %>% 
  select(-date) %>% 
  group_by(colour, timediff) %>% 
  summarize(av_count = mean(count), se_count = sd(count)/sqrt(length(count))) %>% 
  mutate(exp_time = exp(-1/timediff))

balanus_time <- ggplot(aes(x = exp_time, y = av_count, colour = colour), 
                       data = balanus_plot) +
  geom_point(size = 2.7) +
  theme_classic() +
  scale_colour_manual(values = c("grey80", "grey20")) +
  geom_errorbar(aes(ymax = av_count + se_count, ymin = av_count - se_count),
                width = 0.01) +
  xlab("exp(-1/Time since experiment start (weeks))") +
  ylab(expression("Mean number of"~ italic(B. ~glandula))) +
  labs(colour = "Temperature") +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16))
balanus_time
?scale_colour_manual
# need to change dates to be equal (surveys on same date)

ulothrix <- read_csv("ulothrix.csv")

ulo_plot <- ulothrix %>% 
  dplyr::select(-X1, -X1_1, -date) %>% 
  group_by(colour, timediff) %>% 
  summarize(av_cover = mean(percent_cover), 
            se_cover = sd(percent_cover)/sqrt(length(percent_cover))) %>% 
  mutate(exp_time = exp(-1/timediff))

ulo_time <- ggplot(aes(x = timediff, y = av_cover, colour = colour), 
                       data = ulo_plot) +
  geom_point(size = 2.7) +
  theme_classic() +
  scale_colour_manual(values = c("grey80", "grey20")) +
  geom_errorbar(aes(ymax = av_cover + se_cover, ymin = av_cover - se_cover),
                width = 0.7) +
  xlab("Time since experiment start (weeks)") +
  ylab(expression("Mean cover of"~ italic(Ulothrix ~sp.)~ "(%)")) +
  labs(colour = "Temperature") +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  geom_smooth(aes(y = pred, x = timediff), data = ulothrix)
ulo_time


