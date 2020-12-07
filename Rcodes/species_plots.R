## plotting patterns in counts

balanus_plot <- read_csv("./clean_data/SVSHW_survey_clean.csv") %>% 
  filter(date == "2019-10-20" & species == "balanus") %>% 
  select(-date) %>% 
  mutate(trty1 = substring(treatment, 1, 1)) %>% 
  group_by(trty1) %>% 
  summarize(av_count = mean(count), se_count = sd(count)/sqrt(length(count))) %>% 
  rename(Treatment = trty1)

balanus_time <- ggplot(aes(x = Treatment, y = av_count, fill = Treatment), 
                       data = balanus_plot) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("blue", "red")) +
  geom_errorbar(aes(ymax = av_count + se_count, ymin = av_count - se_count),
                width = 0.3, position = "dodge") +
  xlab("Treatment") +
  ylab(expression("Mean abundance of"~ italic(B. ~glandula))) +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16))
balanus_time

# need to change dates to be equal (surveys on same date)

ulo_plot <- read_csv("./clean_data/SVSHW_survey_clean.csv") %>% 
  filter(date == "2019-10-20" & species == "ulothrix") %>% 
  mutate(trty1 = substring(treatment, 1, 1)) %>% 
  group_by(trty1) %>% 
  summarize(av_cover = mean(percent_cover, na.rm = TRUE), 
            se_cover = sd(percent_cover)/sqrt(length(percent_cover))) %>% 
  rename(Treatment = trty1)

ulo_time <- ggplot(aes(x = Treatment, y = av_cover, fill = Treatment), 
                       data = ulo_plot) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("blue","red")) +
  geom_errorbar(aes(ymax = av_cover + se_cover, ymin = av_cover - se_cover),
                width = 0.7) +
  xlab("Treatment") +
  ylab(expression("Mean cover of"~ italic(Ulothrix ~sp.)~ "(%)")) +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) 
ulo_time


