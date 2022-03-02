library(tidyverse)

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

ulo_plot <- read_csv("./plotting_df/ulo_time.csv") %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

treatment_colours <- c("cornflowerblue","tomato3","blue","purple","orange","darkred")

ulo_time <- ggplot(aes(x = date, y = mean_abund, color = treatment), 
                       data = ulo_plot) +
  geom_point() +
  geom_line()+
  theme_classic() +
  scale_color_manual(values = treatment_colours) +
  geom_errorbar(aes(ymax = mean_abund + se_abund, ymin = mean_abund - se_abund),
                width = 7) +
  labs(x = "Date", y = expression("Mean cover of"~ italic(Ulothrix ~sp.)~ "(%)"),
       col = "Treatment") +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) 
ulo_time


# all algae -- y1 driven almost entirely by Ulothrix ephemeral green.
algae_plot <- read_csv("./plotting_df/algal_cover.csv") %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

cover <- ggplot(aes(x = date, y = mean_cover, col = treatment), data = algae_plot) +
  theme_classic() +
  geom_point()+
  geom_line() +
  geom_errorbar(aes(ymax = mean_cover + se_cover, ymin = mean_cover - se_cover)) +
  labs(x = "Date", y ="Algal cover (%)", col = "Treatment") +
  scale_color_manual(values = treatment_colours) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 16))

  
# limpets


limpet_plot <- limpets_y2 %>% group_by(date, treatment) %>% 
  summarize(mean_limpets = mean(total_limpets),
            se_limpets = std.error(total_limpets))

ggplot(limpet_plot, aes(x=date, y=mean_limpets, col = treatment)) +
  geom_point()+
  geom_line() +
  theme_classic() +
  labs(x = "Date", y = "Number of limpets", col = "Treatment") +
  scale_color_manual(values = c("blue","purple","orange","darkred")) +
  geom_errorbar(aes(ymax = mean_limpets + se_limpets, 
                    ymin = mean_limpets - se_limpets), width = 5) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  

litt_plot <- litts_y2 %>% group_by(date, treatment) %>% 
  summarize(mean_litts = mean(total_litts),
            se_litts = std.error(total_litts))

ggplot(litt_plot, aes(x=date, y=mean_litts, col = treatment)) +
  geom_point()+
  geom_line() +
  theme_classic() +
  labs(x = "Date", y = "Number of littorines", col = "Treatment") +
  scale_color_manual(values = c("blue","purple","orange","darkred")) +
  geom_errorbar(aes(ymax = mean_litts + se_litts, 
                    ymin = mean_litts - se_litts), width = 5) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 16))

# diversity

diversity.plot <- read_csv("./plotting_df/diversity.csv") %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))
View(diversity.plot)
richness <- ggplot(diversity.plot, aes(y = mean.rich, x =date, col = treatment)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = treatment_colours) +
  geom_errorbar(aes(ymax = mean.rich + se.rich, ymin = mean.rich - se.rich)) +
  labs(y = "Species richness", x = "Date", col = "Treatment") +
  theme_classic()+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 16))

shannon.i <- ggplot(diversity.plot, aes(y = mean.shann.i, x =date, col = treatment)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = treatment_colours) +
  geom_errorbar(aes(ymax = mean.shann.i + se.shann.i, ymin = mean.shann.i - se.shann.i)) +
  labs(y = "Invertebrate Shannon diversity", x = "Date", col = "Treatment") +
  theme_classic()+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 16))
shannon.i

shannon.a <- ggplot(diversity.plot, aes(y = mean.shann.a, x =date, col = treatment)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = treatment_colours) +
  geom_errorbar(aes(ymax = mean.shann.a + se.shann.a, ymin = mean.shann.a - se.shann.a)) +
  labs(y = "Algal Shannon diversity", x = "Date", col = "Treatment") +
  theme_classic()+
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 16))
shannon.a

