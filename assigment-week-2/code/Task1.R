#############################
# Task 1: Summarizing the data frame and simple ggplot with chatGPT (Fig. 1b) (Chung)
# Author: María José Pino
# Date: 2024-10-03
#############################

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(ggpubr)

# Load the data
data <- read_delim("assigment-week-2/data/Task1/00_data/aCRE.IPSC.info.tsv", delim = "\t")


##### Data formatting #####
summary_data <- data %>%
  group_by(transcription, proximity, genic) %>%
  summarise(count = n(), .groups = "drop") %>%  # count occurrences
  arrange(desc(count))

data_all <- data %>%
  group_by(proximity, genic) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(pct = count / sum(count) * 100,
         scope = "All",
         region_type = ifelse(proximity == "proximal", "proximal",
                  ifelse(proximity == "distal" & genic == "intergenic", "distal_intergenic",
                  ifelse(proximity == "distal" & genic != "intergenic", "distal_genic", NA))))


data_trascription_distal_intergenic <- data %>%
    filter(proximity == "distal" & genic == "intergenic") %>%
    group_by(transcription) %>%
    mutate(transcription_ann = case_when(transcription == "none-Trn" ~ "untranscribed",
                                    transcription == "weak-Trn" ~ "ambiguous",
                                    transcription == "firm-Trn" ~ "transcribed", 
                                    TRUE ~ NA_character_
                                    )) %>%
    group_by(transcription_ann) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(pct = count / sum(count) * 100,
         scope = "Distal Intergenic",
         region_type = "distal_intergenic")


##### Plots ######
## All cREs ##
# create a stack plot with ggplot2
ggplot(data_all, aes(x = scope, y = count, fill = region_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("proximal" = "blue", "distal_genic" = "green", "distal_intergenic" = "orange")) +
  labs(title = "Distribution of aCREs by Region Type",
       x = "Scope",
       y = "Percentage of aCREs",
       fill = "Region Type") +
  theme_minimal() +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
                position = position_fill(vjust = 0.5),
                color = "white",
                size = 5)




plot_all <- ggplot(data_all, aes(x = scope, y = pct, fill = factor(region_type, levels = c("distal_intergenic", "distal_genic", "proximal")))) +  
  geom_bar(stat = "identity") +  
  #scale_y_continuous(labels = scales::percent_format()) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5), size = 3, fontface = "bold", color = "white") +  
  labs(x = NULL, y = "Percentage (%)", fill = "Region Type", title = "Locations") +  
  theme_minimal() +  
  scale_fill_manual(values = c("#fdbf6f", "#ff7f00", "#984ea3")) +  
  theme(  
    plot.title = element_text(hjust = 0.5)  
  )  

ggsave(plot_all, filename = "assigment-week-2/results/Task1/plot_all.pdf", width = 5, height = 5, dpi = 300)


## Distal Intergenic cREs ##
# create a stack plot with ggplot2
ggplot(data_trascription_distal_intergenic, aes(x = scope, y = count, fill = transcription_ann)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("untranscribed" = "red", "ambiguous" = "purple", "transcribed" = "yellow")) +
  labs(title = "Distribution of Distal Intergenic aCREs by Transcription Type",
       x = "Scope",
       y = "Percentage of aCREs",
       fill = "Transcription Type") +
  theme_minimal() +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
                position = position_fill(vjust = 0.5),
                color = "white",
                size = 5)

plot_intergenic <- ggplot(data_trascription_distal_intergenic, aes(x = scope, y = pct, fill = factor(transcription_ann, levels = c( "untranscribed", "ambiguous", "transcribed")))) +  
  geom_bar(stat = "identity") +  
  #scale_y_continuous(labels = scales::percent_format()) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5), size = 3, fontface = "bold", color = "white") +  
  labs(x = NULL, y = "Percentage (%)", fill = "Transcription Type", title = "Transcription") +  
  theme_minimal() +  
  scale_fill_manual(values = c("#377eb8", "#4daf4a", "#e41a1c")) +  
  theme(  
    plot.title = element_text(hjust = 0.5)  
  )

ggsave(plot_intergenic, filename = "assigment-week-2/results/Task1/plot_intergenic.pdf", width = 5, height = 5, dpi = 300)


combined_plot <- ggarrange(plot_all, plot_intergenic, ncol = 2)  
ggsave("assigment-week-2/results/Task1/task1.stack_bar.pdf", plot = combined_plot, width = 8, height = 4)  
