library(dplyr)

entryS <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/Report 9-Volunteer Medical History - Report 9-Volunteer Medical History.csv")
mapping <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/Original Study Mapping - Sheet3.csv")

entryS <- entryS %>%
  left_join(mapping, by = c("Your.Biome.Health.App.ID" = "Biome.Health.App.ID")) %>%
  mutate(Your.Biome.Health.App.ID = STUDY.ID) %>%
  select(-STUDY.ID)
