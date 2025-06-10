################################################################################
## Escape Barrier Score Script
################################################################################

# if libraries are not installed, install them
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
if (!requireNamespace("pracma", quietly = TRUE)) {
  install.packages("pracma")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(readxl)
library(pracma)
library(tidyverse)

# load escape path data and filter out NA
escape_paths <- read_excel("data/Path Coordinates.xlsx") %>% filter(IC50 > 0)

# load escape path frequency data and filter out NA / non filled data
esacpe_freq <- read_excel("data/Path Frequencies.xlsx") %>% filter(Frequency > 0)

# calculate the area under the curve for each pathway
# return as a dataframe with pathway and AUC
escape_paths <- escape_paths %>%
  group_by(PathID) %>%
  summarise(AUC = trapz(IC50, GrowthRate))

# join to escape_freq data
escape_paths <- left_join(escape_paths, esacpe_freq, by = "PathID")

# compute the path specific AUC
escape_paths$NormAUC <- escape_paths$AUC * escape_paths$Frequency

# update the column order
escape_paths <- escape_paths %>%
  select(ID,PathID, AUC, Frequency, NormAUC, Virus, Antibody)

# sort by ID 
escape_paths <- escape_paths %>%
  arrange(ID)

write.csv(escape_paths, "data/PathEscapeBarrierScore.csv", row.names = FALSE)

# group by virus antibody combinations, 
# return only unique antibody virus combinations and the escapability score
escapebarrierscore <- escape_paths %>%
  mutate(NormAUC = AUC *Frequency) %>% 
  group_by(Virus, Antibody) %>%
  mutate(Escape_Barrier_Score = sum(NormAUC)) %>% 
  dplyr::select(ID,Virus, Antibody, Escape_Barrier_Score) %>%
  unique() %>%
  drop_na()

# sort by ID
escapebarrierscore <- escapebarrierscore %>%
  arrange(ID)

write.csv(escapebarrierscore, "data/TotalEscapeBarrierScore.csv", row.names = FALSE)


