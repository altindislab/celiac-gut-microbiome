##For 2 Years

a <- ultimate_data

postsort <- filter(a, a$Sorted== "Postsort")

acols <- select(postsort, ASV, Sample , Abundance ,Age, Disease.status, Sort, Group, Taxa)

gsub_deneme <- acols

#grpdeneme <- gsub("\\.[1-2]$", "", gsub("\\.[1-2]\\.", ".", acols$Sample))
#grpdeneme <- gsub("\\.[1-2]$", "333.", acols$Sample)

library(magrittr)

x <- "a'b c"

gsub_deneme$Sample %<>%
  gsub(".pos", "", .) %>%
  gsub(".neg", "", .) 

heat_data <- gsub_deneme

#heat_data_igpos <- filter(heat_data, heat_data$Sort== "IGpos")
#heat_data_igneg <- filter(heat_data, heat_data$Sort== "IGneg")

hip <- filter(heat_data, heat_data$Sort== "IGpos")
hin <- filter(heat_data, heat_data$Sort== "IGneg")

deneme <-merge(hip, hin, by = c("ASV", "Sample", "Disease.status", "Taxa"), all = TRUE)

deneme23 <- deneme %<>% mutate(ICI = Abundance.x/Abundance.y)

all_asv_ICI_values <- deneme23

hd <- all_asv_ICI_values

write.csv(hd, "ALL_ICI.csv", na = "")


ici_data <-  mutate(ICI = hip$Abundance/hin$Abundance)
ICI = hip$Abundance/hin$Abundance




#hip_celiac <- filter(hip, hip$Disease.status== "Celiac")
#hin_celiac <- filter(hin, hin$Disease.status== "Celiac")

#hip_control <- filter(hip, hip$Disease.status== "Celiac")
#hin_celiac <- filter(hin, hin$Disease.status== "Celiac")


