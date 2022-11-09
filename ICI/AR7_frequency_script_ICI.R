###hd is heatmap all data
#lets filter ICI>=2
library(plyr)
hd_filtered <- filter(hd, hd$ICI>=0)

##hd2 is age 2.5

hd2 <- filter(hd_filtered, hd_filtered$Age.x=="2.5year")

count_age2 <- count(hd2, c("ASV","Disease.status"))
count_age2 <- filter(count_age2, count_age2$freq>=2)

hd2_with_freq <- semi_join(hd2, count_age2 , by = "ASV")#by = c("ASV","Disease.status"))


###Take those ASV from initial dataframe
selected_ASV_from_initial_data <- semi_join(hd, hd2_with_freq, by ="ASV")

###Select AGE 2.5
selected_ASV_from_initial_data_age_2_5 <- filter(selected_ASV_from_initial_data, selected_ASV_from_initial_data$Age.x=="2.5year")


#######AGE 5
##hd5 is age 5

hd5 <- filter(hd_filtered, hd_filtered$Age.x=="5year")

count_age5 <- count(hd5, c("ASV","Disease.status"))
#count_age5 <- filter(count_age5, count_age5$freq>=1)

hd5_with_freq <- semi_join(hd5, count_age5 , by = "ASV")#by = c("ASV","Disease.status"))


###Take those ASV from initial dataframe
selected_ASV_from_initial_data <- semi_join(hd, hd5_with_freq, by ="ASV")

###Select AGE 2.5
selected_ASV_from_initial_data_age_5 <- filter(selected_ASV_from_initial_data, selected_ASV_from_initial_data$Age.x=="5year")

