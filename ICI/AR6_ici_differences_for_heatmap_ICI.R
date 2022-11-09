

age2_data <- inner_join(age_2_celiac_ICI, age_2_control_ICI, by = c('ASV','asv_label','Taxa') , copy =T)
age2_data = age2_data  %>% mutate(ICI_ratio = ICI1/ICI2)

write.csv(age2_data, "age2.5_ICI.csv", na = "")


age5_data <- inner_join(age_5_celiac_ICI, age_5_control_ICI, by = c('ASV','asv_label','Taxa') , copy =T)
age5_data = age5_data  %>% mutate(ICI_ratio = ICI1/ICI2)

write.csv(age5_data, "age5_ICI.csv", na = "")
