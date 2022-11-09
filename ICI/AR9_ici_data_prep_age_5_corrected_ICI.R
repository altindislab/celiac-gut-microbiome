###Compare All ICI with P values
library(data.table)
library(readr)
library(tidyverse)

age2_5_ICI <- read_csv("age2.5_ICI.csv")
age5_ICI <- read_csv("age5_ICI.csv")
df_to_parse <- mtt.df

###Add Row names of ASV as new column
df_to_parse=df_to_parse %>%  rownames_to_column(var="ASV")%>% remove_rownames 

###Prepare Data and Change Column Names
oldnms <- c("ICI1", "ICI2","ICI_ratio")
newnms <- c("ICI_Celiac", "ICI_Control","ICI_Celiac/ICI_Control")

age_5_ICI_renamed <- setnames(age5_ICI, oldnms, newnms)

df5_with_ICI <- inner_join(age_5_ICI_renamed, df_to_parse,  by="ASV")

####ICI Control List
###Filter ICI Control>10
df5_ICI_Control_Over10 <- filter(df5_with_ICI, df5_with_ICI$ICI_Control>5)

###Select Igpos vs Igneg comparison in Control

df5_ICI_filtered1 <- select(df5_ICI_Control_Over10, c(2,4,5,3,6,7,matches("IGpos_vs_IGneg_in_Control_Age5")))

#is.num <- sapply(df5_ICI_filtered1, is.numeric)
#df5_ICI_filtered1[is.num] <- lapply(df5_ICI_filtered1[is.num], round, 4)

df5_ICI_filtered2 <-filter(df5_ICI_filtered1, IGpos_vs_IGneg_in_Control_Age5.p < 0.05)

####ICI Celiac List
###Filter ICI Control>10
df5_ICI_Celiac_Over10 <- filter(df5_with_ICI, df5_with_ICI$ICI_Celiac>5)

###Select Igpos vs Igneg comparison in Control

df5_ICI_filtered3 <- select(df5_ICI_Celiac_Over10, c(2,4,5,3,6,7,matches("IGpos_vs_IGneg_in_Celiac_Age5")))

#is.num <- sapply(df5_ICI_filtered3, is.numeric)
#df5_ICI_filtered3[is.num] <- lapply(df5_ICI_filtered3[is.num], round, 4)

df5_ICI_filtered4 <-filter(df5_ICI_filtered3, IGpos_vs_IGneg_in_Celiac_Age5.p < 0.05)


###Rename Taxa.x to Taxa
df5_ICI_filtered2 <- setnames(df5_ICI_filtered2, "Taxa.x", "Taxa")
df5_ICI_filtered4 <- setnames(df5_ICI_filtered4, "Taxa.x", "Taxa")


age5_ICI_sig <- full_join(df5_ICI_filtered2,df5_ICI_filtered4)


write.csv(df5_ICI_filtered2, "age5_ICI5_pval_Control.csv")
write.csv(df5_ICI_filtered4, "age5_ICI5_pval_Celiac.csv")

