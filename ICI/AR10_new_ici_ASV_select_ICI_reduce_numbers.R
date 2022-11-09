##08.09.22
####NEW ICI Select
####Take top for all comparison
library(dplyr)
age2_p <- read_csv("age2_p.csv")
age5_p <- read_csv("age5_p.csv")


####AGE 2 Filter
##filter out all ICI<1
filter1 <- filter(age2_p, age2_p$`Control Mean`>=1&age2_p$`Celiac Mean`>=1)
filter1sorted <- filter1[order(filter1$pvalue), ]
toptenflt1 = slice_head(filter1sorted, n = 10)

##ICI control < 1 and ICI Celiac >1
filter2 <- filter(age2_p, age2_p$`Celiac Mean`>=1&age2_p$`Control Mean`<=1)
filter2sorted <- filter2[order(filter2$pvalue), ]
toptenflt2 = slice_head(filter2sorted, n = 10)

##ICI control > 1 and ICI Celiac <1
filter3 <- filter(age2_p, age2_p$`Celiac Mean`<=1&age2_p$`Control Mean`>=1)
filter3sorted <- filter3[order(filter3$pvalue), ]
toptenflt3 = slice_head(filter3sorted, n = 10)

##ICI control > 10 or ICI Celiac >10
filter4 <- filter(age2_p, age2_p$`Celiac Mean`>=10|age2_p$`Control Mean`>=10)
filter4sorted <- filter4[order(filter4$pvalue), ]
toptenflt4 = slice_head(filter3sorted, n = 10)

##Porphyromonas asaccharolytica Addition
filter5 <- filter(age2_p, age2_p$Taxa=="Porphyromonas asaccharolytica")
toptenflt5 <- filter5
###Combine Filters
age2_5_top30 <- rbind(toptenflt1, toptenflt2, toptenflt3, toptenflt4, toptenflt5)

###Remove Duplicates
age2_5_top30 <- age2_5_top30 %>%distinct()

###Calculate Ratio
age2_5_top30 = mutate(age2_5_top30, ratio=`Control Mean`/`Celiac Mean`)

###Sort based on ratio
age2_5_top30 <- age2_5_top30[order(age2_5_top30$ratio), ]

#age2_row_order = as.factor(age2_5_top30$ASV)

###Select ASV and ratio for heatmap data order
ratio2_5 <- select(age2_5_top30, "ASV","ratio")



####AGE 5 Filter
##filter out all ICI<1
filter1 <- filter(age5_p, age5_p$`Control Mean`>=1&age5_p$`Celiac Mean`>=1)
filter1sorted <- filter1[order(filter1$pvalue), ]
toptenflt1 = slice_head(filter1sorted, n = 10)

##ICI control < 1 and ICI Celiac >1
filter2 <- filter(age5_p, age5_p$`Celiac Mean`>=1&age5_p$`Control Mean`<=1)
filter2sorted <- filter2[order(filter2$pvalue), ]
toptenflt2 = slice_head(filter2sorted, n = 10)

##ICI control > 1 and ICI Celiac <1
filter3 <- filter(age5_p, age5_p$`Celiac Mean`<=1&age5_p$`Control Mean`>=1)
filter3sorted <- filter3[order(filter3$pvalue), ]
toptenflt3 = slice_head(filter3sorted, n = 10)

##ICI control > 10 or ICI Celiac >10
filter4 <- filter(age5_p, age5_p$`Celiac Mean`>=10|age5_p$`Control Mean`>=10)
filter4sorted <- filter4[order(filter4$pvalue), ]
toptenflt4 = slice_head(filter3sorted, n = 10)

##Porphyromonas asaccharolytica Additiontaken from age 2 and modified to real values and remove ASV342
###Be careful age_5_ICI file might have wrong filename ICI1= Celiac ICI2= Control
filter5 <- filter(age2_p, age2_p$Taxa=="Porphyromonas asaccharolytica")
filter5$`Control Mean` = 0.06121536
filter5$`Celiac Mean` = 1.424882

toptenflt5 <- filter5
###Combine Filters
age5_top30 <- rbind(toptenflt1, toptenflt2, toptenflt3, toptenflt4, toptenflt5)

###Remove Duplicates and remove ASV342
age5_top30 <- age5_top30 %>%distinct()


###Calculate Ratio
age5_top30 = mutate(age5_top30, ratio=`Control Mean`/`Celiac Mean`)

###Sort based on ratio
age5_top30 <- age5_top30[order(age5_top30$ratio), ]

#age5_row_order = as.factor(age2_5_top30$ASV)

###Select ASV and ratio for heatmap data order
ratio5 <- select(age5_top30, "ASV","ratio")

