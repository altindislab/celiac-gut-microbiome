#########
#########
###  AGE 2 CONTROL
#########
#########
deneme1 =aggregate(comp_age2_pos_control$Abundance, list(comp_age2_pos_control$ASV), mean)
colnames(deneme1) <- c('ASV', 'Abundance1')

deneme2 =aggregate(comp_age2_neg_control$Abundance, list(comp_age2_neg_control$ASV), mean)
colnames(deneme2) <- c('ASV', 'Abundance2')

deneme3 = inner_join(deneme1 , deneme2 , by = 'ASV')

deneme5 = deneme3  %>% mutate(ICI2 = Abundance1/Abundance2)

deneme6 =inner_join(deneme5, annot_for_match, by = 'ASV')

deneme7 = select(deneme6, c('ASV','ICI2','Taxa'))

deneme7$asv_label <- paste(deneme7$ASV,deneme7$Taxa)

age_2_control_ICI <- deneme7

###NA list filter
deneme4 = merge(deneme1 , deneme2 , by = 'ASV', all=TRUE)
deneme4NA = filter(deneme4, is.na(deneme4$Abundance1) | is.na(deneme4$Abundance2))
age_2_control_NA_list <- deneme4NA

#########
#########
###  AGE 2 CELIAC
#########
#########

deneme1 =aggregate(comp_age2_pos_celiac$Abundance, list(comp_age2_pos_celiac$ASV), mean)
colnames(deneme1) <- c('ASV', 'Abundance1')

deneme2 =aggregate(comp_age2_neg_celiac$Abundance, list(comp_age2_neg_celiac$ASV), mean)
colnames(deneme2) <- c('ASV', 'Abundance2')

deneme3 = inner_join(deneme1 , deneme2 , by = 'ASV')
#deneme4 = merge(deneme1 , deneme2 , by = 'ASV', all=TRUE)

deneme5 = deneme3  %>% mutate(ICI1 = Abundance1/Abundance2)

deneme6 =inner_join(deneme5, annot_for_match, by = 'ASV')

deneme7 = select(deneme6, c('ASV','ICI1','Taxa'))

deneme7$asv_label <- paste(deneme7$ASV,deneme7$Taxa)

age_2_celiac_ICI <- deneme7

###NA list filter
deneme4 = merge(deneme1 , deneme2 , by = 'ASV', all=TRUE)
deneme4NA = filter(deneme4, is.na(deneme4$Abundance1) | is.na(deneme4$Abundance2))
age_2_celiac_NA_list <- deneme4NA

