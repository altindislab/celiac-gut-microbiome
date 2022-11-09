#### Composite Data Split

comp_age2_pos <- filter(ultimate, ultimate$Age=='2.5year' & ultimate$Sort=='IGpos')
comp_age2_neg <- filter(ultimate, ultimate$Age=='2.5year' & ultimate$Sort=='IGneg')
comp_age5_pos <- filter(ultimate, ultimate$Age=='5year' & ultimate$Sort=='IGpos')
comp_age5_neg <- filter(ultimate, ultimate$Age=='5year' & ultimate$Sort=='IGneg')



comp_age2_pos_control <- filter(ultimate, ultimate$Age=='2.5year' & ultimate$Sort=='IGpos' & ultimate$Disease.status == 'Control')
comp_age2_neg_control <- filter(ultimate, ultimate$Age=='2.5year' & ultimate$Sort=='IGneg' & ultimate$Disease.status == 'Control')
comp_age2_pos_celiac <- filter(ultimate, ultimate$Age=='2.5year' & ultimate$Sort=='IGpos'  & ultimate$Disease.status == 'Celiac')
comp_age2_neg_celiac <- filter(ultimate, ultimate$Age=='2.5year' & ultimate$Sort=='IGneg'  & ultimate$Disease.status == 'Celiac')



comp_age5_pos_control <- filter(ultimate, ultimate$Age=='5year' & ultimate$Sort=='IGpos' & ultimate$Disease.status == 'Control')
comp_age5_neg_control <- filter(ultimate, ultimate$Age=='5year' & ultimate$Sort=='IGneg' & ultimate$Disease.status == 'Control')
comp_age5_pos_celiac <- filter(ultimate, ultimate$Age=='5year' & ultimate$Sort=='IGpos'  & ultimate$Disease.status == 'Celiac')
comp_age5_neg_celiac <- filter(ultimate, ultimate$Age=='5year' & ultimate$Sort=='IGneg'  & ultimate$Disease.status == 'Celiac')

#deneme =tapply(comp_age2_pos$Abundance, comp_age2_pos$ASV, mean)
#deneme1 =tapply(comp_age2_pos_control$Abundance, comp_age2_pos_control$ASV, mean)
#deneme2 =tapply(comp_age2_neg_control$Abundance, comp_age2_neg_control$ASV, mean)

deneme1 =aggregate(comp_age2_pos_control$Abundance, list(comp_age2_pos_control$ASV), mean)
colnames(deneme1) <- c('ASV', 'Abundance1')

deneme2 =aggregate(comp_age2_neg_control$Abundance, list(comp_age2_neg_control$ASV), mean)
colnames(deneme2) <- c('ASV', 'Abundance2')

deneme3 = inner_join(deneme1 , deneme2 , by = 'ASV')
deneme4 = merge(deneme1 , deneme2 , by = 'ASV', all=TRUE)

deneme5 = deneme3  %>% mutate(ICI = Abundance1/Abundance2)

deneme6 =inner_join(deneme5, annot_for_match, by = 'ASV')

deneme7 = select(deneme6, c('ASV','ICI','Taxa'))

deneme7$asv_label <- paste(deneme7$ASV,deneme7$Taxa)

age_2_control_ICI <- deneme7


