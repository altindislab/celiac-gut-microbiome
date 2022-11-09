###Clean Version Of ASV ABUNDANCES
source("fcns/config.r")
#setwd("external/altindis/celiac_disease_16sn/phyloseq")
library(phyloseq)
library(vegan)
library(microbiome)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)


all_abundance <- (ps %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
                    #arrange(OTU) %>% rename(ASV = OTU) %>% 
                    select(OTU,Kingdom, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance))

all_abundance2 <- ps %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>%
  #arrange(OTU) %>% rename(ASV = OTU) %>% 
  select(OTU, Sample, Abundance)

flt <- filter(all_abundance, all_abundance$Abundance>0)

#asv_sequences from differential otu EDGER script which matches OTU sequences to ASV number
colnames(asv_sequences) <- c('ASV', 'OTU')
deneme3 <- inner_join(flt, asv_sequences, by ='OTU',copy=TRUE)

library(tibble)
md_dd <- tibble::rownames_to_column(md, "Sample")
mtt_dd <- tibble::rownames_to_column(mtt, "ASV")
annot_for_match <- tibble::rownames_to_column(annot_combined, "ASV")

asv_list <- mtt_dd$ASV

deneme4 <- inner_join(deneme3, md_dd, by ='Sample',copy=TRUE)

composite_data <-  deneme4

deneme5 <- inner_join(deneme4, mtt_dd, by ='ASV',copy=TRUE)

ultimate <- deneme5

deneme6 <- tibble::rownames_to_column(annot_combined, "ASV")

deneme7 <- inner_join(ultimate, deneme6, by ='ASV',copy=TRUE)

ultimate_data <- deneme7
