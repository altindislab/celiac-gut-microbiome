#knitr::opts_chunk$set(echo = FALSE, include = FALSE)
#setwd("B:/")
source("fcns/config.r")
source("fcns/select_top.r")
source("fcns/select_top10.r")
source("fcns/select_top15.r")
source("fcns/fcns_limma/edger_contrasts.r")
#setwd("external/altindis/celiac_disease_16sn/edger")
library(plyr)
library(reshape2)
library(edgeR)
library(ggpubr)
#```

## Purpose
#Differential OTUs by using edger [1]
#Results are in folder [edger](../edger).

## Data
#OTU table [otu_table.csv](../phyloseq/otu_table.csv), 
#Taxonomy table [taxonomy_table.csv](../phyloseq/taxonomy_table.csv), 
#Metadata table [metadata_table.csv](../phyloseq/metadata_table.csv)

#```{r parse}
counts <- data.matrix(read.csv("../CD_new/phyloseq/otu_table.csv", row.names = 1,  skipNul = TRUE, fileEncoding="UCS-2LE"))
pheno <- read.csv("../CD_new/phyloseq/metadata_table.csv", row.names = 1,  skipNul = TRUE, fileEncoding="UCS-2LE")
annot <- read.csv("../CD_new/phyloseq/taxonomy_table.csv", row.names = 1,  skipNul = TRUE, fileEncoding="UCS-2LE")

oneyeardata <- subset(pheno, Age=="1year")
otheryeardata <- subset(pheno, Age!="1year")
oneyearlist <- row.names(oneyeardata)
otheryearlist<- row.names(otheryeardata)


#counts2 <- counts %>% colnames(counts)!=(oneyearlist)
#excludecounts <- subset(counts, select= (oneyearlist))
#deneme= counts %>% semi_join(excludecounts, by=colnames())

excludecounts <- subset(counts, select= (otheryearlist))

counts <- excludecounts


##DD removed age 2.5 samples (19582)
pheno2 <- subset(pheno,Age=="2.5year")
pheno2_filtered <- subset(pheno2, ABIS.no!="ABIS19582")

#counts <- subset(counts, colnames!=oneyearlist)
#pheno <- subset(pheno, Age!="1year")


##DD removed age 5 samples (19582, 8743, 17974, 22068, 3933)
pheno5 <- subset(pheno,Age=="5year")
pheno5_filtered <- subset(pheno5, ABIS.no!="ABIS19582") ##No sample at this age already
pheno5_filtered <- subset(pheno5_filtered, ABIS.no!="ABIS8743")
pheno5_filtered <- subset(pheno5_filtered, ABIS.no!="ABIS17974")
pheno5_filtered <- subset(pheno5_filtered, ABIS.no!="ABIS22068")
pheno5_filtered <- subset(pheno5_filtered, ABIS.no!="ABIS3933")

##DD Combine filtered pheno2 and pheno5
pheno <- rbind(pheno2_filtered,pheno5_filtered)

##After review
ultimatelist<- row.names(pheno)
ultimatecounts <- subset(counts, select= (ultimatelist))

counts <- ultimatecounts
counts_yedek <- counts



#DD changed "OTU" to "ASV"
#rownames(counts) <- rownames(annot) <- paste0("OTU", 1:nrow(counts))
rownames(counts) <- rownames(annot) <- paste0("ASV", 1:nrow(counts))

##DD ADDITION FOR ICI PART 10.14.22
library(dplyr)
library(tibble)
composite <- counts
composite <- as.data.frame(composite)
counts <- as.data.frame(counts)

composite <- composite %>%
  rownames_to_column(var="ASV")


deneme<- row.names(composite) 
deneme <- mutate(composite, ASV=row.names(composite))

deneme <- cbind(rownames(counts),rownames(counts_yedek))

asv_sequences <-  deneme
#asv_sequences <- as.data.frame(asv_sequences)
#setnames(asv_sequences, c("ASV","OTU"))

#annot$Taxa <- apply(annot, 1, function(v) paste(v[!is.na(v)], collapse = "_"))
annot_species <- annot %>% filter(!is.na(Species))
annot_other <- annot %>% filter(is.na(Species))

annot_species$Taxa <- apply(annot_species,1, function(v) paste((tail(v[!is.na(v)],2) ),collapse = " "))
annot_other$Taxa <- apply(annot_other,1, function(v) paste((tail(v[!is.na(v)],1) ),collapse = " "))

annot_combined <- rbind(annot_species,annot_other)

annot <- annot_combined

#annot <- annot_species

#annot$Taxa <- apply(annot,1, function(v) paste((tail(v[!is.na(v)],2) ),collapse = " "))
pheno$Disease.status <- factor(pheno$Disease.status, levels = c("Control", "Celiac"))
pheno$Group <- factor(pheno$Group, levels = unique(pheno$Group))
```

#species_hit <- annot %>% filter(!is.na(Species))
#deneme <- annot %>% filter(!is.na(Species))
#annot$Taxanew <- apply(annot, 1, function(v) paste(v[!is.na(v)], collapse = " "))
#annot$deneme <- apply(annot, 1, function(v) paste(tail(v[!is.na(v)],2)),collapse = " ")

#annot$deneme123 <- apply(annot,1, function(v) (tail(v[!is.na(v)],2) )
#annot$deneme123 <- apply(annot,1, function(v) paste((tail(v[!is.na(v)],2) ),collapse = " "))

We combine the counts of each two technical replicates.   

## Filtering
```{r filt}
dge <- DGEList(counts = counts)
dge <- dge[rowSums(counts >= 10) > 5, ]
dim(dge)

# donot using  calcNormFactors()
```

To filter out low abundant OTUs, we keep OTUs that have read counts at least 10 in 5 samples. There are 611 OTUs after filtering. 

## Test for differential OUTs
```{r edger}
contr.v <- c(Celiac_vs_Control_in_Presort_Age2.5 = "Presort_Celiac_2.5 - Presort_Control_2.5",
             Celiac_vs_Control_in_IGneg_Age2.5 = "IGneg_Celiac_2.5 - IGneg_Control_2.5",
             Celiac_vs_Control_in_IGpos_Age2.5 = "IGpos_Celiac_2.5 - IGpos_Control_2.5", 
             IGpos_vs_IGneg_in_Control_Age2.5 = "IGpos_Control_2.5 - IGneg_Control_2.5",  
             IGpos_vs_IGneg_in_Celiac_Age2.5 = "IGpos_Celiac_2.5 - IGneg_Celiac_2.5",
             
             Celiac_vs_Control_in_Presort_Age5 = "Presort_Celiac_5 - Presort_Control_5",
             Celiac_vs_Control_in_IGneg_Age5 = "IGneg_Celiac_5 - IGneg_Control_5",
             Celiac_vs_Control_in_IGpos_Age5 = "IGpos_Celiac_5 - IGpos_Control_5", 
             IGpos_vs_IGneg_in_Control_Age5 = "IGpos_Control_5 - IGneg_Control_5",  
             IGpos_vs_IGneg_in_Celiac_Age5 = "IGpos_Celiac_5 - IGneg_Celiac_5")

pdf("edger_qc_metrics.pdf")
res <- edger.contrasts(dge, grp = pheno$Group, contrasts.v = contr.v, plot = TRUE)
dev.off()
mtt <- res$mtt; logcpm <- res$logcpm
rm(res)
mtt.df <- data.frame(signif(mtt, 3), annot[rownames(mtt), ])
source("fcns/fcns_plot/sig_hist.r")
sig.hist(mtt, name = "signif_hist_edgeR", pi0 = TRUE)
write.csv(mtt.df, "otu_stats_edgeR.csv", na = "")
#```

#To discover the differential OTUs, we use edgeR, an R package for differential expression analysis of digital gene expression data [1]. We perform empirical Bayes quasi-likelihood F-tests for the following comparisons: between any 2 disease status at different age and sorting, or between IGpos and IGneg in different age and disease status.    
#The histograms of significance are at [signif_hist.pdf](./signif_hist.pdf). If no OTUs were associated with the phenotype, we would expect the p-value histogram to be flat and all FDRs to be near one. The more associated OTUs there are, the more enrichment there is at low p-values, the lower will be the FDRs.  We also estimate the proportion of the true null hypothesis (i.e. non-significant OTUs) [2].    
#OTU statistics tables for all OTUs [otu_stats.csv](./otu_stats.csv). The table contains the average logCPM of each group, p-values, FDR, log fold-change, fold-change, and taxonomy information.      

## Plots
#```{r plots}
## boxplot
topotu1 <- select.top(mtt, contrasts.v = contr.v[grep("Presort", names(contr.v))])
pdf("top_otus_boxplot_presort.pdf", 9, 3)
for (otu in topotu1){
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort == "Presort"], pheno[pheno$Sort == "Presort", ])
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM)) + theme_bw()
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none")
  ggp <- ggp + geom_boxplot(mapping = aes(fill = Disease.status)) + facet_grid(~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test")
  plot(ggp)
}
dev.off()


#####
source("fcns/select_top.r")
pdf("top_otus_violin_plots_presort.pdf", 9, 3)
for (otu in topotu1){
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort == "Presort"], pheno[pheno$Sort == "Presort", ])
  #dat2p <- data.frame(logCPM = logcpm[otu, ], pheno)
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM, fill = Disease.status))+ theme_bw()
  #ggp <- ggp + ggtitle(annot[otu, "Taxa"]) + theme(legend.position = "none")
  
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none") ## DD ADDITION
  
  ggp <- ggp + geom_violin(trim = FALSE, alpha = 0.3) + geom_jitter(shape = 21, position = position_jitter(0.2), alpha = 0.5)
  ggp <- ggp + stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "black", width = 0.1) ##DD addition for median bar
  ggp <- ggp + facet_grid(~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test") 
  plot(ggp)
}
dev.off()

#####

#topotu2 <- select.top(mtt, contrasts.v = contr.v[-grep("Presort", names(contr.v))])
topotu2 <- select.top(mtt, contrasts.v = contr.v)
pdf("top_otus_boxplot_postsort.pdf", 9, 5)
for (otu in topotu2){
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort != "Presort"], pheno[pheno$Sort != "Presort", ])
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM)) + theme_bw()
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none")
  ggp <- ggp + geom_boxplot(mapping = aes(fill = Disease.status)) + facet_grid(Sort~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test")
  plot(ggp)
}
dev.off()

#####
source("fcns/select_top.r")
pdf("top_otus_violin_plots_postsort.pdf", 9, 3)
for (otu in topotu2){
  #dat2p <- data.frame(logCPM = logcpm[otu, ], pheno)
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort != "Presort"], pheno[pheno$Sort != "Presort", ])
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM, fill = Disease.status)) + theme_bw()
  #ggp <- ggp + ggtitle(annot[otu, "Taxa"]) + theme(legend.position = "none")
  
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none") ## DD ADDITION
  
  ggp <- ggp + geom_violin(trim = FALSE, alpha = 0.3) + geom_jitter(shape = 21, position = position_jitter(0.2), alpha = 0.5)
  
  #ggp <- ggp +  stat_summary(fun.data = "mean_sdl")
  ggp <- ggp + stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "black", width = 0.1)
  ggp <- ggp + facet_grid(~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test") 
  plot(ggp)
}
dev.off()

#####

library(ggstatsplot)
pdf("top_otus_violin_plots_postsort15.pdf", 9, 3)
for (otu in topotu2){
  #dat2p <- data.frame(logCPM = logcpm[otu, ], pheno)
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort != "Presort"], pheno[pheno$Sort != "Presort", ])
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM, fill = Disease.status)) + theme_bw()
  
  #ggp <- ggbetweenstats(data = dat2p, x = Disease.status, y = logCPM)
  #ggp <- ggp + ggtitle(annot[otu, "Taxa"]) + theme(legend.position = "none")
  
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none") ## DD ADDITION
  
  ggp <- ggp + geom_violin(trim = FALSE, alpha = 0.3) + geom_jitter(shape = 21, position = position_jitter(0.2), alpha = 0.5)+
    stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
                 colour = "black", width = 0.1)+stat_compare_means(method = "t.test")
  #ggo <- ggp + stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2 )
  ggp <- ggp + geom_jitter(shape = 21, position = position_jitter(0.2), alpha = 0.5)
  
  #ggp <- ggp +  stat_summary(fun.data = "mean_sdl")
  ggp <- ggp + facet_grid(~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test") 
  plot(ggp)
}
dev.off()

## boxplot for all
pdf("all_otus_boxplot_presort.pdf", 9, 3)
for (otu in rownames(logcpm)){
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort == "Presort"], pheno[pheno$Sort == "Presort", ])
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM)) + theme_bw()
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none")
  ggp <- ggp + geom_boxplot(mapping = aes(fill = Disease.status)) + facet_grid(~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test")
  plot(ggp)
}
dev.off()

#####
source("fcns/select_top.r")
pdf("all_otus_violin_plots_presort.pdf", 9, 3)
for (otu in  rownames(logcpm)){
  #dat2p <- data.frame(logCPM = logcpm[otu, ], pheno)
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort == "Presort"], pheno[pheno$Sort == "Presort", ])
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM, fill = Disease.status)) + theme_bw()
  #ggp <- ggp + ggtitle(annot[otu, "Taxa"]) + theme(legend.position = "none")
  
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none") ## DD ADDITION
  
  ggp <- ggp + geom_violin(trim = FALSE, alpha = 0.3) + geom_jitter(shape = 21, position = position_jitter(0.2), alpha = 0.5)
  ggp <- ggp + stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "black", width = 0.1) ##DD addition
  ggp <- ggp + facet_grid(~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test") 
  plot(ggp)
}
dev.off()

#####

pdf("all_otus_boxplot_postsort.pdf", 9, 5)
for (otu in rownames(logcpm)){
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort != "Presort"], pheno[pheno$Sort != "Presort", ])
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM)) + theme_bw()
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none")
  ggp <- ggp + geom_boxplot(mapping = aes(fill = Disease.status)) + facet_grid(Sort~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test")
  plot(ggp)
}
dev.off()

#####
source("fcns/select_top.r")
pdf("all_otus_violin_plots_postsort.pdf", 9, 3)
for (otu in  rownames(logcpm)){
  #dat2p <- data.frame(logCPM = logcpm[otu, ], pheno)
  dat2p <- data.frame(logCPM = logcpm[otu, pheno$Sort != "Presort"], pheno[pheno$Sort != "Presort", ])
  
  ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM, fill = Disease.status)) + theme_bw()
  #ggp <- ggp + ggtitle(annot[otu, "Taxa"]) + theme(legend.position = "none")
  
  ggp <- ggp + ggtitle(paste(otu, ",", annot[otu, "Taxa"])) + theme(legend.position = "none") ## DD ADDITION
  
  ggp <- ggp + geom_violin(trim = FALSE, alpha = 0.3) + geom_jitter(shape = 21, position = position_jitter(0.2), alpha = 0.5)
  ggp <- ggp + stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "black", width = 0.1) ## DD ADDITION
  ggp <- ggp + facet_grid(~ Age)+scale_fill_manual(values=c( "#00BFC4","#F8766D"))
  #stat_compare_means(label.x=2,label.y=0.5, method = "t.test") 
  plot(ggp)
}
dev.off()

#####

## Clostridium XlVa
otus <- readLines("Clostridium_XlVa.txt")
samp <- rownames(pheno)[pheno$Age == "5year"]
dat2p <- cbind(t(logcpm[otus, samp]), pheno[samp, c("Disease.status","Sort")])
dat2p <- melt(dat2p, id.vars =  c("Disease.status","Sort"),  variable.name = "OTU",value.name = "logCPM")

pdf("Clostridium_XlVa_boxplot.pdf", 6, 12)
ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM)) + theme_bw()
ggp <- ggp + ggtitle("Age: 5 years, Clostridium XlVa") + xlab("") + theme(legend.position = "none")
ggp <- ggp + geom_boxplot(mapping = aes(fill = Disease.status)) + facet_grid(OTU~Sort)
plot(ggp)
dev.off()

####pool 6
otus <- readLines("Clostridium_XlVa.txt")
samp <- rownames(pheno)[pheno$Age == "5year" & pheno$Sort == "Presort"]
dat2p <- cbind(logCPM = colMeans(logcpm[otus, samp]), pheno[samp, "Disease.status", drop = FALSE])
pval <- signif(t.test(logCPM~Disease.status, data = dat2p)$p.value, 3)

pdf("Clostridium_XlVa_boxplot2_presort.pdf", 3.5, 3.5)
ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM)) + theme_bw()
ggp <- ggp + ggtitle(paste("Age: 5 years, Presort, \n6 Clostridium XlVa, pval =", pval)) + xlab("") + theme(legend.position = "none")
ggp <- ggp + geom_boxplot(mapping = aes(fill = Disease.status))
plot(ggp)
dev.off()

###pool all
otus <- rownames(mtt.df)[which(mtt.df$Genus == "Clostridium_XlVa")]
samp <- rownames(pheno)[pheno$Age == "5year" & pheno$Sort == "Presort"]
dat2p <- cbind(logCPM = colMeans(logcpm[otus, samp]), pheno[samp, "Disease.status", drop = FALSE])
pval <- signif(t.test(logCPM~Disease.status, data = dat2p)$p.value, 3)

pdf("Clostridium_XlVa_boxplot3_presort.pdf", 3.5, 3.5)
ggp <- ggplot(data = dat2p, mapping = aes(x = Disease.status, y = logCPM)) + theme_bw()
ggp <- ggp + ggtitle(paste("Age: 5 years, Presort, \nAll Clostridium XlVa, pval =", pval)) + xlab("") + theme(legend.position = "none")
ggp <- ggp + geom_boxplot(mapping = aes(fill = Disease.status))
plot(ggp)
dev.off()

source("fcns/ezheat.R")
source("fcns/prune_mat.R")
## heatmaps
pheno.df1 <- pheno[pheno$Sort == "Presort", c("Disease.status", "Age")]
gaps.col1 <- which(diff(as.numeric(pheno.df1$Age), lag=1) != 0)
ezheat(logcpm[topotu1, pheno$Sort == "Presort"], pheno.df = pheno.df1,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu1, "Taxa"], labcols = "",
       gaps_col = gaps.col1, main = "logCPM", name = "top_otus_heat_presort", height = 8, width = 16, clip = 2)


pheno.df2 <- pheno[pheno$Sort != "Presort", c("Disease.status", "Age", "Sort")]
gaps.col2 <- which(diff(as.data.frame.numeric_version(pheno.df2$Age), lag=1) != 0)
ezheat(logcpm[topotu2, pheno$Sort != "Presort"], pheno.df = pheno.df2,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu2, "Taxa"], labcols = "",
       gaps_col = c(16,31,47,62,75,84,96), main = "logCPM", name = "top_otus_heat_postsort(10.3.22)", height = 8, width = 25, clip = 2)

pheno.df3 <- pheno[pheno$Sort != "Presort" & pheno$Age== "5year", c("Disease.status", "Sort")]
gaps.col3 <- which(diff(as.numeric(pheno.df3$Age), lag=1) != 0)
ezheat(logcpm[topotu2, pheno$Sort != "Presort"& pheno$Age== "5year"], pheno.df = pheno.df3,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu2, "Taxa"], labcols = "",
       gaps_col = gaps.col3, main = "logCPM", name = "top_otus_heat_postsort_5year(9.29.22)", height = 8, width = 12.5, clip = 2)

pheno.df4 <- pheno[pheno$Sort != "Presort" & pheno$Age== "2.5year", c("Disease.status", "Sort")]
gaps.col4 <- which(diff(as.numeric(pheno.df4$Age), lag=1) != 0)
ezheat(logcpm[topotu2, pheno$Sort != "Presort"& pheno$Age== "2.5year"], pheno.df = pheno.df4,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu2, "Taxa"], labcols = "",
       gaps_col = gaps.col4, main = "logCPM", name = "top_otus_heat_postsort_2_5year(7.8.22)", height = 8, width = 12.5, clip = 2)

pheno.df5 <- pheno[pheno$Sort == "Presort", c("Disease.status","Age")]
#gaps.col5 <- which(diff(as.numeric(pheno.df5$Age), lag=1) != 0)
ezheat(logcpm[topotu1, pheno$Age== "2.5year"], pheno.df = pheno.df5,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu1, "Taxa"], labcols = "",
       gaps_col = gaps.col5, main = "logCPM", name = "top_otus_heat_presort_2_5_year(7.8.22)", height = 8, width = 16, clip = 2)

pheno.df6 <- pheno[pheno$Sort == "Presort"& pheno$Age== "5year", c("Disease.status")]
gaps.col6 <- which(diff(as.numeric(pheno.df6$Age), lag=1) != 0)
ezheat(logcpm[topotu1, pheno$Sort == "Presort"& pheno$Age== "5year"], pheno.df = pheno.df6,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu1, "Taxa"], labcols = "",
       gaps_col = gaps.col6, main = "logCPM", name = "top_otus_heat_presort_5_year(7.8.22)", height = 8, width = 16, clip = 2)
```


###DD NEW HEATMAPS!
topotu7 <- select.top15(mtt, contrasts.v = contr.v[grep("Presort", names(contr.v))]) ###top30 otu presort

pheno.df7 <- pheno[pheno$Sort == "Presort", c("Age","Disease.status")]
gaps.col7 <- which(diff(as.numeric(pheno.df7$Age), lag=1) != 0)
ezheat(logcpm[topotu7, pheno$Sort == "Presort"], pheno.df = pheno.df7,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu7, "Taxa"], labcols = "",
       gaps_col = gaps.col7, main = "logCPM", name = "top15_otus_heat_presort", height = 8, width = 16, clip = 2)

Equal numbers of top OTUs (based on p-values) are selected from each comparison. The boxplots for top OUTs are at [top_otus_boxplot_presort.pdf](./top_otus_boxplot_presort.pdf), and [top_otus_boxplot_postsort.pdf](./top_otus_boxplot_postsort.pdf). The same sets of top OTUs are use in heatmaps [top_otus_heat_presort.pdf](./top_otus_heat_presort.pdf) and [top_otus_heat_postsort.pdf](./top_otus_heat_postsort.pdf).


###DD HEATMAPS AFTER REVIEW 10.11.22
## heatmaps
pheno.df1 <- pheno[pheno$Sort == "Presort", c("Disease.status", "Age")]
gaps.col1 <- which(diff(as.numeric(pheno.df1$Disease.status), lag=1) != 0)
ezheat(logcpm[topotu1, pheno$Sort == "Presort"], pheno.df = pheno.df1,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu1, "Taxa"], labcols = "",
       gaps_col = gaps.col1, main = "logCPM", name = "top_otus_heat_presort", height = 8, width = 16, clip = 2)


pheno.df2 <- pheno[pheno$Sort != "Presort", c("Disease.status", "Age", "Sort")]
gaps.col2 <- which(diff(as.numeric(pheno.df1$Disease.status), lag=1) != 0)
ezheat(logcpm[topotu2, pheno$Sort != "Presort"], pheno.df = pheno.df2,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu2, "Taxa"], labcols = "",
       gaps_col = gaps.col2, main = "logCPM", name = "top_otus_heat_postsort(10.11.22)", height = 8, width = 25, clip = 2)

pheno.df3 <- pheno[pheno$Sort != "Presort" & pheno$Age== "5year", c("Disease.status", "Sort")]
gaps.col3 <- which(diff(as.numeric(pheno.df3$Disease.status), lag=1) != 0)
ezheat(logcpm[topotu2, pheno$Sort != "Presort"& pheno$Age== "5year"], pheno.df = pheno.df3,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu2, "Taxa"], labcols = "",
       gaps_col = gaps.col3, main = "logCPM", name = "top_otus_heat_postsort_5year(10.11.22)", height = 8, width = 12.5, clip = 2)

pheno.df4 <- pheno[pheno$Sort != "Presort" & pheno$Age== "2.5year", c("Disease.status", "Sort")]
gaps.col4 <- which(diff(as.numeric(pheno.df4$Disease.status), lag=1) != 0)
ezheat(logcpm[topotu2, pheno$Sort != "Presort"& pheno$Age== "2.5year"], pheno.df = pheno.df4,  sc = "z", reorder_rows = TRUE, labrows = annot[topotu2, "Taxa"], labcols = "",
       gaps_col = gaps.col4, main = "logCPM", name = "top_otus_heat_postsort_2_5year(10.11.22)", height = 8, width = 12.5, clip = 2)

```{r check}
stopifnot(rownames(counts) == rownames(annot))
stopifnot(colnames(counts) == rownames(pheno))
stopifnot(mtt.df[1, "Taxa"] == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Citrobacter")
```

## Reference
[1] Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140.   
[2] Langaas, M, Ferkingstad, E, and Lindqvist, B (2005). Estimating the proportion of true null hypotheses, with application to DNA microarray data. Journal of the Royal Statistical Society Series B 67, 555-572.   

