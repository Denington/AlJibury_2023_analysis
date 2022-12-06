########################################################################################################################
#### Figure 5 stuff
#### Plotting how local contacts change with enhancer numbers/ k27ac atac differences
#for now, just loading 1 rep and comparing- will probs take a mean later
########################################################################################################################

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(ggpubr)
library(here)

cOff.U <- 0.001
cOff.O <- 0.05


########################################################################################################################################
###### functions ###################################################################################
########################################################################################################################################

get_numbers <- function(df,feature){
  down.l <- length(df[df$LC_change == 'Reduction' & df$change == 'Repressed' & df$Feature == feature,]$dif) 
  down.s <- length(df[df$LC_change == 'Unchanged' & df$change == 'Repressed' & df$Feature == feature,]$dif)
  down.g <- length(df[df$LC_change == 'Increase' & df$change == 'Repressed' & df$Feature == feature,]$dif)
  
  up.l <- length(df[df$LC_change == 'Reduction' & df$change == 'Induced' & df$Feature == feature,]$dif) 
  up.s <- length(df[df$LC_change == 'Unchanged' & df$change == 'Induced' & df$Feature == feature,]$dif)
  up.g <- length(df[df$LC_change == 'Increase' & df$change == 'Induced' & df$Feature == feature,]$dif)
  
  
  other.l <- length(df[df$LC_change == 'Reduction' & df$change == 'Unchanged' & df$Feature == feature,]$dif) 
  other.s <- length(df[df$LC_change == 'Unchanged' & df$change == 'Unchanged' & df$Feature == feature,]$dif)
  other.g <- length(df[df$LC_change == 'Increase' & df$change == 'Unchanged' & df$Feature == feature,]$dif)
  
  return(data.frame('Lost'=c(up.l,other.l,down.l),
                    'Neither'=c(up.s,other.s,down.s),
                    'Gained'=c(up.g,other.g,down.g),
                    'Feature_change' = c('Induced','Unchanged','Repressed'),row.names = c('Induced','Unchanged','Repressed')))
  
}



gaincap <- function(x){return(round(100*(as.numeric(x['Gained'])/ sum(c(as.numeric(x['Gained']),as.numeric(x['Neither']),as.numeric(x['Lost'])))),2))}
losscap <- function(x){return(round(100*(as.numeric(x['Lost'])/ sum(c(as.numeric(x['Gained']),as.numeric(x['Neither']),as.numeric(x['Lost'])))),2))}

get_gainloss_cap <- function(df,feature,method){
  df$GainedCaptured <- apply(df,1,gaincap)
  df$LostCaptured <- apply(df,1,losscap)
  
  df$Sensitivity 
  df$Feature <- feature
  df$Method <- method
  df$thresholds <- paste0(cOff.U, '_',cOff.O)
  
  return(df)
}

########################################################################################################################################
######loading/setting up local contacts ##############################################################
########################################################################################################################################

DP <- import.wig('Local_Contacts/CD69negDPWTR1.wig')
SP <- import.wig('Local_Contacts/CD69posCD4SPWTR1.wig')

DP <- DP[unique(queryHits(findOverlaps(DP,SP)))]
SP <- SP[unique(subjectHits(findOverlaps(DP,SP)))]
DP$score <- as.numeric(DP$score)

DP$score[is.na(DP$score)] <- 0
SP$score <- as.numeric(SP$NA.)
SP$score[is.na(SP$score)] <- 0


DP$dif <- SP$score - DP$score
DP$log2FoldDif <- log2(SP$score / DP$score)

DP$LC_change <- 'Unchanged'
DP[DP$dif < quantile(DP$dif,0.1),]$LC_change <- 'Reduction'
DP[DP$dif > quantile(DP$dif,0.9),]$LC_change <- 'Increase' #10th and 90th %iles
DP$LC_change <- factor(DP$LC_change,levels=c('Reduction','Unchanged','Increase'))

######loading enhancers/ genes / atac##############################################################

load('Gene_expression_analysis/mm9_gene_expression.rdata.RData')

enh <- read.table('ChIP_analysis/H3K27acc_diferentiation_310522.bed')
enh.gr <- makeGRangesFromDataFrame(enh,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
enh.gr$Feature_change <- 'Unchanged'; enh.gr$Feature_change[enh.gr$V11 < cOff.U & enh.gr$V9 > 0] <- 'Induced'; enh.gr$Feature_change[enh.gr$V11 < cOff.U & enh.gr$V9 < 0] <- 'Repressed'


atac <- read.table('ChIP_analysis/ATACseq_diferentiation.bed')
atac.gr <- makeGRangesFromDataFrame(atac,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
atac.gr$Feature_change <- 'Unchanged'; atac.gr$Feature_change[atac.gr$V11 < cOff.U & atac.gr$V9 > 0] <- 'Induced'; atac.gr$Feature_change[atac.gr$V11 < cOff.U & atac.gr$V9 < 0] <- 'Repressed'


difexp_to_granges <- function(genechanges){
  changes <- genechanges 
  changes.gr <- mm9_pc_genes.gr[mm9_pc_genes.gr$ens_id %in% genechanges$ensembl_gene_id ]
  changes.gr <- changes.gr[!duplicated(changes.gr$ens_id)]
  changes.gr <- promoters(changes.gr,upstream = 100,downstream = 100 ) #YES!
  
  changes <- changes[changes$ensembl_gene_id %in% changes.gr$ens_id,]
  changes <- changes[!duplicated(changes$ensembl_gene_id),]
  
  changes <- changes[order(changes$ensembl_gene_id),]
  changes.gr <- changes.gr[order(changes.gr$ens_id)]
  changes.gr <- promoters(changes.gr,upstream = 100,downstream = 100)
  
  changes.gr$log2FoldChange <- changes$log2FoldChange
  changes.gr$padj  <- changes$padj
  changes.gr$baseMean  <- changes$baseMean
  
  return(changes.gr)
}


genes <- read.csv('Gene_expression_analysis/04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv')
genes.gr <- difexp_to_granges(genes)
genes.gr <- genes.gr[!is.na(genes.gr$log2FoldChange)]
genes.gr <- genes.gr[!is.na(genes.gr$padj)]
genes.gr$Feature_change <- 'Unchanged'; genes.gr$Feature_change[genes.gr$padj < cOff.U & genes.gr$log2FoldChange > 0] <- 'Induced'; genes.gr$Feature_change[genes.gr$padj < cOff.U & genes.gr$log2FoldChange < 0] <- 'Repressed'

enh.gr <- enh.gr[-unique(queryHits(findOverlaps(enh.gr,genes.gr,maxgap = 2000)))]

#################################################################################################################################################
#################################################################################################################################################
####### A,B
#################################################################################################################################################

enh.up <- enh.gr[enh.gr$Feature_change == 'Induced']
enh.down <- enh.gr[enh.gr$Feature_change == 'Repressed']
enh.stab <- enh.gr[enh.gr$Feature_change == 'Unchanged']

atac.up <- atac.gr[atac.gr$Feature_change == 'Induced']
atac.down <- atac.gr[atac.gr$Feature_change == 'Repressed']
atac.stab <- atac.gr[atac.gr$Feature_change == 'Unchanged']

gene.up <- genes.gr[genes.gr$Feature_change == 'Induced']
gene.down <- genes.gr[genes.gr$Feature_change == 'Repressed']
gene.stab <- genes.gr[genes.gr$Feature_change == 'Unchanged']
#

gene.upLC <- DP[unique(queryHits(findOverlaps(DP,gene.up)))] ; gene.upLC$change <- 'Induced'  ; gene.upLC$Feature <- 'Genes'
gene.downLC <- DP[unique(queryHits(findOverlaps(DP,gene.down)))] ; gene.downLC$change <- 'Repressed'  ; gene.downLC$Feature <- 'Genes'
gene.stabLC <- DP[unique(queryHits(findOverlaps(DP,gene.stab)))]  ; gene.stabLC$change <- 'Unchanged' ; gene.stabLC$Feature <- 'Genes'   

enh.upLC <- DP[unique(queryHits(findOverlaps(DP,enh.up)))]   ; enh.upLC$change <- 'Induced'  ; enh.upLC$Feature <- 'Enhancers'
enh.downLC <- DP[unique(queryHits(findOverlaps(DP,enh.down)))]  ; enh.downLC$change <- 'Repressed' ; enh.downLC$Feature <- 'Enhancers'
enh.stabLC <- DP[unique(queryHits(findOverlaps(DP,enh.stab)))]  ; enh.stabLC$change <- 'Unchanged'  ; enh.stabLC$Feature <- 'Enhancers'   

atac.upLC <- DP[unique(queryHits(findOverlaps(DP,atac.up)))]  ; atac.upLC$change <- 'Induced'  ; atac.upLC$Feature <- 'ATAC'  
atac.downLC <- DP[unique(queryHits(findOverlaps(DP,atac.down)))]  ; atac.downLC$change <- 'Repressed' ; atac.downLC$Feature <- 'ATAC'
atac.stabLC <- DP[unique(queryHits(findOverlaps(DP,atac.stab)))] ; atac.stabLC$change <- 'Unchanged' ; atac.stabLC$Feature <- 'ATAC'

collectedLC.gr <- c(gene.upLC,gene.downLC,gene.stabLC,
                    enh.upLC,enh.downLC,enh.stabLC,
                    atac.upLC,atac.downLC,atac.stabLC)

collectedLC.df <- as.data.frame(collectedLC.gr)
collectedLC.df$change <- factor(collectedLC.df$change,c('Induced','Unchanged','Repressed'))
collectedLC.df$Feature <- factor(collectedLC.df$Feature,c('Enhancers','Genes','ATAC'))

get.wilcoxP <- function(df,fcA,fcB,clustA,clustB){
  
  return(wilcox.test(df$dif[df$Feature == fcA & df$change == clustA],
                     df$dif[df$Feature == fcB & df$change == clustB])$p.value)
  
}

pairwise.grouped <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~change,
  "Genes",   "Enhancers",     signif(get.wilcoxP(collectedLC.df,'Genes','Enhancers','Induced','Induced'),2), 1.0,"Induced",
  "Genes",   "Enhancers", signif(get.wilcoxP(collectedLC.df,'Genes','Enhancers','Repressed','Repressed'),2), 0.9,"Repressed",
  "Enhancers",   "Enhancers",     signif(get.wilcoxP(collectedLC.df,'Enhancers','Enhancers','Induced','Repressed'),2), 1.3,"Induced",
  "Genes",   "Genes", signif(get.wilcoxP(collectedLC.df,'Genes','Genes','Induced','Repressed'),2), 1.4,"Repressed")


ggplot(collectedLC.df[collectedLC.df$Feature != 'ATAC',],aes(Feature,dif,fill=change)) + geom_boxplot(outlier.alpha = 0.7) + theme_classic() + 
  ylab('Change in HiC contacts')  + scale_fill_manual(values=c("#EFC000FF","#868686FF","#0073C2FF")) + ggtitle('H3K27ac Induced/Unchanged/Repressed') +
  xlab('Enhancer fate DP --> SP') +   
  add_pvalue(pairwise.grouped, colour = "black", tip.length = 0.01,
    step.group.by = "change",
    step.increase = 0
  )

print('B')
Bout <- c(dim(collectedLC.df[collectedLC.df$Feature == 'Enhancers',])[1],
  dim(collectedLC.df[collectedLC.df$Feature == 'Enhancers' & collectedLC.df$change == 'Induced',])[1],
  dim(collectedLC.df[collectedLC.df$Feature == 'Enhancers' & collectedLC.df$change == 'Unchanged',])[1],
  dim(collectedLC.df[collectedLC.df$Feature == 'Enhancers' & collectedLC.df$change == 'Repressed',])[1],
  
  dim(collectedLC.df[collectedLC.df$Feature == 'Genes',])[1],
  dim(collectedLC.df[collectedLC.df$Feature == 'Genes' & collectedLC.df$change == 'Induced',])[1],
  dim(collectedLC.df[collectedLC.df$Feature == 'Genes' & collectedLC.df$change == 'Unchanged',])[1],
  dim(collectedLC.df[collectedLC.df$Feature == 'Genes' & collectedLC.df$change == 'Repressed',])[1]) #


pcdf <- rbind(
  get_gainloss_cap(get_numbers(collectedLC.df,'Enhancers'),feature = 'Enhancers',method = 'woof'),
  get_gainloss_cap(get_numbers(collectedLC.df,'Genes'),feature = 'Genes',method = 'woof'),
  get_gainloss_cap(get_numbers(collectedLC.df,'ATAC'),feature = 'ATAC',method = 'woof')
)

pcdf$Feature_change <- factor(pcdf$Feature_change,c('Induced','Unchanged','Repressed'))
pcdf$Feature <- factor(pcdf$Feature,c('Enhancers','Genes','ATAC'))
pcdf$lost_tofuse <- -1*pcdf$LostCaptured
pcdflong <- pcdf %>% pivot_longer(c(GainedCaptured,lost_tofuse),names_to = "DirectionContacts", values_to = "Proportion")

get.FisherstestU <- function(df,fcA,fcB,clustA,clustB){
  
  return(fisher.test(
    matrix(c(df$Gained[df$Feature_change == fcA & df$Feature == clustA], df$Neither[df$Feature_change == fcA & df$Feature == clustA] + df$Lost[df$Feature_change == fcA & df$Feature == clustA],
             df$Gained[df$Feature_change == fcB & df$Feature == clustB], df$Neither[df$Feature_change == fcB & df$Feature == clustB] + df$Lost[df$Feature_change == fcB & df$Feature == clustB]),nrow = 2)
  )$p.value)
  
}

get.FisherstestD <- function(df,fcA,fcB,clustA,clustB){
  
  return(fisher.test(
    matrix(c(df$Lost[df$Feature_change == fcA & df$Feature == clustA], df$Neither[df$Feature_change == fcA & df$Feature == clustA] + df$Gained[df$Feature_change == fcA & df$Feature == clustA],
             df$Lost[df$Feature_change == fcB & df$Feature == clustB], df$Neither[df$Feature_change == fcB & df$Feature == clustB] + df$Gained[df$Feature_change == fcB & df$Feature == clustB]),nrow = 2)
  )$p.value)
  
}

pairwise.grouped.I <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~Feature_change,
  'Enhancers','Genes',     signif(get.FisherstestU(pcdf,'Induced','Induced','Enhancers','Genes'),2), 55,   "Induced",
  'Enhancers','ATAC',    signif(get.FisherstestU(pcdf,'Induced','Induced','Enhancers','ATAC'),2), 60,   "Induced",
  'ATAC','Genes',    signif(get.FisherstestU(pcdf,'Induced','Induced','ATAC','Genes'),2), 55,   "Induced")

pairwise.grouped.R <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~Feature_change,
  'Enhancers','Genes',     signif(get.FisherstestD(pcdf,'Repressed','Repressed','Enhancers','Genes'),2), -65,  "Repressed",
  'Enhancers','ATAC',    signif(get.FisherstestD(pcdf,'Repressed','Repressed','Enhancers','ATAC'),2), -70,  "Repressed",
  'ATAC','Genes',    signif(get.FisherstestD(pcdf,'Repressed','Repressed','ATAC','Genes'),2), -65,  "Repressed"
)


ggplot(pcdflong,aes(Feature,Proportion,fill=Feature_change)) + 
  geom_bar(stat="identity",position = position_dodge(),color='black')  + theme_classic() + 
  scale_fill_manual(values=c("#EFC000FF","#868686FF","#0073C2FF")) + xlab('Number of enhancers') +
  add_pvalue(
    pairwise.grouped.I,
    colour = "black",
    tip.length = 0.01,
    step.group.by = "Feature_change",
    step.increase = 0
  ) +    add_pvalue(
    pairwise.grouped.R,
    colour = "black",
    tip.length = -0.01,
    step.group.by = "Feature_change",
    step.increase = 0
  ) +  ylab('Fraction of enhances with 3D change') #+

print('A')
Aout <- pcdf

sum(Aout[Aout$Feature == 'Enhancers',c('Lost','Neither','Gained')])
sum(Aout[Aout$Feature == 'Genes',c('Lost','Neither','Gained')])
sum(Aout[Aout$Feature == 'ATAC',c('Lost','Neither','Gained')])


########################################################################################################################################
########### (C) & (D))- splitting by numbers of enhancers                    ################################################
########################################################################################################################################

upClust <- makeGRangesFromDataFrame(read.table('~/Downloads/merged_induced_H3K27ac.50kb.bed',header = F),
                                    seqnames.field ='V1',start.field='V2',end.field = 'V3',keep.extra.columns = T)
downClust <- makeGRangesFromDataFrame(read.table('~/Downloads/merged_repressed_H3K27ac.50kb.bed',header = F),
                                      seqnames.field ='V1',start.field='V2',end.field = 'V3',keep.extra.columns = T)
stabClust <- makeGRangesFromDataFrame(read.table('~/Downloads/merged_unchanged_H3K27ac.50kb.bed',header = F),
                                    seqnames.field ='V1',start.field='V2',end.field = 'V3',keep.extra.columns = T)
stabClust <- stabClust[-unique(queryHits(findOverlaps(stabClust,c(upClust,downClust))))]

#
get_av_lcsO <- function(lc,clust,fc){
  
  lc$Feature_change <- fc
  
  c1 <- clust[clust$Count ==1] 
  c2 <- clust[clust$Count ==2]  
  c3 <- clust[clust$Count >=3]  
  
  lc1 <- lc[unique(queryHits(findOverlaps(lc,c1)))] ; lc1$NumEnhancers <- 1
  lc2 <- lc[unique(queryHits(findOverlaps(lc,c2)))] ; lc2$NumEnhancers <- 2
  lc3 <- lc[unique(queryHits(findOverlaps(lc,c3)))] ; lc3$NumEnhancers <- 3
  return(c(lc1,lc2,lc3))
}

get_av_lcs <- function(lc,clust,fc){
  
  lc$Feature_change <- fc
  
  c1 <- clust[clust$V4 =='Singlet'] 
  c2 <- clust[clust$V4 == 'Pairs']  
  c3 <- clust[clust$V4 =='Groups']  
  
  lc1 <- lc[unique(queryHits(findOverlaps(lc,c1)))] ; lc1$NumEnhancers <- 1 ; lc1$Enhancers <- 'Singlet'
  lc2 <- lc[unique(queryHits(findOverlaps(lc,c2)))] ; lc2$NumEnhancers <- 2 ; lc2$Enhancers <- 'Pairs'
  lc3 <- lc[unique(queryHits(findOverlaps(lc,c3)))] ; lc3$NumEnhancers <- 3 ; lc3$Enhancers <- 'Groups'
  return(c(lc1,lc2,lc3))
}


upClust.lc <- get_av_lcs(DP,upClust,'Induced')
downClust.lc <- get_av_lcs(DP,downClust,'Repressed')
stabClust.lc <- get_av_lcs(DP,stabClust,'Unchanged')

clusts <- as.data.frame(c(upClust.lc,downClust.lc,stabClust.lc))
clusts$Feature_change <- factor(clusts$Feature_change,c('Induced','Unchanged','Repressed'))
clusts$Enhancers <- factor(clusts$Enhancers,c('Singlet','Pairs','Groups'))
spare <- ggplot(clusts,aes(Feature_change,dif,fill=Feature_change)) + geom_boxplot(outlier.alpha = 0.7) + theme_classic() + 
  ylab('Difference in local contacts SP - DP') + scale_fill_manual(values=c("#EFC000FF","#868686FF","#0073C2FF")) + 
  facet_wrap(.~Enhancers) + stat_compare_means(comparisons = list(c('Induced','Repressed'))) + 
  xlab('Enhancer fate DP --> SP') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

get.wilcoxP <- function(df,fcA,fcB,clustA,clustB){

  return(wilcox.test(df$dif[df$Feature_change == fcA & df$Enhancers == clustA],
                  df$dif[df$Feature_change == fcB & df$Enhancers == clustB])$p.value)
  
}


pairwise.grouped.I <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~Feature_change,
  "Singlet",   "Pairs",     signif(get.wilcoxP(clusts,'Induced','Induced','Singlet','Pairs'),2), 1,        "Induced",
  "Singlet",   "Groups",    signif(get.wilcoxP(clusts,'Induced','Induced','Singlet','Groups'),2), 0.9,        "Induced",
  "Pairs",     "Groups",    signif(get.wilcoxP(clusts,'Induced','Induced','Pairs','Groups'),2), 1,        "Induced")

pairwise.grouped.R <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~Feature_change,
  "Singlet",   "Pairs",     signif(get.wilcoxP(clusts,'Repressed','Repressed','Singlet','Pairs'),2), -0.7,        "Repressed",
  "Singlet",   "Groups",    signif(get.wilcoxP(clusts,'Repressed','Repressed','Singlet','Groups'),2), -0.8,        "Repressed",
  "Pairs",     "Groups",    signif(get.wilcoxP(clusts,'Repressed','Repressed','Pairs','Groups'),2), -0.7,        "Repressed"
)

ggplot(clusts,aes(Enhancers,dif,fill=Feature_change)) + geom_boxplot(outlier.alpha = 0.7) + theme_classic() + 
  ylab('Difference in local contacts SP - DP') + scale_fill_manual(values=c("#EFC000FF","#868686FF","#0073C2FF")) + 
  xlab('Enhancer fate DP --> SP') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  add_pvalue(
    pairwise.grouped.I,
    colour = "black",
    tip.length = 0.01,
    step.group.by = "Feature_change",
    step.increase = 0
  ) +    add_pvalue(
    pairwise.grouped.R,
    colour = "black",
    tip.length = -0.01,
    step.group.by = "Feature_change",
    step.increase = 0
  )


print('E enhancer signlet/pair numbers')
Eout <- c(
  dim(clusts[clusts$Enhancers == 'Singlet',])[1],
  dim(clusts[clusts$Enhancers == 'Singlet' & clusts$Feature_change == 'Induced',])[1],
  dim(clusts[clusts$Enhancers == 'Singlet' & clusts$Feature_change == 'Unchanged',])[1],
  dim(clusts[clusts$Enhancers == 'Singlet' & clusts$Feature_change == 'Repressed',])[1],
  dim(clusts[clusts$Enhancers == 'Pairs',])[1],
  dim(clusts[clusts$Enhancers == 'Pairs' & clusts$Feature_change == 'Induced',])[1],
  dim(clusts[clusts$Enhancers == 'Pairs' & clusts$Feature_change == 'Unchanged',])[1],
  dim(clusts[clusts$Enhancers == 'Pairs' & clusts$Feature_change == 'Repressed',])[1],
  dim(clusts[clusts$Enhancers == 'Groups',])[1],
  dim(clusts[clusts$Enhancers == 'Groups' & clusts$Feature_change == 'Induced',])[1],
  dim(clusts[clusts$Enhancers == 'Groups' & clusts$Feature_change == 'Unchanged',])[1],
  dim(clusts[clusts$Enhancers == 'Groups' & clusts$Feature_change == 'Repressed',])[1]
)
  



#######
clusts$Feature <- as.character(clusts$Enhancers)
clusts$change <- clusts$Feature_change
pcdf <- rbind(
  get_gainloss_cap(get_numbers(clusts,'Singlet'),feature = 'Singlet',method = 'woof'),
  get_gainloss_cap(get_numbers(clusts,'Pairs'),feature = 'Pairs',method = 'woof'),
  get_gainloss_cap(get_numbers(clusts,'Groups'),feature = 'Groups',method = 'woof')
)

#pcdf[pcdf$Feature == 3,]$Feature <- '3+'
pcdf$Feature_change <- factor(pcdf$Feature_change,c('Induced','Unchanged','Repressed'))
pcdf$Feature <- factor(pcdf$Feature,c('Singlet','Pairs','Groups'))
pcdf$lost_tofuse <- -1*pcdf$LostCaptured
pcdflong <- pcdf %>% pivot_longer(c(GainedCaptured,lost_tofuse),names_to = "DirectionContacts", values_to = "Proportion")

get.FisherstestU <- function(df,fcA,fcB,clustA,clustB){
  
  return(fisher.test(
    matrix(c(df$Gained[df$Feature_change == fcA & df$Feature == clustA], df$Neither[df$Feature_change == fcA & df$Feature == clustA] + df$Lost[df$Feature_change == fcA & df$Feature == clustA],
             df$Gained[df$Feature_change == fcB & df$Feature == clustB], df$Neither[df$Feature_change == fcB & df$Feature == clustB] + df$Lost[df$Feature_change == fcB & df$Feature == clustB]),nrow = 2)
    )$p.value)
  
}

get.FisherstestD <- function(df,fcA,fcB,clustA,clustB){
  
  return(fisher.test(
    matrix(c(df$Lost[df$Feature_change == fcA & df$Feature == clustA], df$Neither[df$Feature_change == fcA & df$Feature == clustA] + df$Gained[df$Feature_change == fcA & df$Feature == clustA],
             df$Lost[df$Feature_change == fcB & df$Feature == clustB], df$Neither[df$Feature_change == fcB & df$Feature == clustB] + df$Gained[df$Feature_change == fcB & df$Feature == clustB]),nrow = 2)
  )$p.value)
  
}

pairwise.grouped.I <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~Feature_change,
  "Singlet",   "Pairs",     signif(get.FisherstestU(pcdf,'Induced','Induced','Singlet','Pairs'),2), 80,        "Induced",
  "Singlet",   "Groups",    signif(get.FisherstestU(pcdf,'Induced','Induced','Singlet','Groups'),2), 90,        "Induced",
  "Pairs",     "Groups",    signif(get.FisherstestU(pcdf,'Induced','Induced','Pairs','Groups'),2), 80,        "Induced")

pairwise.grouped.R <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~Feature_change,
  "Singlet",   "Pairs",     signif(get.FisherstestD(pcdf,'Repressed','Repressed','Singlet','Pairs'),2), -80,        "Repressed",
  "Singlet",   "Groups",    signif(get.FisherstestD(pcdf,'Repressed','Repressed','Singlet','Groups'),2), -90,        "Repressed",
  "Pairs",     "Groups",    signif(get.FisherstestD(pcdf,'Repressed','Repressed','Pairs','Groups'),2), -80,        "Repressed"
)


ggplot(pcdflong,aes(Feature,Proportion,fill=Feature_change)) + 
  geom_bar(stat="identity",position = position_dodge(),color='black')  + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values=c("#EFC000FF","#868686FF","#0073C2FF")) + xlab('Number of enhancers') +
  add_pvalue(
    pairwise.grouped.I,
    colour = "black",
    tip.length = 0.01,
    step.group.by = "Feature_change",
    step.increase = 0
  ) +    add_pvalue(
    pairwise.grouped.R,
    colour = "black",
    tip.length = -0.01,
    step.group.by = "Feature_change",
    step.increase = 0
  ) +  ylab('Fraction of enhances with 3D change') #+

print('D enhancer signlet/pair numbers')
Dout <- pcdf


########################################################################################################################################
########### (E) & (F)- splitting k27ac and atac depending on whether changing separately                    ################################################
########################################################################################################################################

max.dist <- 2000 #max distance between atac & k27ac peaks to be considered separate or not


enh.and.atac.up <- DP[unique(queryHits(findOverlaps(DP,enh.up[unique(queryHits(findOverlaps(enh.up,atac.up,maxgap = max.dist)))])))]  
enh.and.atac.up$multiD <- 'Up' ; enh.and.atac.up$Feature <- 'H3K27ac & ATAC Up'; enh.and.atac.up$change <- 'Induced'
enh.and.atac.up$otherChange <- 'Induced'

enh.and.atac.down <- DP[unique(queryHits(findOverlaps(DP,enh.down[unique(queryHits(findOverlaps(enh.down,atac.down,maxgap = max.dist)))])))]  
enh.and.atac.down$multiD <- 'Down' ; enh.and.atac.down$Feature <- 'H3K27ac & ATAC Down' ; enh.and.atac.down$change <- 'Repressed'
enh.and.atac.down$otherChange <- 'Repressed'

enh.only.up <- DP[unique(queryHits(findOverlaps(DP,enh.up[-unique(queryHits(findOverlaps(enh.up,c(atac.up,atac.down),maxgap = max.dist)))])))] 
enh.only.up$multiD <- 'Up'; enh.only.up$Feature <- 'H3K27ac only Up' ; enh.only.up$change <- 'Induced'
enh.only.up$otherChange <- 'Unchanged'

enh.only.down <- DP[unique(queryHits(findOverlaps(DP,enh.down[-unique(queryHits(findOverlaps(enh.down,c(atac.up,atac.down),maxgap = max.dist)))])))] 
enh.only.down$multiD <- 'Down' ; enh.only.down$Feature <- 'H3K27ac only Down' ; enh.only.down$change <- 'Repressed'
enh.only.down$otherChange <- 'Unchanged'

enhUp.atacDown <- DP[unique(queryHits(findOverlaps(DP,enh.up[unique(queryHits(findOverlaps(enh.up,atac.down,maxgap = max.dist)))])))]  
enhUp.atacDown$multiD <- 'Up' ; enhUp.atacDown$Feature <- 'H3K27ac Up, ATAC Down'; enhUp.atacDown$change <- 'Induced'
enhUp.atacDown$otherChange <- 'Repressed'

enhDown.atacUp <- DP[unique(queryHits(findOverlaps(DP,enh.down[unique(queryHits(findOverlaps(enh.down,atac.up,maxgap = max.dist)))])))]  
enhDown.atacUp$multiD <- 'Down' ; enhDown.atacUp$Feature <- 'H3K27ac Down, ATAC Up' ; enhDown.atacUp$change <- 'Repressed'
enhDown.atacUp$otherChange <- 'Induced'

#

atac.and.enh.up <- DP[unique(queryHits(findOverlaps(DP,atac.up[unique(subjectHits(findOverlaps(enh.up,atac.up,maxgap = max.dist)))])))] 
atac.and.enh.up$multiD <- 'Up' ; atac.and.enh.up$Feature <- 'H3K27ac & ATAC Up' ; atac.and.enh.up$change <- 'Induced'
atac.and.enh.up$otherChange <- 'Induced'

atac.and.enh.down <- DP[unique(queryHits(findOverlaps(DP,atac.down[unique(subjectHits(findOverlaps(enh.down,atac.down,maxgap = max.dist)))])))] 
atac.and.enh.down$multiD <- 'Down' ; atac.and.enh.down$Feature <- 'H3K27ac & ATAC Down' ; atac.and.enh.down$change <- 'Repressed'
atac.and.enh.down$otherChange <- 'Repressed'

atac.only.up <- DP[unique(queryHits(findOverlaps(DP,atac.up[-unique(subjectHits(findOverlaps(enh.up,c(atac.up,atac.down),maxgap = max.dist)))])))] 
atac.only.up$multiD <- 'Up' ; atac.only.up$Feature <- 'ATAC only Up' ; atac.only.up$change <- 'Induced'
atac.only.up$otherChange <- 'Unchanged'

atac.only.down <- DP[unique(queryHits(findOverlaps(DP,atac.down[-unique(subjectHits(findOverlaps(enh.down,c(atac.up,atac.down),maxgap = max.dist)))])))] 
atac.only.down$multiD <- 'Down' ; atac.only.down$Feature <- 'ATAC only Down' ; atac.only.down$change <- 'Repressed'
atac.only.down$otherChange <- 'Unchanged'

atacUp.enhDown <- DP[unique(queryHits(findOverlaps(DP,atac.up[unique(queryHits(findOverlaps(atac.up,enh.down,maxgap = max.dist)))])))]  
atacUp.enhDown$multiD <- 'Up' ; atacUp.enhDown$Feature <- 'H3K27ac Down, ATAC Up'; atacUp.enhDown$change <- 'Induced'
atacUp.enhDown$otherChange <- 'Repressed'

atacDown.enhUp <- DP[unique(queryHits(findOverlaps(DP,atac.down[unique(queryHits(findOverlaps(atac.down,enh.up,maxgap = max.dist)))])))]  
atacDown.enhUp$multiD <- 'Down' ; atacDown.enhUp$Feature <- 'H3K27ac Up, ATAC Down' ; atacDown.enhUp$change <- 'Repressed'
atacDown.enhUp$otherChange <- 'Induced'

#getting stab
atac.enh.stab <- DP[unique(queryHits(findOverlaps(DP,atac.stab[-unique(subjectHits(findOverlaps(c(enh.up,enh.down),atac.stab,maxgap = max.dist)))]
)))] 
atac.enh.stab$multiD <- 'Unchanged' ; atac.enh.stab$Feature <- 'Both Unchanged' ; atac.enh.stab$change <- 'Unchanged'
atac.enh.stab$otherChange <- 'Unchanged'

enh.atac.stab <- DP[unique(queryHits(findOverlaps(DP,enh.stab[-unique(subjectHits(findOverlaps(c(atac.up,atac.down),enh.stab,maxgap = max.dist)))]
)))] 
enh.atac.stab$multiD <- 'Unchanged' ; enh.atac.stab$Feature <- 'Both Unchanged' ; enh.atac.stab$change <- 'Unchanged'
enh.atac.stab$otherChange <- 'Unchanged'

atac.stab.enhUp <- DP[unique(queryHits(findOverlaps(DP,atac.stab[unique(subjectHits(findOverlaps(enh.up,atac.stab,maxgap = max.dist)))])))] 
atac.stab.enhUp$multiD <- 'Unchanged' ; atac.stab.enhUp$Feature <- 'H3K27ac Up, ATAC Unchanged' ; atac.stab.enhUp$change <- 'Unchanged'
atac.stab.enhUp$otherChange <- 'Induced'

atac.stab.enhDown <- DP[unique(queryHits(findOverlaps(DP,atac.stab[unique(subjectHits(findOverlaps(enh.down,atac.stab,maxgap = max.dist)))])))] 
atac.stab.enhDown$multiD <- 'Unchanged' ; atac.stab.enhDown$Feature <- 'H3K27ac Down, ATAC Unchanged' ; atac.stab.enhDown$change <- 'Unchanged'
atac.stab.enhDown$otherChange <- 'Repressed'

enh.stab.atacUp <- DP[unique(queryHits(findOverlaps(DP,enh.stab[unique(subjectHits(findOverlaps(atac.up,enh.stab,maxgap = max.dist)))])))] 
enh.stab.atacUp$multiD <- 'Unchanged' ; enh.stab.atacUp$Feature <- 'H3K27ac Unchanged, ATAC Up' ; enh.stab.atacUp$change <- 'Unchanged'
enh.stab.atacUp$otherChange <- 'Induced'

enh.stab.atacDown <- DP[unique(queryHits(findOverlaps(DP,enh.stab[unique(subjectHits(findOverlaps(atac.down,enh.stab,maxgap = max.dist)))])))] 
enh.stab.atacDown$multiD <- 'Unchanged' ; enh.stab.atacDown$Feature <- 'H3K27ac Unchanged, ATAC Down' ; enh.stab.atacDown$change <- 'Unchanged'
enh.stab.atacDown$otherChange <- 'Repressed'


#
atac.only.up$Feature <- 'ATAC only' ; atac.only.down$Feature <- 'ATAC only'
atac.stab.enhUp$Feature <- 'H3K27ac only' ; atac.stab.enhDown$Feature <- 'H3K27ac only'
atac.and.enh.up$Feature <- 'H3K27ac & ATAC' ; atac.and.enh.down$Feature <- 'H3K27ac & ATAC'

atac.stab.enhDown$change <- 'Repressed'; atac.stab.enhUp$change <- 'Induced'
#
enh.only.up$Feature <- 'H3K27ac only' ; enh.only.down$Feature <- 'H3K27ac only'
enh.stab.atacUp$Feature <- 'ATAC only' ; enh.stab.atacDown$Feature <- 'ATAC only'
enh.and.atac.up$Feature <- 'H3K27ac & ATAC' ; enh.and.atac.down$Feature <- 'H3K27ac & ATAC'

enh.stab.atacDown$change <- 'Repressed'; enh.stab.atacUp$change <- 'Induced'
####

#this is different as is both enhancer and atac centered, depending. probably the best
enh.atac.ss <- c(enh.only.up,enh.only.down,atac.only.up,atac.only.down,atac.and.enh.up,atac.and.enh.down)
enh.atac.ss <- as.data.frame(enh.atac.ss)
enh.atac.ss$change <- factor(enh.atac.ss$change,c('Induced','Repressed'))
enh.atac.ss$Feature <- factor(enh.atac.ss$Feature,c('H3K27ac only','ATAC only','H3K27ac & ATAC'))


pairwise.grouped.I <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~change,
  "H3K27ac only",   "ATAC only",     signif(get.wilcoxP(enh.atac.ss,'H3K27ac only','ATAC only','Induced','Induced'),2), 1,        "Induced",
  "H3K27ac only",   "H3K27ac & ATAC", signif(get.wilcoxP(enh.atac.ss,'H3K27ac only','H3K27ac & ATAC','Induced','Induced'),2), 0.9,        "Induced",
  "ATAC only",     "H3K27ac & ATAC", signif(get.wilcoxP(enh.atac.ss,'ATAC only','H3K27ac & ATAC','Induced','Induced'),2), 1,        "Induced")

pairwise.grouped.R <- tibble::tribble(
  ~group1, ~group2, ~p.adj,  ~y.position, ~change,
  "H3K27ac only",   "ATAC only",     signif(get.wilcoxP(enh.atac.ss,'H3K27ac only','ATAC only','Repressed','Repressed'),2), -1,        "Repressed",
  "H3K27ac only",   "H3K27ac & ATAC", signif(get.wilcoxP(enh.atac.ss,'H3K27ac only','H3K27ac & ATAC','Repressed','Repressed'),2), -0.9,        "Repressed",
  "ATAC only",     "H3K27ac & ATAC",  signif(get.wilcoxP(enh.atac.ss,'ATAC only','H3K27ac & ATAC','Repressed','Repressed'),2), -1,        "Repressed")

ggplot(enh.atac.ss,aes(Feature,dif,fill=change)) + geom_boxplot(outlier.alpha = 0.7) + theme_classic() + 
  ylab('Change in HiC contacts')  + ggtitle('H3K27ac Induced/Unchanged/Repressed') +
  xlab('Enhancer fate DP --> SP') + 
  scale_fill_manual(values=c("#EFC000FF","#0073C2FF")) +
  coord_cartesian(ylim = c(-1,1)) +   add_pvalue(
    pairwise.grouped.I,
    colour = "black",
    tip.length = 0.01,
    step.group.by = "change",
    step.increase = 0
  ) +   add_pvalue(
    pairwise.grouped.R,
    colour = "black",
    tip.length = -0.01,
    step.group.by = "change",
    step.increase = 0
  )

print('C')
Cout <- c(dim(enh.atac.ss[enh.atac.ss$Feature == 'H3K27ac only',])[1],
  dim(enh.atac.ss[enh.atac.ss$Feature == 'H3K27ac only' & enh.atac.ss$change == 'Induced',])[1],
  dim(enh.atac.ss[enh.atac.ss$Feature == 'H3K27ac only' & enh.atac.ss$change == 'Repressed',])[1],
  
  dim(enh.atac.ss[enh.atac.ss$Feature == 'ATAC only',])[1],
  dim(enh.atac.ss[enh.atac.ss$Feature == 'ATAC only' & enh.atac.ss$change == 'Induced',])[1],
  dim(enh.atac.ss[enh.atac.ss$Feature == 'ATAC only' & enh.atac.ss$change == 'Repressed',])[1],
  
  dim(enh.atac.ss[enh.atac.ss$Feature == 'H3K27ac & ATAC',])[1],
  dim(enh.atac.ss[enh.atac.ss$Feature == 'H3K27ac & ATAC' & enh.atac.ss$change == 'Induced',])[1],
  dim(enh.atac.ss[enh.atac.ss$Feature == 'H3K27ac & ATAC' & enh.atac.ss$change == 'Repressed',])[1])

