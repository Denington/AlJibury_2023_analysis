library(rtracklayer)
library(ggplot2)
library(ggsci)
library(GenomicInteractions)
library(genomation)
################
################

#cOff.U <- 0.001
#cOff.O <- 0.05
#mults <- T
#need to add the mults option to the standard measures, and add a distance thing
###############
###############

get_numbers <- function(df,feature){
  down.l <- length(df[df$bdif == 'Reduction' & df$Direction == 'Down' & df$feature == feature,]$dif) 
  down.s <- length(df[df$bdif == 'No change' & df$Direction == 'Down' & df$feature == feature,]$dif)
  down.g <- length(df[df$bdif == 'Increase' & df$Direction == 'Down' & df$feature == feature,]$dif)
  
  up.l <- length(df[df$bdif == 'Reduction' & df$Direction == 'Up' & df$feature == feature,]$dif) 
  up.s <- length(df[df$bdif == 'No change' & df$Direction == 'Up' & df$feature == feature,]$dif)
  up.g <- length(df[df$bdif == 'Increase' & df$Direction == 'Up' & df$feature == feature,]$dif)
  
  
  other.l <- length(df[df$bdif == 'Reduction' & df$Direction == 'Other' & df$feature == feature,]$dif) 
  other.s <- length(df[df$bdif == 'No change' & df$Direction == 'Other' & df$feature == feature,]$dif)
  other.g <- length(df[df$bdif == 'Increase' & df$Direction == 'Other' & df$feature == feature,]$dif)
  
  return(data.frame('Lost'=c(up.l,other.l,down.l),
                    'Neither'=c(up.s,other.s,down.s),
                    'Gained'=c(up.g,other.g,down.g),
                    'Feature_change' = c('Induced','No Change','Repressed'),row.names = c('Induced','No Change','Repressed')))
  
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


getGainStats <- function(df,method,feature){
  TP <- as.numeric(df['Induced','Gained'])
  FP <- as.numeric(df['No Change','Gained'] + df['Repressed','Gained'])
  TN <- as.numeric(df['No Change','Lost'] + df['No Change','Neither'] + df['Repressed','Lost'] + df['Repressed','Neither']) 
  FN <- as.numeric(df['Induced','Neither'] + df['Induced','Lost'])
  
  mccA <- ((TP * TN) - (FP * FN))
  mccB <- (TP + FP)*( TP + FN )*(TN + FP)*(TN + FN) 
  mcc <- round(mccA / sqrt(mccB),3)
  sens <- round(TP / (TP + FN),3)
  spec <- round(TN / (TN + FP),3)
  
  DOR <- round((TP*TN) / (FP*FN),3)
  Jac_Ind <- round(TP /(TP + FN + FP),3)
  
  return(c('Method'=method,'Feature'=feature,Direction='Up','Sensitivity'=sens,'Specificity'=spec,'MCC'=mcc,'Diagnostic_odds_ratio'=DOR,'Jaccard_Index'=Jac_Ind))
}

getLossStats <- function(df,method,feature){
  TP <- as.numeric(df['Repressed','Lost'])
  FP <- as.numeric(df['No Change','Lost'] + df['Induced','Lost'])
  TN <- as.numeric(df['No Change','Gained'] + df['No Change','Neither'] + df['Induced','Gained'] + df['Induced','Neither']) 
  FN <- as.numeric(df['Repressed','Neither'] + df['Repressed','Gained'])
  
  mccA <- ((TP * TN) - (FP * FN))
  mccB <- (TP + FP)*( TP + FN )*(TN + FP)*(TN + FN) 
  mcc <- round(mccA / sqrt(mccB),3)
  
  sens <- round(TP / (TP + FN),3)
  spec <- round(TN / (TN + FP),3)
  
  DOR <- round((TP*TN) / (FP*FN),3)
  Jac_Ind <- round(TP /(TP + FN + FP),3)
  return(c('Method'=method,'Feature'=feature,Direction='Down','Sensitivity'=sens,'Specificity'=spec,'MCC'=mcc,'Diagnostic_odds_ratio'=DOR,'Jaccard_Index'=Jac_Ind))
}

############

setwd('~/Downloads/Ins')

load('~/Downloads/hic_scripts-main/mm9_gene_expression.rdata.RData')
DP <- import.bw('CD69nDPWT.10kb.insul_score_100000.KRnorm.bw')
SP <- import.bw('CD69pSPWT.10kb.insul_score_100000.KRnorm.bw')

enh <- read.table('../H3K27acc_diferentiation_310522.bed')
enh.gr <- makeGRangesFromDataFrame(enh,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)

atac <- read.table('../ATACseq_diferentiation.bed')
atac.gr <- makeGRangesFromDataFrame(atac,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)


difexp_to_granges <- function(genechanges){
  changes <- genechanges #[wtwt$padj < 0.05,]
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


genes <- read.csv('~/Downloads/04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv')
genes.gr <- difexp_to_granges(genes)
genes.gr <- genes.gr[!is.na(genes.gr$log2FoldChange)]
genes.gr <- genes.gr[!is.na(genes.gr$padj)]


##############


DP <- DP[unique(queryHits(findOverlaps(DP,SP)))]
SP <- SP[unique(subjectHits(findOverlaps(DP,SP)))]

DP$score[is.na(DP$score)] <- 0
SP$score[is.na(SP$score)] <- 0


DP$dif <- SP$score - DP$score
dp <- DP; sp <- SP

down <- dp[(queryHits(findOverlaps(dp,enh.gr[enh.gr$V11 < cOff.U & enh.gr$V9 < 0])))]
up <- dp[(queryHits(findOverlaps(dp,enh.gr[enh.gr$V11 < cOff.U & enh.gr$V9 > 0])))]
other <- dp[(queryHits(findOverlaps(dp,enh.gr[enh.gr$V11 >= cOff.O])))]

up$Direction <- 'Up'
down$Direction <- 'Down'
other$Direction <- 'Other'

if(mults == F){
  other <- other[-queryHits(findOverlaps(other,down))] 
  other <- other[-queryHits(findOverlaps(other,up))]
}

woof <- rbind(as.data.frame(up),as.data.frame(down),as.data.frame(other)); woof$feature <- 'Enhancer'
#
down.a <- dp[(queryHits(findOverlaps(dp,atac.gr[atac.gr$V11 < cOff.U & atac.gr$V9 < 0])))]
up.a <- dp[(queryHits(findOverlaps(dp,atac.gr[atac.gr$V11 < cOff.U & atac.gr$V9 > 0])))]
other.a <- dp[(queryHits(findOverlaps(dp,atac.gr[atac.gr$V11 >= cOff.O])))]

up.a$Direction <- 'Up'
down.a$Direction <- 'Down'
other.a$Direction <- 'Other'

if(mults == F){
  other.a <- other.a[-queryHits(findOverlaps(other.a,down.a))] 
  other.a <- other.a[-queryHits(findOverlaps(other.a,up.a))]
}

woof.a <- rbind(as.data.frame(up.a),as.data.frame(down.a),as.data.frame(other.a)); woof.a$feature <- 'ATAC'

#
up.g <- dp[(queryHits(findOverlaps(dp,genes.gr[genes.gr$padj < cOff.U & genes.gr$log2FoldChange > 0])))]
down.g <- dp[(queryHits(findOverlaps(dp,genes.gr[genes.gr$padj < cOff.U & genes.gr$log2FoldChange < 0])))]
other.g <- dp[(queryHits(findOverlaps(dp,genes.gr[genes.gr$padj > cOff.O])))]

up.g$Direction <- 'Up'
down.g$Direction <- 'Down'
other.g$Direction <- 'Other'

if(mults == F){
  other.g <- other.g[-queryHits(findOverlaps(other.g,down.g))] 
  other.g <- other.g[-queryHits(findOverlaps(other.g,up.g))]
}

woof.g <- rbind(as.data.frame(up.g),as.data.frame(down.g),as.data.frame(other.g)); woof.g$feature <- 'Gene'

###

woof3 <- rbind(woof.g[,c('dif','Direction','feature')],woof[,c('dif','Direction','feature')],woof.a[,c('dif','Direction','feature')])

woof3$bdif <- 'No change'
woof3[woof3$dif < -0.2,]$bdif <- 'Reduction'
woof3[woof3$dif > 0.3,]$bdif <- 'Increase' #10th and 90th %iles
woof3$bdif <- factor(woof3$bdif,levels=c('Reduction','No change','Increase'))
woof3$difd <- 1

woof3$feature <- factor(woof3$feature,levels = c('Gene','Enhancer','ATAC'))


genes.ins <- get_numbers(woof3,'Gene')
enh.ins <- get_numbers(woof3,'Enhancer')
atac.ins <- get_numbers(woof3,'ATAC')


#####


genes.ins <- get_gainloss_cap(genes.ins,'Genes','Insulation')
enh.ins <- get_gainloss_cap(enh.ins,'Enhancers','Insulation')
atac.ins <- get_gainloss_cap(atac.ins,'ATAC','Insulation')

###

ins_stats <- rbind(
  getGainStats(genes.ins,'Insulation','Genes'),
  getLossStats(genes.ins,'Insulation','Genes'),
  
  getGainStats(enh.ins,'Insulation','Enhancers'),
  getLossStats(enh.ins,'Insulation','Enhancers'),
  
  getGainStats(atac.ins,'Insulation','ATAC'),
  getLossStats(atac.ins,'Insulation','ATAC')
)

##################### local con


setwd('~/Downloads/Local_Contacts')

#load('~/Downloads/hic_scripts-main/mm9_gene_expression.rdata.RData')
DP <- import.wig('CD69negDPWTR1.wig')
SP <- import.wig('CD69posCD4SPWTR1.wig')

DP <- DP[unique(queryHits(findOverlaps(DP,SP)))]
SP <- SP[unique(subjectHits(findOverlaps(DP,SP)))]
DP$score <- as.numeric(DP$score)

DP$score[is.na(DP$score)] <- 0
SP$score <- as.numeric(SP$NA.)
SP$score[is.na(SP$score)] <- 0


DP$dif <- SP$score - DP$score
dp <- DP; sp <- SP



down <- dp[(queryHits(findOverlaps(dp,enh.gr[enh.gr$V11 < cOff.U & enh.gr$V9 < 0])))]
up <- dp[(queryHits(findOverlaps(dp,enh.gr[enh.gr$V11 < cOff.U & enh.gr$V9 > 0])))]
other <- dp[(queryHits(findOverlaps(dp,enh.gr[enh.gr$V11 >= cOff.O])))]

up$Direction <- 'Up'
down$Direction <- 'Down'
other$Direction <- 'Other'

if(mults == F){
  other <- other[-queryHits(findOverlaps(other,down))] 
  other <- other[-queryHits(findOverlaps(other,up))]
}

woof <- rbind(as.data.frame(up),as.data.frame(down),as.data.frame(other)); woof$feature <- 'Enhancer'
#
down.a <- dp[(queryHits(findOverlaps(dp,atac.gr[atac.gr$V11 < cOff.U & atac.gr$V9 < 0])))]
up.a <- dp[(queryHits(findOverlaps(dp,atac.gr[atac.gr$V11 < cOff.U & atac.gr$V9 > 0])))]
other.a <- dp[(queryHits(findOverlaps(dp,atac.gr[atac.gr$V11 >= cOff.O])))]

up.a$Direction <- 'Up'
down.a$Direction <- 'Down'
other.a$Direction <- 'Other'

if(mults == F){
  other.a <- other.a[-queryHits(findOverlaps(other.a,down.a))] 
  other.a <- other.a[-queryHits(findOverlaps(other.a,up.a))]
}

woof.a <- rbind(as.data.frame(up.a),as.data.frame(down.a),as.data.frame(other.a)); woof.a$feature <- 'ATAC'

#
up.g <- dp[(queryHits(findOverlaps(dp,genes.gr[genes.gr$padj < cOff.U & genes.gr$log2FoldChange > 0])))]
down.g <- dp[(queryHits(findOverlaps(dp,genes.gr[genes.gr$padj < cOff.U & genes.gr$log2FoldChange < 0])))]
other.g <- dp[(queryHits(findOverlaps(dp,genes.gr[genes.gr$padj > cOff.O])))]

up.g$Direction <- 'Up'
down.g$Direction <- 'Down'
other.g$Direction <- 'Other'

if(mults == F){
  other.g <- other.g[-queryHits(findOverlaps(other.g,down.g))] 
  other.g <- other.g[-queryHits(findOverlaps(other.g,up.g))]
}

woof.g <- rbind(as.data.frame(up.g),as.data.frame(down.g),as.data.frame(other.g)); woof.g$feature <- 'Gene'


###

woof3 <- rbind(woof.g[,c('dif','Direction','feature')],woof[,c('dif','Direction','feature')],woof.a[,c('dif','Direction','feature')])

woof3$bdif <- 'No change'
woof3[woof3$dif < -0.075,]$bdif <- 'Reduction'
woof3[woof3$dif > 0.23,]$bdif <- 'Increase' #10th and 90th %iles
woof3$bdif <- factor(woof3$bdif,levels=c('Reduction','No change','Increase'))
woof3$difd <- 1


woof3$feature <- factor(woof3$feature,levels = c('Gene','Enhancer','ATAC'))


genes.locCon <- get_numbers(woof3,'Gene')
enh.locCon <- get_numbers(woof3,'Enhancer')
atac.locCon <- get_numbers(woof3,'ATAC')


genes.locCon <- get_gainloss_cap(genes.locCon,'Genes','Local Contacts')
enh.locCon <- get_gainloss_cap(enh.locCon,'Enhancers','Local Contacts')
atac.locCon <- get_gainloss_cap(atac.locCon,'ATAC','Local Contacts')

###

locc_stats <- rbind(
  getGainStats(genes.locCon,'Local Contacts','Genes'),
  getLossStats(genes.locCon,'Local Contacts','Genes'),
  
  getGainStats(enh.locCon,'Local Contacts','Enhancers'),
  getLossStats(enh.locCon,'Local Contacts','Enhancers'),
  
  getGainStats(atac.locCon,'Local Contacts','ATAC'),
  getLossStats(atac.locCon,'Local Contacts','ATAC')
)
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

setwd('~/Downloads/')

DPonlyCDs <- makeGRangesFromDataFrame(read.table('ContactDomains/DP.only CDs.bedpe',sep='\t',header=F) ,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3')
SPonlyCDs <- makeGRangesFromDataFrame(read.table('ContactDomains/SP.only.CDs.bedpe',sep='\t',header=F) ,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3')
sharedCDs <- makeGRangesFromDataFrame(read.table('ContactDomains/shared.CDs.bedpe',sep='\t',header=F) ,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3')


load('~/Downloads/hic_scripts-main/mm9_gene_expression.rdata.RData')
enh <- read.table('H3K27acc_diferentiation_310522.bed')
enh.gr <- makeGRangesFromDataFrame(enh,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
enh.atac <- read.table('ATACseq_diferentiation.bed')
atac.gr <- makeGRangesFromDataFrame(enh.atac,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)

atac.gr$padj <- atac.gr$V11
atac.gr$log2FoldChange <- atac.gr$V9

enh.gr$padj <- enh.gr$V11
enh.gr$log2FoldChange <- enh.gr$V9
#

get_numbers <- function(gained,lost,shared,feature,cutoffUnder=0.001,cutoffOver= 0.05){
  down.g <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj <= cutoffUnder & feature$log2FoldChange < 0 ],gained))]))
  up.g <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj <= cutoffUnder & feature$log2FoldChange > 0 ],gained))]))
  other.g <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj > cutoffOver  ],gained))]))
  
  down.l <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj <= cutoffUnder & feature$log2FoldChange < 0 ],lost))]))
  up.l <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj <= cutoffUnder & feature$log2FoldChange > 0 ],lost))]))
  other.l <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj > cutoffOver  ],lost))]))
  
  down.s <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj <= cutoffUnder & feature$log2FoldChange < 0 ],shared))]))
  up.s <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj <= cutoffUnder & feature$log2FoldChange > 0 ],shared))]))
  other.s <- length(unique(feature[queryHits(findOverlaps(feature[feature$padj > cutoffOver  ],shared))]))
  
  #return(c( up.l, length(feature[feature$padj <= cutoffUnder & feature$log2FoldChange > 0 ]) - up.g - up.l, up.g, 
  #          other.l, length(feature[feature$padj > cutoffOver ]) - other.g - other.l, other.g, 
  #          down.l, length(feature[feature$padj <= cutoffUnder & feature$log2FoldChange < 0 ]) - down.g - down.l, down.g)) 
  return(data.frame('Lost'=c(up.l,other.l,down.l),
                    'Neither'=c(
                      (length(feature[feature$padj <= cutoffUnder & feature$log2FoldChange > 0 ]) - up.g - up.l),
                      (length(feature[feature$padj > cutoffOver ]) - other.g - other.l),
                      (length(feature[feature$padj <= cutoffUnder & feature$log2FoldChange < 0 ]) - down.g - down.l)),
                    'Gained'=c(up.g,other.g,down.g),
                    'Feature_change' = c('Induced','No Change','Repressed'),row.names = c('Induced','No Change','Repressed')))
  
  #return(c( up.l, up.s, up.g, other.l, other.s, other.g, down.l, down.s, down.g)) 
}


genes.CDs <- get_numbers(SPonlyCDs,DPonlyCDs,sharedCDs,genes.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)
enh.CDs <- get_numbers(SPonlyCDs,DPonlyCDs,sharedCDs,enh.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)
atac.CDs <- get_numbers(SPonlyCDs,DPonlyCDs,sharedCDs,atac.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)



##########
#for tads, I'm saying they're *shared* if their mutual length overlap is >80% or so

DPtads <- import.bed('~/Downloads/Ya.DP.TADs.5kb.bed')
SPtads <- import.bed('~/Downloads/Ya.SP.TADs.5kb.bed')

hits <- findOverlaps(DPtads, SPtads,ignore.strand=T)
overlaps <- pintersect(DPtads[queryHits(hits)], SPtads[subjectHits(hits)],ignore.strand=T)
percentOverlap <- width(overlaps) / width(SPtads[subjectHits(hits)])
percentOverlap2 <- width(overlaps) / width(DPtads[queryHits(hits)])
overlaps$SPol <- percentOverlap
overlaps$DPol <- percentOverlap2
pO <- data.frame(percentOverlap,percentOverlap2)
overlaps$minO <- apply(pO,1,min)

DPtads.only <- overlaps[overlaps$minO < 0.6]

#mutualTADs <- hits[percentOverlap > 0.9 & percentOverlap2 > 0.9]
sharedTADs <- overlaps[overlaps$minO > 0.8]




#DPtads.only <- DPtads.only[-unique(queryHits(findOverlaps(DPtads.only,sharedTADs,minoverlap = 0.9,ignore.strand=T)))] 
#SPtads.only <- SPtads.only[-unique(queryHits(findOverlaps(SPtads.only,sharedTADs,minoverlap = 0.9,ignore.strand=T)))] 

DPtads.only <- unique(c(DPtads.only,DPtads[-queryHits(findOverlaps(DPtads,SPtads,ignore.strand=T))]))

#SP
hits <- findOverlaps(SPtads, DPtads,ignore.strand=T)
overlaps <- pintersect(SPtads[queryHits(hits)], DPtads[subjectHits(hits)],ignore.strand=T)
percentOverlap <- width(overlaps) / width(DPtads[subjectHits(hits)])
percentOverlap2 <- width(overlaps) / width(SPtads[queryHits(hits)])
overlaps$DPol <- percentOverlap
overlaps$SPol <- percentOverlap2
pO <- data.frame(percentOverlap,percentOverlap2)
overlaps$minO <- apply(pO,1,min)

SPtads.only <- overlaps[overlaps$minO < 0.6]


SPtads.only <- unique(c(SPtads.only,SPtads[-subjectHits(findOverlaps(DPtads,SPtads,ignore.strand=T))]))


genes.TADs <- get_numbers(SPtads.only,DPtads.only,sharedTADs,genes.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)
enh.TADs <- get_numbers(SPtads.only,DPtads.only,sharedTADs,enh.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)
atac.TADs <- get_numbers(SPtads.only,DPtads.only,sharedTADs,atac.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)

##########Loops

load.loops <- function(loopfile,header=T,seqrename=F,skipl=F){
  if(skipl == T){
    x <- read.table(loopfile, header=header,skip = 1)
  }
  else{
    x <- read.table(loopfile, header=header)}
  if(header == T){
    if(seqrename==T){
      x$chr1 <- paste0('chr',x$chr1)
      x$chr2 <- paste0('chr',x$chr2)
    }
    a <- makeGRangesFromDataFrame(x,seqnames.field = 'chr1',start.field = 'x1',end.field = 'x2')
    b <- makeGRangesFromDataFrame(x,seqnames.field = 'chr2',start.field = 'y1',end.field = 'y2')
  }
  else{
    if(seqrename==T){
      x$V1 <- paste0('chr',x$V1)
      x$V4 <- paste0('chr',x$V4)
    }
    a <- makeGRangesFromDataFrame(x,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3')
    b <- makeGRangesFromDataFrame(x,seqnames.field = 'V4',start.field = 'V5',end.field = 'V6')
  }
  #x <- pairs(GInteractions(a,b))
  #first(x) <- renameSeqlevels(first(x), paste0('chr',seqlevels(first(x))))
  #second(x) <- renameSeqlevels(second(x), paste0('chr',seqlevels(second(x))))
  #outloops <-  makeGInteractionsFromGRangesPairs(x)
  outloops <- GInteractions(a,b)
  return(outloops)
}

DPonlyLoops <- load.loops('~/Downloads/Loops/DP.difloops.bedpe',F,T) 
SPonlyLoops <- load.loops('~/Downloads/Loops/SP.difloops.bedpe',F,T) 

DPLoops <- load.loops('~/Downloads/Loops/DP.loops.bedpe',F,T,skipl = T) 
SPLoops <- load.loops('~/Downloads/Loops/SP.loops.bedpe',F,T,skipl = T) 

sharedLoops <- DPLoops[unique((linkOverlaps(DPLoops,SPLoops)$query))]

genes.loops <- get_numbers(SPonlyLoops,DPonlyLoops,sharedLoops,genes.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)
enh.loops <- get_numbers(SPonlyLoops,DPonlyLoops,sharedLoops,enh.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)
atac.loops <- get_numbers(SPonlyLoops,DPonlyLoops,sharedLoops,atac.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)


#######Subcompartments

subcompOlaps <- function(genes,A,B,mark, stab=F,cutoffUnder=0.001,cutoffOver = 0.05){
  both <- A
  both$cluster <- both$score
  both$cluster2 <- B$score
  both$A <- NaN ; both$B <- NaN
  cont = 1 # i know this is bad, i'm lazy 
  for(i in c('A1','A2','A3','B1','B2','B3')){
    both[both$cluster == i,]$A <- cont
    both[both$cluster2 == i,]$B <- cont
    cont <- cont + 1
  }
  both$dif <- both$A - both$B
  both$change <- 'None'
  both[both$dif < 0,]$change <- 'Decreased Activity'
  both[both$dif > 0,]$change <- 'Increased Activity'
  
  genes <- resize(genes,5000,'center')
  upgenes <- genes[genes$padj < cutoffUnder & genes$log2FoldChange > 0,]
  downgenes <- genes[genes$padj < cutoffUnder & genes$log2FoldChange < 0,]
  stabgenes <- genes[genes$padj > cutoffOver,]
  
  a <- both[unique(queryHits(findOverlaps(both,upgenes)))]
  b <- both[unique(queryHits(findOverlaps(both,downgenes)))]
  a$GeneChange <- 'Upregulated' 
  b$GeneChange <- 'Downregulated'  
  
  if(stab){
    c <- both[unique(queryHits(findOverlaps(both,stabgenes)))]
    c$GeneChange <- 'None'
    return(c(a,b,c))
  }
  #not <- genes[~(genes$ens_id %in% a$ens_id | genes$ens_id %in% b$ens_id | genes$ens_id %in% c$ens_id)]
  return(c(a,b))
  #return(data.frame(f=c(mark,mark,mark,mark),direction=c(direction,direction,direction,direction), Overlaps=c('CD69-','CD69+','both','neither'),n=c(length(a),length(b),length(c),length(not))))
}
#
DPsubs <- readBed('~/Downloads//testDPsubcomps.bed')
SPsubs <- readBed('~/Downloads//testSPsubcomps.bed')


subcomp.changes <- subcompOlaps(genes.gr,DPsubs,SPsubs,'woof')
#

genes.subc <- get_numbers(subcomp.changes[subcomp.changes$dif > 0],subcomp.changes[subcomp.changes$dif < 0],subcomp.changes[subcomp.changes$dif == 0],genes.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)
enh.subc <- get_numbers(subcomp.changes[subcomp.changes$dif > 0],subcomp.changes[subcomp.changes$dif < 0],subcomp.changes[subcomp.changes$dif == 0],enh.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)
atac.subc <- get_numbers(subcomp.changes[subcomp.changes$dif > 0],subcomp.changes[subcomp.changes$dif < 0],subcomp.changes[subcomp.changes$dif == 0],atac.gr,cutoffUnder = cOff.U,cutoffOver = cOff.O)


######################################################################################################################################################
######################################################################################################################################################

genes.CDs <- get_gainloss_cap(genes.CDs,'Genes','CDs')
enh.CDs <- get_gainloss_cap(enh.CDs,'Enhancers','CDs')
atac.CDs <- get_gainloss_cap(atac.CDs,'ATAC','CDs')

genes.TADs <- get_gainloss_cap(genes.TADs,'Genes','TADs')
enh.TADs <- get_gainloss_cap(enh.TADs,'Enhancers','TADs')
atac.TADs <- get_gainloss_cap(atac.TADs,'ATAC','TADs')

genes.loops <- get_gainloss_cap(genes.loops,'Genes','loops')
enh.loops <- get_gainloss_cap(enh.loops,'Enhancers','loops')
atac.loops <- get_gainloss_cap(atac.loops,'ATAC','loops')

genes.subc <- get_gainloss_cap(genes.subc,'Genes','subcompartments')
enh.subc <- get_gainloss_cap(enh.subc,'Enhancers','subcompartments')
atac.subc <- get_gainloss_cap(atac.subc,'ATAC','subcompartments')

###

standard_stats <- rbind(
  getGainStats(genes.CDs,'CDs','Genes'),
  getLossStats(genes.CDs,'CDs','Genes'),
  getGainStats(enh.CDs,'CDs','Enhancers'),
  getLossStats(enh.CDs,'CDs','Enhancers'),
  getGainStats(atac.CDs,'CDs','ATAC'),
  getLossStats(atac.CDs,'CDs','ATAC'),
  
  getGainStats(genes.TADs,'TADs','Genes'),
  getLossStats(genes.TADs,'TADs','Genes'),
  getGainStats(enh.TADs,'TADs','Enhancers'),
  getLossStats(enh.TADs,'TADs','Enhancers'),
  getGainStats(atac.TADs,'TADs','ATAC'),
  getLossStats(atac.TADs,'TADs','ATAC'),
  
  getGainStats(genes.loops,'loops','Genes'),
  getLossStats(genes.loops,'loops','Genes'),
  getGainStats(enh.loops,'loops','Enhancers'),
  getLossStats(enh.loops,'loops','Enhancers'),
  getGainStats(atac.loops,'loops','ATAC'),
  getLossStats(atac.loops,'loops','ATAC'),
  
  getGainStats(genes.subc,'subcompartments','Genes'),
  getLossStats(genes.subc,'subcompartments','Genes'),
  getGainStats(enh.subc,'subcompartments','Enhancers'),
  getLossStats(enh.subc,'subcompartments','Enhancers'),
  getGainStats(atac.subc,'subcompartments','ATAC'),
  getLossStats(atac.subc,'subcompartments','ATAC')
)  
  
outstats <- data.frame(rbind(ins_stats,locc_stats,standard_stats))
outstats$thresholds <- paste0(cOff.U, '_',cOff.O)


outtab <- rbind(
  genes.locCon,
  enh.locCon,
  atac.locCon,
  genes.ins,
  enh.ins,
  atac.ins,
  genes.CDs,
  enh.CDs,
  atac.CDs,
  genes.TADs,
  enh.TADs,
  atac.TADs,  
  genes.loops,
  enh.loops,
  atac.loops,
  genes.subc,
  enh.subc,
  atac.subc
)


#write.table(outstats,paste0('~/Desktop/feature_3Dchange_stats.', paste0(cOff.U, '_',cOff.O) ,'.tsv'),sep='\t',
#            quote = F,row.names = F)

#write.table(outtab,paste0('~/Desktop/feature_3Dchange_counts.', paste0(cOff.U, '_',cOff.O) ,'.tsv'),sep='\t',
#            quote = F,row.names = F)