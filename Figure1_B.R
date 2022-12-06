###Fig 1
#need to get stats (fishers tests presumably for some of these things)

library(EnhancedVolcano)
library(cowplot)
library(Genomation)
library(here)

source('general_scripts/general_functions.R')
load('Gene_expression_analysis/mm9_gene_expression.rdata.RData')

################################################


genes.gr <- conds.gr.list$dev
genes <- as.data.frame(genes.gr)


enhancers.h3k27ac <- import.bedGraph('ChIP_analysis/H3K27acc_diferentiation_310522.bed')
enhancers.h3k27ac$padj <- enhancers.h3k27ac$NA.6 ;  enhancers.h3k27ac$log2FoldChange <- enhancers.h3k27ac$NA.4

enhancers.h3k27ac <- enhancers.h3k27ac[-unique(queryHits(findOverlaps(enhancers.h3k27ac,genes.gr)))]
enhancers.h3k27ac[-unique(queryHits(findOverlaps(enhancers.h3k27ac,genes.gr,maxgap = 2000)))]

################################################
#A. Volcano plot of gene expression change
################################################


keyvals <- ifelse(
  genes$log2FoldChange < 0 & genes$padj < 0.001, 'royalblue',
  ifelse(genes$log2FoldChange > 0 & genes$padj < 0.001, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'black'] <- 'Unchanged'
names(keyvals)[keyvals == 'red'] <- 'Induced'
names(keyvals)[keyvals == 'royalblue'] <- 'Repressed'

A <- EnhancedVolcano(genes,
                     lab = genes$name,
                     x = 'log2FoldChange',
                     y = 'padj',
                     pointSize = 0.5,
                     labSize=0,
                     pCutoff=0.001,
                     FCcutoff=0,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE, 
                     cutoffLineType = 'blank',  
                     vline = 0,
                     vlineCol = 'black',
                     vlineType = 'dotted',
                     vlineWidth = 1.0,
                     title = 'RNA-seq DP to CD4 SP', colAlpha = 1,
                     selectLab = rownames(genes)[which(names(keyvals) %in% c('Repressed', 'Induced'))],
                     colCustom = keyvals,
                     subtitle = '',caption = '',titleLabSize = 0,
                     legendPosition = 'bottom',cutoffLineWidth = 0.1,
                     axisLabSize = 10,legendLabSize = 7,legendIconSize = 3
                     
                     
)




################################################
#A. Heatmatrices for comparing H3K27ac at differential regions DP --> SP
################################################

pal = c('#FFFFFF', brewer.pal(9,"Reds"))
windowsize <- 2000

DPchip1 <- '/ChIP_analysis/bams/14_K27acCD69negDPWTR1_mapped_sorted_RemoveDuplicates.bam'
DPchip2 <- '/ChIP_analysis/bams/14_K27acCD69negDPWTR2_mapped_sorted_RemoveDuplicates.bam'
SPchip2 <- '/ChIP_analysis/bams/ChIPseq008_K27acCD69posCD4SPWTR2_AGTCAA_mapped_sorted_RemoveDuplicates.bam'
SPchip1 <- '/ChIP_analysis/bams/ChIPseq007_K27acCD69posCD4SPWTR1_GATCAG_mapped_sorted_RemoveDuplicates.bam'

upenh <- enhancers.h3k27ac[enhancers.h3k27ac$padj < 0.001 & enhancers.h3k27ac$log2FoldChange >0]
downenh <- enhancers.h3k27ac[enhancers.h3k27ac$padj < 0.001 & enhancers.h3k27ac$log2FoldChange <0]
stabenh <- enhancers.h3k27ac[enhancers.h3k27ac$padj > 0.05]

#####
upsmDP1 <- ScoreMatrixBin(DPchip1,resize(upenh,windowsize,fix = 'center'),bin.num = windowsize/10)
downsmDP1 <- ScoreMatrixBin(DPchip1,resize(downenh,windowsize,fix = 'center'),bin.num = windowsize/10)
upsmDP2 <- ScoreMatrixBin(DPchip2,resize(upenh,windowsize,fix = 'center'),bin.num = windowsize/10)
downsmDP2 <- ScoreMatrixBin(DPchip2,resize(downenh,windowsize,fix = 'center'),bin.num = windowsize/10)

upsmSP1 <- ScoreMatrixBin(SPchip1,resize(upenh,windowsize,fix = 'center'),bin.num = windowsize/10)
downsmSP1 <- ScoreMatrixBin(SPchip1,resize(downenh,windowsize,fix = 'center'),bin.num = windowsize/10)
upsmSP2 <- ScoreMatrixBin(SPchip2,resize(upenh,windowsize,fix = 'center'),bin.num = windowsize/10)
downsmSP2 <- ScoreMatrixBin(SPchip2,resize(downenh,windowsize,fix = 'center'),bin.num = windowsize/10)


ups <- ScoreMatrixList(list(upsmDP1,upsmDP2,upsmSP2,upsmSP1))
downs <- ScoreMatrixList(list(downsmDP1,downsmDP2,downsmSP2,downsmSP1))


#get input values over same ranges

upsmDP1_C <- ScoreMatrixBin('/ChIP_analysis/bams/ChIPseq019_InputCD69negDPWTR1_GATCAG_mapped_sorted_RemoveDuplicates.bam',resize(upenh,windowsize,fix = 'center'),bin.num = windowsize/10)
downsmDP1_C <- ScoreMatrixBin('/ChIP_analysis/bams//ChIPseq019_InputCD69negDPWTR1_GATCAG_mapped_sorted_RemoveDuplicates.bam',resize(downenh,windowsize,fix = 'center'),bin.num = windowsize/10)
upsmDP2_C <- ScoreMatrixBin('/ChIP_analysis/bams//ChIPseq020_InputCD69negDPWTR2_GTGAAA_mapped_sorted_RemoveDuplicates.bam',resize(upenh,windowsize,fix = 'center'),bin.num = windowsize/10)
downsmDP2_C <- ScoreMatrixBin('/ChIP_analysis/bams/ChIPseq020_InputCD69negDPWTR2_GTGAAA_mapped_sorted_RemoveDuplicates.bam',resize(downenh,windowsize,fix = 'center'),bin.num = windowsize/10)

upsmSP2_C <- ScoreMatrixBin('/ChIP_analysis/bams/ChIPseq026_InputCD69posCD4SPWTR2_ATGTCA_mapped_sorted_RemoveDuplicates.bam',resize(upenh,windowsize,fix = 'center'),bin.num = windowsize/10)
downsmSP2_C <- ScoreMatrixBin('/ChIP_analysis/bams/ChIPseq026_InputCD69posCD4SPWTR2_ATGTCA_mapped_sorted_RemoveDuplicates.bam',resize(downenh,windowsize,fix = 'center'),bin.num = windowsize/10)
upsmSP1_C <- ScoreMatrixBin('/ChIP_analysis/bams/ChIPseq025_InputCD69posCD4SPWTR1_GTGGCC_mapped_sorted_RemoveDuplicates.bam',resize(upenh,windowsize,fix = 'center'),bin.num = windowsize/10)
downsmSP1_C <- ScoreMatrixBin('/ChIP_analysis/bams/ChIPseq025_InputCD69posCD4SPWTR1_GTGGCC_mapped_sorted_RemoveDuplicates.bam',resize(downenh,windowsize,fix = 'center'),bin.num = windowsize/10)

#normalise values to input
DP1em.UP <- enrichmentMatrix(upsmDP1,upsmDP1_C)
DP1em.Down <- enrichmentMatrix(downsmDP1,downsmDP1_C)
DP2em.UP <- enrichmentMatrix(upsmDP2,upsmDP2_C)
DP2em.Down <- enrichmentMatrix(downsmDP2,downsmDP2_C)

SP1em.UP <- enrichmentMatrix(upsmSP1,upsmSP1_C)
SP1em.Down <- enrichmentMatrix(downsmSP1,downsmSP1_C)
SP2em.UP <- enrichmentMatrix(upsmSP2,upsmSP2_C)
SP2em.Down <- enrichmentMatrix(downsmSP2,downsmSP2_C)

#clean lower end, setting a minimum value for plotting
DP1em.UP[DP1em.UP@.Data < 0] <- 0; DP1em.Down[DP1em.Down@.Data < 0] <- 0 ; DP2em.UP[DP2em.UP@.Data < 0] <- 0; DP2em.Down[DP2em.Down@.Data < 0] <- 0
SP1em.UP[SP1em.UP@.Data < 0] <- 0; SP1em.Down[SP1em.Down@.Data < 0] <- 0 ; SP2em.UP[SP2em.UP@.Data < 0] <- 0; SP2em.Down[SP2em.Down@.Data < 0] <- 0


ups <- ScoreMatrixList(list(DP1em.UP,DP2em.UP,SP1em.UP,SP2em.UP))

multiHeatMatrix(ups,col=pal,winsorize=c(0.1,99.9),
                matrix.main = c('DP R1','DP R2','SP R1','SP R2'),order = T,common.scale = T,xcoords = c(-windowsize/2000,windowsize/2000),column.scale = T)

downs <- ScoreMatrixList(list(DP1em.Down,DP2em.Down,SP1em.Down,SP2em.Down))

multiHeatMatrix(downs,col=pal,winsorize=c(0.1,99.9),
                matrix.main = c('DP R1','DP R2','SP R1','SP R2'),order = T,common.scale = T,xcoords = c(-windowsize/2000,windowsize/2000),column.scale = T)
