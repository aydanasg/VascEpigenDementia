library(motifbreakR)
library(MotifDb)
library(TFBSTools)
library(universalmotif)
library(BSgenome)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(S4Vectors)
library(dplyr)
library(tidyverse)
library(jsonlite)
library(rtracklayer)

#Load and filter GWAS summary stats

GWAS <- read_tsv('GWAS.tsv')
GWAS <- GWAS[GWAS$P<5e-8,]
GWAS.mb <- snps.from.rsid(rsid = GWAS$SNP,
                         dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                         search.genome = BSgenome.Hsapiens.UCSC.hg19)

#Load transcription factor binding sites and subset SNPs with overlap
CnT <- import('CnT.bed')

GWAS_CnT <- subsetByOverlaps(GWAS.mb,CnT)

#Run motifbreakR and get p values

mb_results <- motifbreakR(
  snpList = GWAS_CnT,
  filterp = T,
  pwmList = subset(MotifDb, 
                   dataSource %in% c("HOCOMOCOv11-core-A", "HOCOMOCOv11-core-B", "HOCOMOCOv11-core-C","jaspar2022")),
  threshold = 1,
  method = "ic",
  bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
  BPPARAM = BiocParallel::bpparam("SerialParam")
)

mb_results <- calculatePvalue(mb_results)

mb_results <- as.data.frame(mb_results, row.names = 1:length(mb_results))

write_csv(mb_results,'mb_results.csv')