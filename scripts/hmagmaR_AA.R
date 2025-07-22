#Loading required packages
#remotes::install_github("aydanasg/hmagmaR", auth_token = "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", force = TRUE)
library(data.table)
library(hmagmaR)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#creating annotation files 

#creating annotation files 
cell=c("pu1", "neun",  "olig2", "erg", "notch3", "rfx4")

type<-c("pu1"="pu1", "erg"="endo_prim", "notch3"="peri_prim", "rfx4"="astro_prim", "neun" = "neun", "olig2" = "olig2")


#required files 
snpgeneexon<-fread("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/required_files/exonic_promoter_snps_chipSeeker.txt")
snps<-fread("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/required_files/g1000_snps.txt")
annotated_genes<-fread("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/required_files/annotated_genes_TxDb.Hsapiens.UCSC.hg19.knownGene.org.Hs.eg.db.txt")

#Creating AnnotationFileHmagma for cells and running GeneLevelAnalysis_hmagma

for(c in 1:length(cell)){

    print(paste(cell[c]))

    hic<-fread(paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/hic_data/hic_hg19/", type[[cell[c]]], "_interactome.hg19.bed"))
    hic<-hic[,c(1:6)]
    colnames(hic)<-c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
    
    promoterRegions<-get(paste0(cell[c], "_promoters"))
    enhancerRegions<-get(paste0(cell[c], "_enhancers"))
    
        
    hmagmaR::AnnotationFileHmagma(hic = hic, 
                                  promoterRegions = promoterRegions, 
                                  enhancerRegions = enhancerRegions,
                                  snps = snps, annotated_genes = annotated_genes, snpgeneexon = snpgeneexon, 
                                  AnnotationFile = paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/annotation_files/ADvas_AAA_20250305_pm_hg19/ADvas_AAA_20240422_20240910_20240918_20250305_resected_pm_hg19.", cell[c]))
}

#Running HMAGMA analysis 
cell=c("pu1", "neun",  "olig2", "erg", "notch3", "rfx4")
disease<-c('AD_Jansen2019', 'AD_Kunkle2019', 'AD_Bellenguez2022', 'Longevity_Deelan2019_90th', "PD_GP22025",
         'MS_Andlauer2016', 'ALS_Rheenen2021', 'SCZ_Trubetskoy2022',
         'DBP_Evangelou2018', 'SBP_Evangelou2018', 'StrokeAIS_Malik2018', 
         'CAD_Aragam2022', 'AF_Nielsen2018', 'SVD_Sargurupremraj2022', "PVS_Duperron2023", 
         "Covid_Castineira2023"


for(d in 1:length(disease)){
  for(c in 1:length(cell)){

    print(paste("ADvas_AAA_20250305_pm_hg19",disease[d], cell[c]))

    hmagmaR::GeneLevelAnalysis_hmagma(magma = "~/../home/HMAGMA/HMAGMA_system/magma", 
                        g1000 = "~/../home/HMAGMA_Protocol/required_files/g1000/g1000_files/g1000_eur", 
                        gwas = paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/gwas_studies/munged_files/", disease[d], "_munged.tsv"), AnnotationFile =  paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/annotation_files/ADvas_AAA_20250305_pm_hg19/ADvas_AAA_20240422_20240910_20240918_20250305_resected_pm_hg19.", cell[c], ".transcript.annot"), 
                        output = paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/magma_output/ADvas_AAA_20250305_pm_hg19/ADvas_AAA_20240422_20240910_20240918_20250305_resected_pm_hg19.", disease[d], ".", cell[c]))
}}
