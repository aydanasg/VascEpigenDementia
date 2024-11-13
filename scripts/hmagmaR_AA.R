#Loading required packages
#remotes::install_github("aydanasg/hmagmaR", auth_token = "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", force = TRUE)
library(data.table)
library(hmagmaR)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#creating annotation files 

cell<-c("pu1", "erg", "neun", "olig2", "notch3", "rfx4")
disease<-c('AD_Jansen2019', 'AD_Kunkle2019', 'AD_Bellenguez2022', 
         'PD_Nalls2019_proxy', 'MS_Andlauer2016', 'ALS_Rheenen2021', 'SCZ_Trubetskoy2022',
         'DBP_Evangelou2018', 'SBP_Evangelou2018', 'Stroke_Malik2018', 'Stroke_Mishra2022', 
         'CAD_Aragam2022', 'AF_Nielsen2018', 'SVD_Sargurupremraj2022')

#required files 
snpgeneexon<-fread("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/required_files/exonic_promoter_snps_chipSeeker.txt")
snps<-fread("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/required_files/g1000_snps.txt")
annotated_genes<-fread("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/required_files/annotated_genes_TxDb.Hsapiens.UCSC.hg19.knownGene.org.Hs.eg.db.txt")

#Creating AnnotationFileHmagma for cells and running hmagmaR::GeneLevelAnalysis_hmagma

for(c in 1:length(cell)){

    print(paste("ADvas_AAA_20240422_ac_hg19", cell[c]))

    hic<-fread(paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/hic_data/hic_hg19/", cell[c], "_interactome.hg19.bed"))
    regulatoryRegions<-fread(paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/nextflow_ct/analysed_data/ADvas_AAA_20240422_ac_hg19/bedfiles/filtered_peaks_afterQC/", cell[c], "_ADvas_AAA_20240422_ac_hg19.bed"))

    hmagmaR::AnnotationFileHmagma(hic = hic, 
                                  regulatoryRegions = regulatoryRegions, 
                                  snps = snps, annotated_genes = annotated_genes, snpgeneexon = snpgeneexon, 
                                  AnnotationFile = paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/annotation_files/ADvas_AAA_20240422_ac_hg19/ADvas_AAA_20240422_ac_hg19.", cell[c]))
}

for(d in 1:length(disease)){
  for(c in 1:length(cell)){

    print(paste("ADvas_AAA_20240422_ac_hg19",disease[d], cell[c]))

    hmagmaR::GeneLevelAnalysis_hmagma(magma = "~/../home/HMAGMA/HMAGMA_system/magma", 
                        g1000 = "~/../home/HMAGMA_Protocol/required_files/g1000/g1000_files/g1000_eur", 
                        gwas = paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/gwas_studies/munged_files/", disease[d], "_munged.tsv"), 
                        AnnotationFile = paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/annotation_files/ADvas_AAA_20240422_ac_hg19/ADvas_AAA_20240422_ac_hg19.", cell[c], ".transcript.annot"), 
                        output = paste0("~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/magma_output/ADvas_AAA_20240422_ac_hg19.", cell[c]))
}}
