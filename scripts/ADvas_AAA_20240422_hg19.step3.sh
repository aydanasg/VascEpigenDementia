#!/bin/bash

#
#PBS -N ldsc_step12_ADvas_AAA_20240422_hg19
#PBS -o ldsc_step12__ADvas_AAA_20240422_hg19.log

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=5:mem=100gb

## Create conda env (DO THIS STEP BEFORE ATTEMPTING TO RUN LDSC) 
#conda env create --file environment.yml

module load anaconda3/personal
source activate ldsc

#paths to the directories 
echo "Creating paths to directories"
LDSC="/rds/general/user/aa19618/projects/epinott/live/scripts/ldsc/ldsc"
REQUIRED_FILES="/rds/general/user/aa19618/projects/epinott/live/scripts/ldsc/required_files"
ANNOT_FILES="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/ldsc/annot_file"
#Change the path in HERITABILITY_OUTPUT to the directory where your results will be stored 
HERITABILITY_OUTPUT="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/ldsc/heritability"
# Directory where the files are located
BEDFILES="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/nextflow_ct/analysed_data/${ID}/bedfiles/filtered_peaks"
CHAIN="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/ldsc/hg38ToHg19.over.chain.gz"

## Step 3: Running ldsc - disease enrichment  analysis 

#list your cell types for analysis EG:
#only names of cells should match their names in the filenames 
#list of diseases, only keep the disease you are interested in (select from the list below)

DISEASE="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/gwas_studies/ldsc_files/sumstats_ldsc"
export disease=("AD_Jansen2019" "AD_Kunkle2019" "PD_Nalls2019_proxy" "MS_Andlauer2016" "ALS_Rheenen2021" "SCZ_Trubetskoy2022" "CAD_Aragam2022" "AF_Nielsen2018" "DBP_Evangelou2018" "SBP_Evangelou2018" "Stroke_Malik2018" "Stroke_Mishra2022" "Longevity_Deelan2019_90th" "SVD_Sargurupremraj2022" "PVS_Duperron2023" "Covid_Castineira2023")

#For acetylaytion ADvas_AAA_20240422_ac

ID=ADvas_AAA_20240422_ac_hg19

regulation=('all')

for disease_type in "${disease[@]}"
do
  for regulation_type in "${regulation[@]}"
  do
  echo "Running ldsc for: ${disease_type}"
python ${LDSC}/ldsc.py \
--h2 ${DISEASE}/${disease_type}_ldsc.sumstats.gz \
--ref-ld-chr ${REQUIRED_FILES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ANNOT_FILES}/pu1_${ID}.,${ANNOT_FILES}/olig2_${ID}.,${ANNOT_FILES}/rfx4_${ID}.,${ANNOT_FILES}/neun_${ID}.,${ANNOT_FILES}/erg_${ID}.,${ANNOT_FILES}/notch3_${ID}. \
--out ${HERITABILITY_OUTPUT}/${ID}.${regulation_type}.${disease_type} \
--overlap-annot  \
--frqfile-chr ${REQUIRED_FILES}/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${REQUIRED_FILES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients
done
done

regulation=('promoters' 'enhancers')

for disease_type in "${disease[@]}"
do
  for regulation_type in "${regulation[@]}"
  do
  echo "Running ldsc for: ${disease_type}"
python ${LDSC}/ldsc.py \
--h2 ${DISEASE}/${disease_type}_ldsc.sumstats.gz \
--ref-ld-chr ${REQUIRED_FILES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ANNOT_FILES}/pu1_ac_hg19.${regulation_type}.,${ANNOT_FILES}/olig2_ac_hg19.${regulation_type}.,${ANNOT_FILES}/rfx4_ac_hg19.${regulation_type}.,${ANNOT_FILES}/neun_ac_hg19.${regulation_type}.,${ANNOT_FILES}/erg_ac_hg19.${regulation_type}.,${ANNOT_FILES}/notch3_ac_hg19.${regulation_type}. \
--out ${HERITABILITY_OUTPUT}/${ID}.${regulation_type}.${disease_type} \
--overlap-annot  \
--frqfile-chr ${REQUIRED_FILES}/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${REQUIRED_FILES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients
done
done

ID=ADvas_AAA_20240422_tf_hg19

regulation=('tf')

for disease_type in "${disease[@]}"
do
  for regulation_type in "${regulation[@]}"
  do
  echo "Running ldsc for: ${disease_type}"
python ${LDSC}/ldsc.py \
--h2 ${DISEASE}/${disease_type}_ldsc.sumstats.gz \
--ref-ld-chr ${REQUIRED_FILES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.,${ANNOT_FILES}/pu1_${ID}.,${ANNOT_FILES}/olig2_${ID}.,${ANNOT_FILES}/erg_${ID}. \
--out ${HERITABILITY_OUTPUT}/${ID}.${regulation_type}.${disease_type} \
--overlap-annot  \
--frqfile-chr ${REQUIRED_FILES}/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${REQUIRED_FILES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients
done
done

exit 0 


