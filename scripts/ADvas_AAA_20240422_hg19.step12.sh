#!/bin/bash

#
#PBS -N ldsc_step12_ADvas_AAA_20240422_hg19
#PBS -o ldsc_step12__ADvas_AAA_2024422_hg19

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=20:mem=400gb

## Create conda env (DO THIS STEP BEFORE ATTEMPTING TO RUN LDSC) 
#conda env create --file environment.yml

module load anaconda3/personal
source activate ldsc

#making folders ]
echo "Making directories"
mkdir -p ~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/ldsc/annot_file
mkdir -p ~/../projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/ldsc/heritability

#paths to the directories 
echo "Creating paths to directories"
ID=ADvas_AAA_20240422_hg19

LDSC="/rds/general/user/aa19618/projects/epinott/live/scripts/ldsc/ldsc"
REQUIRED_FILES="/rds/general/user/aa19618/projects/epinott/live/scripts/ldsc/required_files"
ANNOT_FILES="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/ldsc/annot_file"
#Change the path in HERITABILITY_OUTPUT to the directory where your results will be stored 
HERITABILITY_OUTPUT="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/ldsc/heritability"
# Directory where the files are located
BEDFILES="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/ldsc/bedfiles/${ID}/all"

#Creating an array for sample_type and cell_type
echo "Creating an array for cell_type"

#An array for cell names 
cell=()

# Iterate through files in the specified directory
for file in "$BEDFILES/"*.bed; do
    if [ -e "$file" ]; then
        # Extract the base filename without the path and the trailing ".bed" pattern
        filename=$(basename "$file")
        cell+=("${filename%.bed}")
    fi
done

# Print the new array
for name in "${cell[@]}"; do
    echo "$name"
done

## Step 1: Creating an annot file 

for cell_type in "${cell[@]}"
do
    for chrom in {1..22}
    do
      echo "Creating an annot file for: ${cell_type} and chr ${chrom}"
## Step 1: Creating an annot file
  python ${LDSC}/make_annot.py \
  --bed-file ${BEDFILES}/${cell_type}.bed \
  --bimfile ${REQUIRED_FILES}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
  --annot-file ${ANNOT_FILES}/${cell_type}.${chrom}.annot.gz
done
done


## Step 2: Generating ldsc scores 

for cell_type in "${cell[@]}"
do
    for chrom in {1..22}
    do 
      echo "Generating ldsc scores for: ${cell_type} and chr ${chrom}"
python ${LDSC}/ldsc.py \
  --l2 \
  --bfile ${REQUIRED_FILES}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
  --ld-wind-cm 1 \
  --annot ${ANNOT_FILES}/${cell_type}.${chrom}.annot.gz  \
  --thin-annot \
  --out ${ANNOT_FILES}/${cell_type}.${chrom} \
  --print-snps ${REQUIRED_FILES}/list.txt 
done
done

exit 0 


