## Step 1: Creating an annotation file for each cell type using 1) functional annotations in .bed format; and 2) bim file from 1000 Genome Phase 3 reference for european population 

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


## Step 2: Generating ldsc scores files for each type using 1) PLINK files for 1000 Genome Phase 3 reference for european population; 2) anotation files for each functional annotation in .annot.gz format; and 3) retaining only HapMap3 SNPs in linst.txt

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

## Step 3: Running ldsc - disease enrichment  analysis using 1) GWAS summary statistics formatted using munge_sumstats.py script; 2) Baseline model LD scores, standard regression weights, and allele frequencies which were based on the 1000 Genomes Phase 3 European population data; and 3)ld scores for our functional annotations

DISEASE="/rds/general/user/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/gwas_studies/ldsc_files/sumstats_ldsc"
export disease=("AD_Jansen2019" "AD_Kunkle2019" "PD_Nalls2019_proxy" "MS_Andlauer2016" "ALS_Rheenen2021" "SCZ_Trubetskoy2022" "CAD_Aragam2022" "AF_Nielsen2018" "DBP_Evangelou2018" "SBP_Evangelou2018" "Stroke_Malik2018" "Stroke_Mishra2022" "Longevity_Deelan2019_90th" "SVD_Sargurupremraj2022" "PVS_Duperron2023" "Covid_Castineira2023")

ID=ADvas_AAA_20240422_ac_hg19

regulation=('all' 'promoters' 'enhancers' 'tf')

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

exit 0 


