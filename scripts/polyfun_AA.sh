
#The pipeline for PolyFun:

	#1. Create a munged summary statistics file in a PolyFun-friendly parquet format.

python ~/polyfun/munge_polyfun_sumstats.py \
--sumstats munged_files/${disease_type}.tsv \
--out gwas_studies/AD_Jansen2019_munged.polyfun.parquet 

#! Make sure that the 0s dont switch to +0 as these are read as characters in many tools including polyfun, ldsc 
#Can use options(scipen = 999) to get around +0


	#2. Creating your own annotations (see script scripts/functional_annotations_polyfun.Rmd)

#Annotation is a matrix which contains SNPs from GWAS that you are interested in and its presence/absence (binary 1/0) within each functional annotation + functional annotations #for ~19 million UK Biobank imputed SNPs with MAF>0.1%, based on the baseline-LF 2.2.UKB annotations (187 annotations) https://github.com/omerwe/polyfun/issues/198 - each #functional annotation is column 

#For each chromosome you need to create two or three files:

	#1. An annotations file containing columns for CHR, BP, SNP, A1, A2 and arbitrary other columns representing your annotations. These files can be either .parquet or .gz files (we recommend using .parquet files).
	#2. An file with extension .l2.M containing a single line with the sums of the columns of each annotation (whitespace-delimited).
	#3. (optional) A file with extension l2.M_5_50 that is similar to the .M file but considers only common SNPs (with MAF between 5% and 50%). By default PolyFun will not use #these files, but you can use them when running S-LDSC to estimate enrichment of common SNP heritability.

	#3. Computing LD-scores with pre-computed UK Biobank LD matrices

python ~/polyfun/compute_ldscores_from_ld.py \
  --annot annot_files/AD_Jansen2019_ADvas_AAA_20240422_hg19.${chrom}.annot.parquet \
  --ukb \
  --ld-dir UKBiobank_LD \
  --out annot_files/AD_Jansen2019_ADvas_AAA_20240422_hg19.${chrom}.l2.ldscore.parquet


	#4. Run PolyFun with L2-regularized S-LDSC

python ~/polyfun/polyfun.py \
    --compute-h2-L2 \
    --no-partitions \
    --output-prefix polyfun_output/AD_Jansen2019_ADvas_AAA_20240422_hg19 \
    --sumstats gwas_studies/polyfun_files/AD_Jansen2019_munged.polyfun.parquet \
    --ref-ld-chr annot_files/AD_Jansen2019_ADvas_AAA_20240422_hg19. \
    --w-ld-chr required_files/baselineLF_2.2.UKB_annot/baselineLF2.2.UKB/weights.UKB. \
    --allow-missing

	#5. Genome-wide fine-mapping step 1: Creating region-specific jobs

#--n using median of a gwas meta analysis as susie or finemap do not support per-SNP sample size -  https://github.com/omerwe/polyfun/issues/161

python ~/polyfun/create_finemapper_jobs.py \
    --sumstats polyfun_output/SVD_Sargurupremraj2022_ADvas_AAA_20240422_hg19.${chrom}.snpvar_ridge_constrained.gz \
    --n 48065 \
    --method susie \
    --max-num-causal 5 \
    --out-prefix SVD_Sargurupremraj2022_ADvas_AAA_20240422_hg19.${chrom} \
    --jobs-file create_finemapper_jobs/SVD_Sargurupremraj2022_ADvas_AAA_20240422_hg19.${chrom}.finemapper_jobs.txt


	#6. Genome-wide fine-mapping step 2: Invoking the region-specific jobs
#In this step the user should submit each of the commands listed in the output file of create_finemapper_jobs.py(specified via the argument --jobs-file) to their computer cluster.

	#7. Genome-wide fine-mapping step 3: Aggregating the results

python ~/polyfun/aggregate_finemapper_results.py \
 --out-prefix susie_output/SVD_Sargurupremraj2022_ADvas_AAA_20240422_hg19 \
 --sumstats  polyfun_output/SVD_Sargurupremraj2022_ADvas_AAA_20240422_hg19.1.snpvar_ridge_constrained.gz \
 --out aggregate_results/SVD_Sargurupremraj2022_ADvas_AAA_20240422_hg19.1.agg.txt.gz \
 --allow-missing-jobs 
done
