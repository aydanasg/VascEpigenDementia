#mtCOJO

ad=("AD_Bellenguez2022")
svd=('SVD_Sargurupremraj2022' 'DBP_Evangelou2018' 'SBP_Evangelou2018' 'Stroke_Malik2018' 'Stroke_Mishra2022' 'CAD_Aragam2022' 'AF_Nielsen2018' 'PVS_Duperron2023')

for ad_type in "${ad[@]}"
do
    for svd_type in "${svd[@]}"
    do
/rds/general/user/aa19618/home/mtCOJO/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 \
--mbfile 1000G_EUR_Phase3.mtcojo_ref_data.txt \
--mtcojo-file summary_list/${ad_type}_${svd_type}.mtCOJO_summary_data.list \
--ref-ld-chr LDscore/ \
--w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/ \
--out results/${ad_type}_${svd_type}.mtcojo_result
done
done