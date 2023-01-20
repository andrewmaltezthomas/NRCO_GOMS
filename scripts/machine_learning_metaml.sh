# Download metaML 
git pull https://github.com/SegataLab/metaml

# Run LODO
declare -a datasets=("DerosaL_2020" "DerosaL_2022" "PRIMM_NL_LeeK_2022" "PRIMM_UK_LeeK_2022" "Manchester_LeeK_2022" "FrankelAE_2017" "GopalakrishnanV_2018" "McCullochJA_2022" "RoutyB_2022_Renal" "RoutyB_2022_Lung")
for i in "${datasets[@]}"
do
python classification_thomas-manghi.py MetaML_ORR_PDvsSDPRCR.txt ${i}_ORR_PDvsSDPRCR_rf -d 1:ORR:R -t Cohort:${i} -nt 1000 -nc 10 -nsl 5 -mf auto -hv 1 -df -r 1 -p 1 -z t__
done

# Run cross_validation
python classification_thomas-manghi.py MetaML_ORR_PDvsSDPRCR.txt CV_ORR_PDvsSDPRCR_rf -d 1:ORR:R -z t__ -nt 1000 -nc 10 -nsl 5 -mf auto -hv 1 -df -r 10 -p 5