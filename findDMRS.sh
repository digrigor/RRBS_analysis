#Bash script to use the combine-p software to perform DMR discovery on the DM analysis output.

#Load the python virtual environment.
source ~/software/combined-pvalues/venv/bin/activate
export PATH=$PATH:~/software/bedtools2/bin/
#Perform DMR discovery on the wilcoxon pvalues. python ~/software/combined-pvalues/cpv/comb-p pipeline -c 4 --dist 300 --step 60 --seed 0.05 --threshold 0.05 -p /store/LevelFour/ARTISTIC/METHYLATION/RRBS/dio01/pcomb_DMR_analysis/Output/wcox --annotate hg19 /store/LevelFour/ARTISTIC/METHYLATION/RRBS/dio01/pcomb_DMR_analysis/Input/bed_for_pcomb_wcox.bed


#Perform DMR discovery on the F-test with age and hpv type as covariates pvalues.
python ~/software/combined-pvalues/cpv/comb-p pipeline -c 4 --dist 300 --step 60 --seed 0.05 --threshold 0.05 -p /store/LevelFour/ARTISTIC/METHYLATION/RRBS/dio01/pcomb_DMR_analysis/Output/ftest --annotate hg19 /store/LevelFour/ARTISTIC/METHYLATION/RRBS/dio01/pcomb_DMR_analysis/Input/bed_for_pcomb_f.bed

#Deactivate the virtual environment.
deactivate
