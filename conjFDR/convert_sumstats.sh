
cd $save_folder
module load anaconda/Anaconda2
python $s_python_convert/sumstats.py csv \
--auto \
--sumstats $trait1 \
--n-val $trait1_n \
--out $(basename "$trait1" | sed "s/\..*$//").csv \
--force

python $s_python_convert/sumstats.py csv \
--auto \
--sumstats $trait2 \
--n-val $trait2_n \
--out $(basename "$trait2" | sed "s/\..*$//").csv \
--force

python $s_python_convert/sumstats.py mat \
--sumstats $(basename "$trait1" | sed "s/\..*$//").csv \
--ref /home/tengli/BQ_GWAS/pleiofdr-master/9545380.ref \
--out $(basename "$trait1" | sed "s/\..*$//").mat

python $s_python_convert/sumstats.py mat \
--sumstats $(basename "$trait2" | sed "s/\..*$//").csv \
--ref /home/tengli/BQ_GWAS/pleiofdr-master/9545380.ref \
--out $(basename "$trait2" | sed "s/\..*$//").mat
