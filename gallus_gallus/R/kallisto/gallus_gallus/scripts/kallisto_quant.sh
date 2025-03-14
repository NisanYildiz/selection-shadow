#!/bin/bash -l
#SBATCH -p bonobo
#SBATCH -n 40
#SBATCH -t 5-00:00:00
#SBATCH -J kallisto_quant
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

direct=/mnt/NEOGENE1/projects/aneuploidy_2016/Second/scWGA/MA_2021/GSE114129

for SRR_name in `cat ${direct}/file_name.txt`

do 

mkdir /mnt/NEOGENE4/projects/melih_2020/kallisto/gallus_gallus/quant/${SRR_name}

/usr/local/sw/kallisto/build/src/kallisto quant --threads=40 \
-i /mnt/NEOGENE4/projects/melih_2020/kallisto/gallus_gallus/gallus_gallus.idx \
-o /mnt/NEOGENE4/projects/melih_2020/kallisto/gallus_gallus/quant/${SRR_name} \
${direct}/data/processed/trimmed/T"$SRR_name"_1.fastq ${direct}/data/processed/trimmed/T"$SRR_name"_2.fastq

done
