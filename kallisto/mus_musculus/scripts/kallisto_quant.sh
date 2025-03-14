#!/bin/bash -l
#SBATCH -p chimp
#SBATCH -n 40
#SBATCH -t 5-00:00:00
#SBATCH -J kallisto_index
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

for file in /mnt/NEOGENE1/projects/aneuploidy_2016/Second/scWGA/MA_2021/GSE99791/data/processed/trimmed/T*.fastq.gz;

do 

filename=$(basename "$file")

mkdir /mnt/NEOGENE4/projects/melih_2020/kallisto/mus_musculus/quant/${filename}

/usr/local/sw/kallisto/build/src/kallisto quant --threads=40 --single -l 50 -s 1 \
-i /mnt/NEOGENE4/projects/melih_2020/kallisto/mus_musculus/mus_musculus.idx \
-o /mnt/NEOGENE4/projects/melih_2020/kallisto/mus_musculus/quant/${filename} \
${file}


done

