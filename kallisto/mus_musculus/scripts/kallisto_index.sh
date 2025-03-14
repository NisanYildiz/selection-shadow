#!/bin/bash -l
#SBATCH -p chimp
#SBATCH -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J kallisto_index
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

cd /mnt/NEOGENE4/projects/melih_2020/kallisto/mus_musculus

/usr/local/sw/kallisto/build/src/kallisto index -i mus_musculus.idx --threads=16 \
/mnt/NEOGENE4/projects/melih_2020/db/transcriptome/Mus_musculus.GRCm38.cdna.all.fa.gz