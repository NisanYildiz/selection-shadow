#!/bin/bash -l
#SBATCH -p bonobo
#SBATCH -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J kallisto_index
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

cd /mnt/NEOGENE4/projects/melih_2020/kallisto/gallus_gallus

/usr/local/sw/kallisto/build/src/kallisto index -i gallus_gallus.idx --threads=16 \
/mnt/NEOGENE4/projects/melih_2020/db/transcriptome/Gallus_gallus.GRCg6a.cdna.all.fa.gz