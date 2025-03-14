#!/bin/bash -l
#SBATCH -p chimp
#SBATCH -n 8
#SBATCH -t 5-00:00:00
#SBATCH -J tajimasD
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

cd /mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/vcf

for sample_vcf in *samples.vcf.filled 
do 

pop=${sample_vcf//_samples.vcf.filled/}
outFile=/mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/${pop}_mus_musculus_tajimasD

#/usr/local/sw/bcftools-1.18/bin/bcftools +fill-tags ${vcfFile} > ${vcfFile}.filled

vk tajima 10000 10000 ${sample_vcf} > ${outFile}

done