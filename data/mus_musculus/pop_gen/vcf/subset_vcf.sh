#!/bin/bash
#Subset vcf based on populations

cd /mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/vcf

for samplelist in *samplelist.txt 
do

pop=${samplelist//_samplelist.txt/}

/usr/local/sw/bcftools-1.18/bin/bcftools view -S ${samplelist} /mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/vcf/NA_East_GER_FRA_EXOME_merged_filtered_mm50_mq30_noIndels_noHWEoutliers_autosomes_LD50_pruned.2.16.22.vcf.filled -o ${pop}_samples.vcf.filled

done