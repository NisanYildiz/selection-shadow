#!/bin/bash -l
#SBATCH -p chimp
#SBATCH -n 8
#SBATCH -t 5-00:00:00
#SBATCH -J bedmap
#SBATCH -o slurm-%j-%N-%u.out
#SBATCH -e slurm-%J-%N-%u.err

cd /mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen
cdsBed=/mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/vcf/Mus_musculus.GRCm38.99.onlyCDS.bed

/usr/local/sw/bedops-2.4.41/bin/sort-bed ${cdsBed} \
| awk -vOFS="\t" '{ print $1, $2, $3, $4 }' - > /mnt/NEOGENE4/projects/melih_2020/mus_musculus/signal.bed

for TajimaFile in *_mus_musculus_tajimasD
do

outFile=/mnt/NEOGENE4/projects/melih_2020/data/mus_musculus/pop_gen/${TajimaFile}_cds_map.bed

#removing first line from the tajima file
sed -i -e "1d" ${TajimaFile}

/usr/local/sw/bedops-2.4.41/bin/bedmap --echo --echo-map --delim '\n' \
--multidelim '\n' ${TajimaFile} \
/mnt/NEOGENE4/projects/melih_2020/mus_musculus/signal.bed \
> ${outFile}

done

rm /mnt/NEOGENE4/projects/melih_2020/mus_musculus/signal.bed
