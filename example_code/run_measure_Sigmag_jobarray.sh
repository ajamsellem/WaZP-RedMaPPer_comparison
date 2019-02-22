#!/bin/bash -l

#SBATCH --job-name=splashback_%j
#SBATCH --time=05:00:00
##SBATCH --nodes=1
#SBATCH --partition=sandyb
#SBATCH --account=pi-sdodelso
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=32000
#SBATCH --output=log/splashback-%A_%a.out
#SBATCH --error=log/splashback-%A_%a.err
#SBATCH --array=0-99%100
#SBATCH --exclusive

#module load python

main_dir=/project/kicp/chihway/brutus/splashback/
jkid=${SLURM_ARRAY_TASK_ID}

mkdir ${main_dir}analysis/results/l20_100_fiducial_bin20/

if ! [ -f "${main_dir}analysis/results/l20_100_fiducial_bin20/Sigmag_${jkid}_6.npz" ]
then

for ((i=0;i<=6;i++))
do

python measure_Sigmag.py 0.2 0.55 7 -100 -20.2356590685 20 100 ${main_dir}data/des_y1_redmapper_jk_fiducial.fits ${main_dir}data/des_y1_redmapper_randoms_jk_fiducial.fits ${main_dir}data/des_y1_gold_jk_fiducial.fits ${main_dir}data/des_y1_gold_randoms_jk_fiducial.fits $jkid 1392.13304984  ${main_dir}analysis/results/l20_100_fiducial_bin20/Sigmag_${jkid}_${i}.npz 20 ${i}

done

fi



