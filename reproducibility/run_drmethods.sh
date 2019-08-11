#!/bin/bash

#SBATCH --time=13-50:00:00
#SBATCH --job-name=dr
#SBATCH --mem=5000
#SBATCH --partition=nomosix
#SBATCH --array=1-260%100
#SBATCH --output=~/scRNAseqDRComparison/code/out/dr%a.out
#SBATCH --error=~/scRNAseqDRComparison/code/err/dr%a.err
#SBATCH --workdir=~/scRNAseqDRComparison/code

bash

let k=0

for id in 1; do
for ((im=1; im<=13; im++)); do
for ((ip=1; ip<=4; ip++)); do
for ((irpt=1; irpt<=5; irpt++)); do
  let k=${k}+1
  if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
  Rscript --verbose ./run_drmethods.R ${id} ${im} ${ip} ${irpt}
  fi	
done
done
done
done
