#!/bin/bash

#SBATCH --time=3-24:00:00
#SBATCH --job-name=scscope
#SBATCH --mem=10000
#SBATCH --partition=nomosix
#SBATCH --array=1-100%50
#SBATCH --output=~/scRNAseqDRComparison/code/out/scscope%a.out
#SBATCH --error=~/scRNAseqDRComparison/code/err/scscope%a.err
#SBATCH --workdir=~/scRNAseqDRComparison/code

bash

let k=0

for idata in 'Baron'; do
for ip in 8 20 38 58; do
for ((irpt=1; irpt<=5; irpt++)); do
  GFILE=../data/sce_${idata}_scScope.csv
  CODEPATH=../algorithm
  let k=${k}+1
  if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
  cd ${CODEPATH}
  python3.6 call_scScope.py --counts_file ${GFILE} --num_pc ${ip} --out ~/scRNAseqDRComparison/results/${idata}/res.${idata}.nPC${ip}.rpt${irpt}.scScope.txt
  fi	
done
done
done
