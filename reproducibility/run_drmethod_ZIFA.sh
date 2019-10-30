#!/bin/bash

#SBATCH --time=3-24:00:00
#SBATCH --job-name=zifa
#SBATCH --mem=10000
#SBATCH --partition=nomosix
#SBATCH --array=1-20
#SBATCH --output=~/scRNAseqDRComparison/code/out/zifa%a.out
#SBATCH --error=~/scRNAseqDRComparison/code/err/zifa%a.err
#SBATCH --workdir=~/scRNAseqDRComparison/code

bash

let k=0

for idata in 'Baron'; do
for ip in 8 20 38 58; do
for ((irpt=1; irpt<=5; irpt++)); do
  GFILE=../data/sce_${idata}_ZIFA.txt
  CODEPATH=../algorithm/ZIFA
  let k=${k}+1
  if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
  cd ${CODEPATH}
  if [ ! -d "~/scRNAseqDRComparison/results/${idata}/res.${idata}.nPC${ip}.rpt${irpt}.ZIFA.txt" ]; then
  ~/anaconda2/bin/python call_ZIFA.py --normcounts_file ${GFILE} --num_pc ${ip} --out ~/scRNAseqDRComparison/results/${idata}/res.${idata}.nPC${ip}.rpt${irpt}.ZIFA.txt
  fi
  fi	
done
done
done
