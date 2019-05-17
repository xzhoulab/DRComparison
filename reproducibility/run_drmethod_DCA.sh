#!/bin/bash

#SBATCH --time=3-24:00:00
#SBATCH --job-name=dca
#SBATCH --mem=10000
#SBATCH --partition=nomosix
#SBATCH --array=1-330%50
#SBATCH --output=~/scRNAseqDRComparison/code/out/dca%a.out
#SBATCH --error=~/scRNAseqDRComparison/code/err/dca%a.err
#SBATCH --workdir=~/scRNAseqDRComparison/code

bash

let k=0

for idata in 'Baron'; do
for ip in 8 20 38 58; do
for ((irpt=1; irpt<=5; irpt++)); do
  GFILE=../data/sce_full_${idata}_DCA.tsv
  OUTPATH=~/scRNAseqDRComparison/results/${idata}
  let k=${k}+1
  if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
  cd ${OUTPATH}
  if [ ! -d "res.nPC${ip}.rpt${irpt}" ]; then
  ~/miniconda3/bin/dca ${GFILE} -s 512,${ip},512 res.nPC${ip}.rpt${irpt}
  fi
  fi	
done
done
done
