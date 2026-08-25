#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --job-name=iDEA
#SBATCH --output=iDEA.out
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --no-requeue

module load R/4.4.1
Rscript iDEA.analysis.cl.R GO_anno_iDEA.RDS DEG_iDEA.RDS idea_GO
