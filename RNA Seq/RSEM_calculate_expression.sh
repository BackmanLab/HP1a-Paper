#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=tiffanykuo2027@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH -t 48:00:00
#SBATCH --job-name="HP1a RSEM"

# load modules you need to use
module load rsem/1.3.0

# Set your working directory
cd /projects/b1042/BackmanLab/Tiffany/HP1a/seq_raw_output

# A command you actually want to execute:
genomefold=/home/tkx1376/Reference_Genomes/ref_hg38ensembl

for file in *.toTranscriptome.out.bam;
do rsem-calculate-expression -p 10 --bam $file "$genomefold/hg38_ensembl" ${file%.*.*.*};
done
