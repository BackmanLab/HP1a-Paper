#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=tiffanykuo2027@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH -t 24:00:00
#SBATCH --job-name="HP1a Alignment"

# Load modules
module load STAR/2.5.2

# Set working directory
cd /projects/b1042/BackmanLab/Tiffany/HP1a/seq_raw

# Commands to execute
genomefold=/home/tkx1376/Reference_Genomes/ref_hg38ensembl
outfold="$(dirname "$(pwd)")/${PWD##*/}_output"
[ ! -d $outfold ]&&mkdir $outfold
for file in *.fastq.gz;
do echo "$file";
STAR --runThreadN 16 --quantMode TranscriptomeSAM --genomeDir $genomefold --readFilesIn <(gunzip -c "${file}") --outFileNamePrefix "$outfold/${file%.*.*}";
done
