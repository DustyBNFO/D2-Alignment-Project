#!/bin/bash

#PBS -S /bin/bash

#

# FastP Test Script

# This script is designed to run the STAR genome index builder

#

# Submit to queue

#PBS -q workq

#

# Request 8 processors on 1 nodes

#PBS -l nodes=1:ppn=8

#

# Request 24 hours of walltime

#PBS -l walltime=24:00:00

#

# Request 36 GB memory per process

#PBS -l pmem=36gb

#

# Contact zeliffdj@vcu.edu

#PBS -M zeliffdj@vcu.edu

#

# Send updates on abort, begin, and end

#PBS -m abe

#

# Name job

#PBS -N D2_Alignment_Index_STAR

#

# Return log

#PBS -o D2_Alignment_Index_Log.out

#

# Return errors

#PBS -e D2_Alignment_Index_error.err

#

# Set working directory to same as this script

cd $PBS_O_WORKDIR

#
# Source conda from user's home folder install
source /home/zeliffdj/miniconda3/etc/profile.d/conda.sh

# Activate fastp conda environment
conda activate /home/projects/MilesLab/teamshare/DZ_D2_Alignment/envs/star

# Run FastQC (and record time) and store to output folder

#Genome BUILDER

THREADS=8

​

STAR --runThreadN 11 --runMode genomeGenerate \
--genomeDir /home/projects/MilesLab/teamshare/DZ_D2_Alignment/genomes \
--genomeFastaFiles /home/projects/MilesLab/teamshare/DZ_D2_Alignment/genomes/GCA_9219983152_FASTA_Converted_DZ_3_23_23.fasta \
--sjdbGTFfile /home/projects/MilesLab/teamshare/DZ_D2_Alignment/genomes/DBA_2J_v3.2.gff3 \
--sjdbOverhang 150 --limitGenomeGenerateRAM 25000000000 



