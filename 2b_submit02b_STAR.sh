#!/bin/bash
#PBS -S /bin/bash
#
# STAR Alignment Script
# This script runs STAR on DeepSeq RNA samples for alignment to a D2 reference genome. This genome should already be indexed before running this script.
#
# Submit to working queue
#PBS -q workq
#
# Request 12 processors on 1 node
#PBS -l nodes=1:ppn=12
#
# Request 12 hours of walltime
#PBS -l walltime=12:00:00
#
# Request 128 GB memory per process
#PBS -l pmem=128gb
#
# Give access to environment variables
#PBS -V
#
# Contact zeliffdj@vcu.edu
#PBS -M zeliffdj@vcu.edu
#
# Send updates on abort and end
#PBS -m ae
#
# Name job batch
#PBS -N D2_STAR_align
#
# Run a job array for all 10 samples
#
#
# Folder to put output and error logs in
#PBS -o /home/projects/MilesLab/teamshare/DZ_D2_Alignment/logs/
#
# Put output and error logs in one single file per job
#PBS -j oe

# Set working directory to wherever script was submitted from
cd $PBS_O_WORKDIR

# Set number of threads to be used by STAR (should be same as number of processors above)
THREADS=8

# Input folder for where trimmed RNA-seq data are
DIRIN=/home/projects/MilesLab/teamshare/B6D2_DeepSeq/B_fastp/

# Name of output folder to put aligned samples
DIROUT=/home/projects/MilesLab/teamshare/DZ_D2_Alignment/outputs/02_STAR-align

# Set location of mouse genome files
DIRINDEX=/home/projects/MilesLab/teamshare/DZ_D2_Alignment/genomes

# Set name of mouse fasta file
DNAINDEX=$DIRINDEX"/GCA_9219983152_FASTA_Converted_DZ_3_23_23.fasta"

# Set name of mouse gtf annotations file
GENINDEX=$DIRINDEX"/DBA_2J_v3.2.gff3"

# Create the output folder
mkdir -p $DIROUT
chmod 777 -R ${DIROUT}

# Create a variable to store the reads for the sample matching the job's index number
#R1=${DIRIN}*_S${PBS_ARRAY_INDEX}_*_R1_001.fastq.gz
#R2=${DIRIN}*_S${PBS_ARRAY_INDEX}_*_R2_001.fastq.gz
R1="/home/projects/MilesLab/teamshare/DZ_D2_Alignment/B6D2_DeepSeq_B_fastp_copy/B_fastp/D11N_S1_R1_001.fastq.gz"
R2="/home/projects/MilesLab/teamshare/DZ_D2_Alignment/B6D2_DeepSeq_B_fastp_copy/B_fastp/D11N_S1_R2_001.fastq.gz"


# Create a variable for output file name
prefix=${DIROUT}/`basename $R2 _L001_R2_001`

#Add Miniconda lines here
# Source conda from user's home folder install
source /home/zeliffdj/miniconda3/etc/profile.d/conda.sh

# Activate fastp conda environment
conda activate /home/projects/MilesLab/teamshare/DZ_D2_Alignment/envs/star

# Run STAR on paired-end data
STAR \
  --runThreadN $THREADS \
  --genomeDir $DIRINDEX \
  --readFilesIn $R1 $R2 \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMorder Paired \
  --outReadsUnmapped Fastx \
  --sjdbGTFfile $GENINDEX \
  --quantMode GeneCounts \
  --outFilterMultimapNmax 10 \
  --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
  --outFileNamePrefix ${prefix}.

echo "Completed STAR alignment on "`date`

# Duplicate removal: https://github.com/alexdobin/STAR/issues/175
