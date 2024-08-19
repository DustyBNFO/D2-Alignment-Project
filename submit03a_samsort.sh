#!/bin/bash
#PBS -S /bin/bash
#
# Samsort Script
# This script runs samtools sort on aligned DeepSeq RNA samples.
#
# Submit to working queue
#PBS -q workq
#
# Request 8 processors on 1 node
#PBS -l nodes=1:ppn=8
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
#PBS -N D2_samsort
#
# Run a job array: all 1-400 samples
#PBS -J 11-50
#
# Folder to put output and error logs in
#PBS -o /home/projects/MilesLab/teamshare/D2_Alignment/logs/
#
# Put output and error logs in one single file per job
#PBS -j oe

# Set working directory to wherever script was submitted from
cd $PBS_O_WORKDIR

# Set number of threads (should be same as number of processors above)
THREADS=8

# Input folder for where aligned RNA-seq data are
DIRIN=/home/projects/MilesLab/teamshare/D2_Alignment/outputs/02_STAR-align

# Name of output folder to put sorted samples
DIROUT=/home/projects/MilesLab/teamshare/D2_Alignment/output/03_STAR-sorted

# Create the output folder
mkdir -p ${DIROUT}

# Create a variable to store an individual alignment's name in
R1=${DIRIN}*_S${PBS_ARRAY_INDEX}.Aligned.sortedByCoord.out.bam

# Run samtools sort
samtools sort \
  -o ${DIROUT}/`basename ${R1}` \
  -O bam \
  --threads ${THREADS} \
  $R1

# Run samtools index
samtools index ${DIROUT}/`basename ${R1}`

# for file in `find ${DIRIN} -type f -name "*.bam" | sort`; do
#	FILEOUT=${DIROUT}/`basename $file`_sorted.bam
#	samtools sort -o ${FILEOUT} -O bam --threads ${THREADS} ${file}
#	samtools index ${FILEOUT}
#	rm $file
#done
