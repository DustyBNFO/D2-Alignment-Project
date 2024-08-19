#!/bin/bash
#PBS -S /bin/bash
#
# Qualimap Script
# This script runs Qualimap on DO Mouse NAc RNA samples for both QC and trimming.
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
# Request 36 GB memory per process
#PBS -l pmem=36gb
#
# Give access to environment variables
#PBS -V
#
# Contact tatomz@vcu.edu
#PBS -M tatomz@vcu.edu
#
# Send updates on abort and end
#PBS -m ae
#
# Name job
#PBS -N D2_featureCounts
#
# Folder to put output and error logs in
#PBS -o /home/projects/MilesLab/teamshare/D2_Alignment/logs/
#
# Put output and error logs in one single file per job
#PBS -j oe

# Set working directory to wherever script was submitted from
cd $PBS_O_WORKDIR

# Source conda from user's home folder install
source /home/tatomz/miniconda3/etc/profile.d/conda.sh

# Activate subread conda environment
conda activate /home/projects/MilesLab/teamshare/ZT_DO_NAc/envs/subread

# Set number of threads to be used by featureCounts (should be same as number of processors above)
THREADS=8

# Set location of mouse genome files
DIRINDEX=/home/projects/MilesLab/teamshare/D2_Alignment/genomes/

# Set name of mouse gtf annotations file
INDEXGENE=Mus_musculus.GRCm39.108.gtf

# Output file with gene counts
COUNTS=counts.txt

# Input folder for where sorted RNA-seq alignment data are
DIRIN=/home/projects/MilesLab/teamshare/D2_Alignment/output/03_STAR-sorted

# Output folder to store counts
DIROUT=/home/projects/MilesLab/teamshare/D2_Alignment/output/04_featureCounts/

# Create the output folder
mkdir -p $DIROUT

# Create a variable to store an individual alignment's name in
R1=${DIRIN}/*_S${PBS_ARRAY_INDEX}_*.Aligned.sortedByCoord.out.bam

# Run featureCounts
featureCounts \
  -T ${THREADS} \
  -p \
  -t exon \
  -g gene_id \
  -a ${DIRINDEX}/${INDEXGENE} \
  -o ${DIROUT}/"all_samples.txt" \
  -s 2 \
  $(ls ${DIRIN}/*.Aligned.sortedByCoord.out.bam)

# Summarize counts per gene_id.
# For single-end, remove "-p" option
# Process all files at once
# ls -p 02_bwa_trimmed/*.bam | tr '\n' ','
# featureCounts -a ${DIRINDEX}/${INDEXGENE} -o ${COUNTS} -p -T ${THREADS} 02_bwa_trimmed/E0771-CT1-B_bwa.bam,02_bwa_trimmed/E0771-CT2-A_bwa.bam,02_bwa_trimmed/E0771-CT3-A_bwa.bam,02_bwa_trimmed/E0771-DT1-B_bwa.bam,02_bwa_trimmed/E0771-DT2-B_bwa.bam,02_bwa_trimmed/E0771-DT3-B_bwa.bam,02_bwa_trimmed/PyMT-CT1-A_bwa.bam,02_bwa_trimmed/PyMT-CT2-A_bwa.bam,02_bwa_trimmed/PyMT-CT3-A_bwa.bam,02_bwa_trimmed/PyMT-DT1-A_bwa.bam,02_bwa_trimmed/PyMT-DT2-A_bwa.bam,02_bwa_trimmed/PyMT-DT3-A_bwa.bam


# Process each file individually
#for file in `find $DIRIN -name "*.bam" -type f | sort`; do
#	featureCounts -T ${THREADS} -p -t exon -g gene_id \
#		-a $DIRINDEX"/"$INDEXGENE \
#		-o $DIROUT"/"`basename $file .Aligned.out.bam`".txt" $file;
#done
