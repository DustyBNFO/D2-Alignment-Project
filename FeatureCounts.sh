#!/bin/bash
#PBS -S /bin/bash
#
# Qualimap Script
# This script runs Qualimap on DO Mouse NAc RNA samples for both QC and trimming.
#
# Submit to working queue
#PBS -q serial
#
# Request 8 processors on 1 node
#PBS -l nodes=1:ppn=2
#
# Request 12 hours of walltime
#PBS -l walltime=24:00:00
#
# Request 2 GB memory per process
#PBS -l pmem=48000MB
#
# Give access to environment variables
#PBS -V
#PBS -M zeliffdj@vcu.edu
#PBS -m ae
#PBS -N featureCounts
#PBS -j oe
#PBS -o /home/projects/MilesLab/teamshare/DZ_D2_Alignment/logs/


cd $PBS_O_WORKDIR

# Source conda from user's home folder install
source /home/zeliffdj/miniconda3/etc/profile.d/conda.sh

# Activate subread conda environment
conda activate /home/projects/MilesLab/teamshare/DZ_D2_Alignment/envs/subread

# Path to genome annotation files
# DIRINDEX=/home/mdozmorov/sequencing/data/ExtData/UCSC/hg38gdc
# INDEXGENE=gencode.v22.annotation.gtf.gz

DIRINDEX=/home/projects/MilesLab/teamshare/DZ_D2_Alignment/genomes
INDEXGENE=/DBA_2J_v3.2_3_14_23.gtf

# DIRINDEX=/home/sequencing/data/ExtData/UCSC/rn6
# INDEXGENE=Rattus_norvegicus.Rnor_6.0.85_chr.gtf

THREADS=24

# Output file with gene counts
COUNTS=counts.txt

DIRIN=/home/projects/MilesLab/teamshare/DZ_D2_Alignment/D_samsort_STAR
DIROUT=/home/projects/MilesLab/teamshare/DZ_D2_Alignment/F_feature_counts

mkdir -p $DIROUT

# Summarize counts per gene_id.
# For single-end, remove "-p" option
# Process all files at once
# ls -p 02_bwa_trimmed/*.bam | tr '\n' ','
# featureCounts -a ${DIRINDEX}/${INDEXGENE} -o ${COUNTS} -p -T ${THREADS} 02_bwa_trimmed/E0771-CT1-B_bwa.bam,02_bwa_trimmed/E0771-CT2-A_bwa.bam,02_bwa_trimmed/E0771-CT3-A_bwa.bam,02_bwa_trimmed/E0771-DT1-B_bwa.bam,02_bwa_trimmed/E0771-DT2-B_bwa.bam,02_bwa_trimmed/E0771-DT3-B_bwa.bam,02_bwa_trimmed/PyMT-CT1-A_bwa.bam,02_bwa_trimmed/PyMT-CT2-A_bwa.bam,02_bwa_trimmed/PyMT-CT3-A_bwa.bam,02_bwa_trimmed/PyMT-DT1-A_bwa.bam,02_bwa_trimmed/PyMT-DT2-A_bwa.bam,02_bwa_trimmed/PyMT-DT3-A_bwa.bam


# Process each file individually
for file in `find $DIRIN -name "*.bam" -type f | sort`; do
	featureCounts -T ${THREADS} -p -s 2 -t exon -g gene_id \
		-a $DIRINDEX"/"$INDEXGENE \
		-o $DIROUT"/"`basename $file .Aligned.out.bam`".txt" $file;
done