# Request 8 processors on 1 node
#PBS -l nodes=1:ppn=8
#
# Request 12 hours of walltime
#PBS -l walltime=24:00:00
#
# Request 128 GB memory per process
#PBS -l pmem=128gb
#
# Give access to environment variables
#PBS -V
#
# Contact gnatowskie@vcu.edu
#PBS -M zeliffdj@vcu.edu
#
# Send updates on abort and end
#PBS -m ae
#
# Name job batch
#PBS -N B6D2_samsort
#
#
# Put output and error logs in one single file per job
#PBS -j oe

# Set working directory to wherever script was submitted from
cd $PBS_O_WORKDIR

# Set number of threads (should be same as number of processors above)
THREADS=8

# Input folder for where aligned RNA-seq data are
DIRIN=/home/projects/MilesLab/teamshare/DZ_D2_Alignment/All_D2_Out_Copy/All_D2_Out/02_STAR-align

# Name of output folder to put sorted samples
DIROUT=/home/projects/MilesLab/teamshare/DZ_D2_Alignment/D_samsort_STAR/

# Create the output folder
mkdir -p ${DIROUT}

 for file in `find ${DIRIN} -type f -name "*.bam" | sort`; do
	FILEOUT=${DIROUT}/`basename $file`
	samtools sort -o ${FILEOUT} -O bam --threads ${THREADS} ${file}
	samtools index ${FILEOUT}
	rm $file
done
