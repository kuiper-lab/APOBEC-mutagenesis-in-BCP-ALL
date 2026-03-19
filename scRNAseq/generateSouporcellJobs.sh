####
SCHEDULER="SLURM"

#### Parameters to set
BATCH_DIR="APOBEC_LX343"
PROJECT_DIR="/hpc/pmc_kuiper/HypermutatedALL_project"
SCRNASEQ_DIR="${PROJECT_DIR}/RESULTS/scRNAseq"
SOUPORCELL_DIR="${PROJECT_DIR}/RESULTS/souporcell/${BATCH_DIR}"
JOB_DIR="${PROJECT_DIR}/CODE/jobs/souporcell/${BATCH_DIR}"


# Analysis-specific variables
k=5       # Define the number of k's
OUTPUT_DIR="${SOUPORCELL_DIR}/k${k}"
BARCODES_FILE="${SCRNASEQ_DIR}/LX343_nonCellLineBarcodes.tsv.gz"
BAM_FILE="${SOUPORCELL_DIR}/possorted_genome_bam.bam"

mkdir -p ${JOB_DIR}
mkdir -p ${OUTPUT_DIR}


# Souporcell parameters
DOCKER_IMAGE_DIR="/hpc/local/Rocky8/pmc_kuiper/software/souporcell"
DOCKER_IMAGE="${DOCKER_IMAGE_DIR}/souporcell_latest.sif"

# Reference genome parameters
REF_DIR="/hpc/local/CentOS7/gen/data/10X-refdata/refdata-gex-GRCh38-2020-A/fasta"
REF_GENOME="${REF_DIR}/genome.fa"

# Parse paths of input directories and files to singularity syntax
REF_GENOME_PARSED="$(echo "/$(basename "$(dirname "${REF_GENOME}")")/$(basename "${REF_GENOME}")")"
BAM_FILE_PARSED="$(echo "/$(basename "$(dirname "${BAM_FILE}")")/$(basename "${BAM_FILE}")")"
BARCODES_FILE_PARSED="$(echo "/$(basename "$(dirname "${BARCODES_FILE}")")/$(basename "${BARCODES_FILE}")")"
SCRNASEQ_DIR_PARSED="$(echo "/$(basename "${SCRNASEQ_DIR}")")"


####
# Set run specific variables
DATE=$(date --iso-8601=seconds)
NUMCPUS=4
WALLTIME="23:59:00"
MEM="48G"
TMPSPACE="100G"


# Create job file
echo "Generating job for: ${BATCH_DIR}-k${k}"  

##########
####CREATE JOB HEADER
JOB_NAME="souporcell.${BATCH_DIR}-k${k}.sh"
JOB_FILE="${JOB_DIR}/${JOB_NAME}"
OUTPUT_LOG="${JOB_FILE}.out"
ERROR_LOG="${JOB_FILE}.err"

if [[ ${SCHEDULER} == "SLURM" ]]
then
    cat <<- EOF > ${JOB_FILE}
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${OUTPUT_LOG}
#SBATCH --error=${ERROR_LOG}
#SBATCH --partition=cpu
#SBATCH --time=${WALLTIME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task ${NUMCPUS}
#SBATCH --mem=${MEM}
#SBATCH --gres=tmpspace:${TMPSPACE}
#SBATCH --nodes=1
#SBATCH --open-mode=append

EOF

elif [[ ${SCHEDULER} == "SGE" ]]
then

    cat <<- EOF > ${JOBFILE}
#$ -S /bin/bash
#$ -N ${JOBNAME}
#$ -o ${OUTPUTLOG}
#$ -e ${ERRORLOG}
#$ -l h_rt=${WALLTIME}
#$ -l h_vmem=${MEM}
#$ -l tmpspace=${TMPSPACE}
#$ -cwd


EOF

else
    echo "Type of scheduler not known: ${SCHEDULER}"
    exit
fi


echo -e """

set -e # exit if any subcommand or pipeline returns a non-zero status
set -u # exit if any uninitialised variable is used


startTime=\$(date +%s)
echo \"startTime: \$startTime\"

cd ${OUTPUT_DIR}

singularity -d exec \\
--home ${OUTPUT_DIR} \\
--bind \$TMPDIR:\$TMPDIR \\
--bind ${REF_DIR}:/fasta \\
--bind ${SCRNASEQ_DIR}:${SCRNASEQ_DIR_PARSED} \\
--bind ${SOUPORCELL_DIR}:/${BATCH_DIR} \\
${DOCKER_IMAGE} \\
souporcell_pipeline.py \\
-i ${BAM_FILE_PARSED} \\
-b ${BARCODES_FILE_PARSED} \\
-f ${REF_GENOME_PARSED} \\
-t ${NUMCPUS} \\
-k ${k} \\
-o \$PWD


#Retrieve and check return code
returnCode=\$?
echo \"Return code \${returnCode}\"

if [ \"\${returnCode}\" -eq \"0\" ]
then
	
	echo -e \"Return code is zero, process was succesfull\n\n\"
	
else
  
	echo -e \"\nNon zero return code not making files final. Existing temp files are kept for debugging purposes\n\n\"
	#Return non zero return code
	exit 1
	
fi


#Write runtime of process to log file
endTime=\$(date +%s)
echo \"endTime: \$endTime\"


#Source: http://stackoverflow.com/questions/12199631/convert-seconds-to-hours-minutes-seconds-in-bash

num=\$endTime-\$startTime
min=0
hour=0
day=0
if((num>59));then
    ((sec=num%60))
    ((num=num/60))
    if((num>59));then
        ((min=num%60))
        ((num=num/60))
        if((num>23));then
            ((hour=num%24))
            ((day=num/24))
        else
            ((hour=num))
        fi
    else
        ((min=num))
    fi
else
    ((sec=num))
fi
echo \"Running time: \${day} days \${hour} hours \${min} mins \${sec} secs\"

""" >> ${JOB_FILE}

