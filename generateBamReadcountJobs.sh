####
SCHEDULER="SLURM"

#### Parameters to set
BATCHDIR="DDOST"
PROJECTDIR="/hpc/pmc_kuiper/HypermutatedALL_project/"
CODEDIR="${PROJECTDIR}/CODE/SHELL/bam-readcount/"
DATADIR="${PROJECTDIR}/DATA/CRAM/"
JOBDIR="${PROJECTDIR}/CODE/jobs/bam-readcount/${BATCHDIR}/"
OUTPUTDIR="${PROJECTDIR}/ANALYSIS/bam-readcount/${BATCHDIR}/"
TARGET_REGIONS="${CODEDIR}/DDOST_hotspot_location.bed"
EXT="_RNA-Seq.cram"


mkdir -p ${JOBDIR}
mkdir -p ${OUTPUTDIR}


# Bam-readcount parameters
DOCKER_IMAGE_DIR="/hpc/local/Rocky8/pmc_kuiper/software/bam-readcount/"
DOCKER_IMAGE="${DOCKER_IMAGE_DIR}/bam-readcount_latest.sif"

# Reference genome parameters
REF_DIR="/hpc/pmc_gen/references/RNA-Seq/references/"
REF_GENOME="${REF_DIR}/ref_genome_GRCh38_gencode_v31_CTAT_lib_Oct012019.fa"


####
# Set run specific variables
DATE=$(date --iso-8601=seconds)
NUMCPUS=1          
WALLTIME="00:30:00"
MEM="20G"
TMPSPACE="50G"


# Read sample file and create job file
for CRAM in $(ls ${DATADIR}/*${EXT})
do
    # Define the sample ID
	SAMPLE_ID=$(basename "${CRAM}" "${EXT}")

    echo "Generating job for sample: ${SAMPLE_ID}"
  

##########
####CREATE JOB HEADER
JOBNAME="bam-readcount.${SAMPLE_ID}.sh"
JOBFILE="${JOBDIR}/${JOBNAME}"
OUTPUTLOG="${JOBDIR}/bam-readcount.${SAMPLE_ID}.sh.out"
ERRORLOG="${JOBDIR}/bam-readcount.${SAMPLE_ID}.sh.err"

if [[ ${SCHEDULER} == "SLURM" ]]
then
    cat <<- EOF > ${JOBFILE}
#!/bin/bash
#SBATCH --job-name=${JOBNAME}
#SBATCH --output=${OUTPUTLOG}
#SBATCH --error=${ERRORLOG}
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


## Run bam-readcount
apptainer -d exec --writable-tmpfs --containall \\
--bind \$TMPDIR:\$TMPDIR \\
--bind ${REF_DIR}:${REF_DIR} \\
--bind ${DATADIR}:${DATADIR} \\
--bind ${CODEDIR}:${CODEDIR} \\
--bind ${DOCKER_IMAGE_DIR}:${DOCKER_IMAGE_DIR} \\
--bind ${OUTPUTDIR}:${OUTPUTDIR} \\
${DOCKER_IMAGE} \\
bam-readcount -w 10 \\
-f ${REF_GENOME} \\
-l ${TARGET_REGIONS} \\
${CRAM} >> ${OUTPUTDIR}/${SAMPLE_ID}-readcount.variants.txt

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

""" >> ${JOBFILE}
    
done

