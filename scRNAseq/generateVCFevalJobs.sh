# Global parameters
LIBRARY="LX448"
N_SAMPLES="4"
PROJECT_DIR="/hpc/pmc_kuiper/HypermutatedALL_project/"
JOB_DIR="${PROJECT_DIR}/CODE/jobs/vcfeval/${LIBRARY}/"
OUTPUT_DIR="${PROJECT_DIR}/RESULTS/vcfeval/${LIBRARY}/k${N_SAMPLES}/"

# VCF parameters
VCF_INPUT_DIR="${PROJECT_DIR}/DATA/VCF/"
WGS_VCF_EXT="_WGS.HyperExomeFiltered.DPgte20.ALTgte3.filtered.uniqueEntries.vcf.gz"
SC_VCF="${VCF_INPUT_DIR}/k${N_SAMPLES}_cluster_genotypes.vcf.gz"

# VCF-eval parameters
# VCFeval is a function of the RTG tool
#
# RTG-tools is downloaded with miniconda as installing from source resulted in error messages
# conda install -c bioconda rtg-tools
RTG_TOOLS_PATH="/hpc/local/Rocky8/pmc_kuiper/software/miniconda-24.9.2/share/rtg-tools-3.13-0/rtg"
REF_GENOME_SDF="${PROJECT_DIR}/CODE/SHELL/scRNAseq/GRCh38.sdf"


mkdir -p "${OUTPUT_DIR}"
mkdir -p "${JOB_DIR}"


# Parse each VCF in the input directory
for WGS_VCF in $( ls ${VCF_INPUT_DIR}/*${WGS_VCF_EXT} )
do
SAMPLE=$( basename ${WGS_VCF} ${WGS_VCF_EXT} )
SAMPLE_VCF_HEADER=$(grep "^#CHROM" "$WGS_VCF" | cut -f11)

echo "Generating job for sample ${SAMPLE} ..";


# Specify job requirements
JOBNAME="${SAMPLE}-vcfeval"
JOBFILE="${JOB_DIR}/${JOBNAME}.sh"
OUTPUTLOG="${JOBFILE}.out"
ERRORLOG="${JOBFILE}.err"
WALLTIME="00:15:00"
NUMTASKS="1"
NUMCPUS="1"
MEM="5G"
TMPSPACE="5G"

# Create job file
cat <<- EOF > ${JOBFILE}
#!/bin/bash
#SBATCH --job-name=${JOBNAME}
#SBATCH --output=${OUTPUTLOG}
#SBATCH --error=${ERRORLOG}
#SBATCH --partition=cpu
#SBATCH --time=${WALLTIME}
#SBATCH --ntasks=${NUMTASKS}
#SBATCH --cpus-per-task ${NUMCPUS}
#SBATCH --mem=${MEM}
#SBATCH --gres=tmpspace:${TMPSPACE}
#SBATCH --nodes=1
#SBATCH --open-mode=append

set -e # exit if any subcommand or pipeline returns a non-zero status
set -u # exit if any uninitialised variable is used

startTime=\$(date +%s)
echo "startTime: \$startTime"

# Extract tumor sample ID from VCF file
SAMPLE_ID_VCF=\$(zgrep "^#CHROM" "$WGS_VCF" | cut -f10)


echo "Running VCF-eval.."

# Execute command for each cluster in the sc-library
for ((K=0; K<${N_SAMPLES}; K++))
do
    ${RTG_TOOLS_PATH} vcfeval \\
    -t ${REF_GENOME_SDF} \\
    -b ${SC_VCF} \\
    -c ${WGS_VCF} \\
    --sample \${K},\${SAMPLE_ID_VCF} \\
    -o ${OUTPUT_DIR}/K${N_SAMPLES}.\${K}_${SAMPLE}
done

echo "Finished running VCF-eval.."

echo "Combining VCF-eval results.."
cd ${OUTPUT_DIR}
for f in K${N_SAMPLES}*${SAMPLE}; do echo \${f}; cat \${f}/summary.txt; done > \\
${OUTPUT_DIR}/${LIBRARY}-vcfeval-summaries-${SAMPLE}-K${N_SAMPLES}.txt

# Remove tmp results
rm -r K${N_SAMPLES}*${SAMPLE}

echo "Finished all steps.."

# Retrieve and check return code
returnCode=\$?
echo "Return code \${returnCode}"

if [ \${returnCode} -eq 0 ]
then
        echo -e "Return code is zero, process was succesfull\n\n"
else
        echo -e "\nNon zero return code not making files final. Existing temp files are kept for debugging purposes\n\n"
        exit 1
fi

# Write runtime of process to log file
# Source: http://stackoverflow.com/questions/12199631/convert-seconds-to-hours-minutes-seconds-in-bash
endTime=\$(date +%s)
echo "endTime: \$endTime"

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
echo "Running time: \${day} days \${hour} hours \${min} mins \${sec} secs"

EOF

done
