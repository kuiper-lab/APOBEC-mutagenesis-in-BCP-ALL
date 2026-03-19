####
SCHEDULER="SLURM"

#### Parameters to set
BATCHDIR="APOBEC"
EXT=".modelFinal.seg"
PROJECTDIR="/hpc/pmc_kuiper/HypermutatedALL_project/"
INPUT_DIR="$PROJECTDIR/ANALYSIS/CNVsingleSampleWorkflow/$BATCHDIR/"
OUTPUT_DIR="$PROJECTDIR/ANALYSIS/annotateCNVsegments/$BATCHDIR/"
JOB_DIR="$PROJECTDIR/CODE/jobs/annotateCNVsegments/$BATCHDIR/"
FASTA_FILE="/hpc/pmc_kuiper/References/hg38bundle/v0/Homo_sapiens_assembly38.fasta"

### BED file to filter on
PANELBEDFILE="$PROJECTDIR/REFERENCE/Homo_sapiens.GRCh38.98_protein_coding_genes.bed"
AUXFILE="$PROJECTDIR/CODE/functions_aux.R"
CODINGGENESBED="$PROJECTDIR/REFERENCE/Homo_sapiens.GRCh38.98_protein_coding_genes.bed"

### Annotation specific vars
REFERENCE_DIR="/hpc/pmc_kuiper/References/"
NUM_FORKS=1 

#### General settings
mkdir -p "$OUTPUT_DIR"
mkdir -p "$JOB_DIR"

for SEGFILE in $(ls $INPUT_DIR/*${EXT})
do	

BASE=$( basename $SEGFILE .seg )
SAMPLE=$( echo $BASE | awk '{print $1}' FS="_")

#Output variable
OUT_FILE="$OUTPUT_DIR/$BASE.annotated.txt"

echo -e "Creating job for sample $SAMPLE"

##########
####CREATE JOB HEADER

JOBNAME="$SAMPLE-annotateCNVsegments"
OUTPUTLOG="$JOB_DIR/$SAMPLE-annotateCNVsegments.sh.out"
ERRORLOG="$JOB_DIR/$SAMPLE-annotateCNVsegments.sh.err"
WALLTIME="00:59:00"
NUMTASKS="1"
NUMCPUS=$NUM_FORKS
MEM="20G"
TMPSPACE="10G"
JOBFILE="$JOB_DIR/$SAMPLE-annotateCNVsegments.sh"

if [[ $SCHEDULER == "SLURM" ]]
then
    cat <<- EOF > $JOBFILE
#!/bin/bash
#SBATCH --job-name=$JOBNAME
#SBATCH --output=$OUTPUTLOG
#SBATCH --error=$ERRORLOG
#SBATCH --partition=cpu
#SBATCH --time=$WALLTIME
#SBATCH --ntasks=$NUMTASKS
#SBATCH --cpus-per-task $NUMCPUS
#SBATCH --mem=$MEM
#SBATCH --gres=tmpspace:$TMPSPACE
#SBATCH --nodes=1
#SBATCH --open-mode=append

EOF

elif [[ $SCHEDULER == "SGE" ]]
then

    cat <<- EOF > $JOBFILE
#$ -S /bin/bash
#$ -cwd
#$ -o $OUTPUTLOG
#$ -e $ERRORLOG
#$ -N $JOBNAME
#$ -l h_rt=$WALLTIME
#$ -l h_vmem=$MEM
#$ -l tmpspace=$TMPSPACE
#$ -pe threaded $NUMCPUS

EOF

else
    echo "Type of scheduler not known: $SCHEDULER"
    exit
fi

echo -e """

set -e # exit if any subcommand or pipeline returns a non-zero status
set -u # exit if any uninitialised variable is used


startTime=\$(date +%s)
echo \"startTime: \$startTime\"


# load the required modules
module load R/3.6.1

echo \"Starting annotation ..\"

# Execute R script
R --slave --no-save --no-restore --no-environ \\
--args $SEGFILE \\
$PANELBEDFILE \\
$AUXFILE \\
$CODINGGENESBED \\
$OUT_FILE \\
< "$PROJECTDIR/CODE/annotateSCNVoutput.R"

echo \"Finished annotation\"

echo \"Starting md5summing..\"
cd $OUTPUT_DIR
md5sum ${BASE}.annotated.txt > ${BASE}.annotated.txt.md5
cd -
echo \"Finished md5summing\"

echo \"Finished all steps\"


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

""" >> $JOBFILE

done
