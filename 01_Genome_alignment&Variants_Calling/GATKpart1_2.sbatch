#!/bin/bash
### Gatk4 Germline SNPs+Indes Pipeline, part I: Data pre-processing (single sample)
###
### Inherit all current environment variables
#SBATCH --get-user-env
### Job name
#SBATCH --job-name 1.surban
### Keep Output and Error
#SBATCH -o job-%x-%j.out
### Queue name
#SBATCH -p MarkDupEBS
### Specify the number of nodes and thread (ppn) for your job.
#SBATCH -n 1
##SBATCH --nodelist=
#SBATCH --cpus-per-task=12
#SBATCH --mem=56G

echo Pipeline partI MarkDup Start time is `date +%Y/%m/%d--%H:%M`
work_start_time=$(date +%s)
hostname

# Sample ID
sample=surban
dataset=NyuWa
WORKDIR=/data
#INPUTDIR=$input_dir
threads=$(nproc)
#mem=64
batch=


# Input files
## get input reads from command line parameters
#bam=$bam

### need to edit the next 1 line

IMGDIR="/shared/data/images"
TMPDIR="/data/tmp"
[ -e $TMPDIR ] && : || mkdir -p $TMPDIR 

# genome and annotations
#ANNODIR=/wgs/anno/hg38/v0
#REFERENCE=${ANNODIR}/Homo_sapiens_assembly38.fasta
#bwa_genome=${ANNODIR}/Homo_sapiens_assembly38.fasta.64
#path2adaptor=/wgs/database/adaptors/adaptor.fa

# Log files
C_LOGS=/shared/data/nw8k/log/$batch$sample.complete.log
U_LOGS=/shared/data/nw8k/log/$batch$sample.uncomplete.log

####++++++++++++++++ Download Reads +++++++++++++++++
#
TOSUTIL=/shared/home/chenlin/500-samples-sbatch/scripts/tosutil
TOSCP="$TOSUTIL cp -p 2 -j 5 -conf /shared/home/zhpn/.tosutilconfig "

echo Download Start Time:`date +%Y/%m/%d--%H:%M`

mappingout=${WORKDIR}/03_BwaOut/${sample}
[ -e ${mappingout} ] && : || mkdir -p ${mappingout}
bwaout=${mappingout}/${sample}.bam
if [ ! -e $bwaout ];then
  ${TOSCP} tos://nw8k-results/bams/$dataset/$sample/$sample.bam $bwaout >  /tmp/bamdownload.log 2>&1
fi
echo Download End Time:`date +%Y/%m/%d--%H:%M`

#++++++++++++++++ MarkDuplicatesSpark +++++++++++++++++
gatkoutdir=${WORKDIR}/04_GatkOut/${sample}
[ -e ${gatkoutdir} ] && : || mkdir -p ${gatkoutdir}
gatkout=${gatkoutdir}/${sample}

echo "#######-- MarkDuplicatesSpark Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)
if [ ! -e ${gatkout}.Mdup.bam ];then
  singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp \
    ${IMGDIR}/gatk4_4.5.0.sif gatk --java-options "-Xmx40G -XX:ParallelGCThreads=12" \
    MarkDuplicatesSpark --spark-master local[12]  --conf 'spark.local.dir=/tmp/' -I ${bwaout} \
    -O ${gatkout}.Mdup.bam  --remove-sequencing-duplicates --read-validation-stringency LENIENT \
    --tmp-dir /tmp/ --conf spark.network.timeout=200h --conf spark.executor.heartbeatInterval=100h --spark-verbosity WARN
  if [ $? -ne 0 ]; then
    echo "MarkDuplicatesSpark failed"
    echo "${sample} MarkDuplicates UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    exit -1;
  fi
fi
## upload markdup bam
${TOSCP} ${gatkout}.Mdup.bam tos://nw8k-results/mdupbams/$dataset/$sample/ >> /tmp/upload.log 2>&1

if [ $? -ne 0 ]; then
  echo "MarkDuplicatesSpark Upload failed"
  echo "${sample} MarkDuplicates_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  exit -1;
fi

${TOSCP} ${gatkout}.Mdup.bam.bai tos://nw8k-results/mdupbams/$dataset/$sample/ >> /tmp/upload.log 2>&1

echo "${sample} MarkDuplicates Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS

rm -f ${bwaout}
rm -rf ${gatkoutdir}
end_time=$(date +%s)
echo "#######-- MarkDuplicatesSpark End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- MarkDuplicatesSpark Elapsed Time:$elapsed_time sec --#######"
