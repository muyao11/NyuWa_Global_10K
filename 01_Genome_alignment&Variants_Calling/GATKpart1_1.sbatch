#!/bin/bash
### Gatk4 Germline SNPs+Indes Pipeline, part I: Data pre-processing (single sample)
###
### Inherit all current environment variables
#SBATCH --get-user-env
### Job name
#SBATCH --job-name 0.surban
### Keep Output and Error
#SBATCH -o job-%x-%j.out
### Queue name
#SBATCH -p C16G64
### Specify the number of nodes and thread (ppn) for your job.
#SBATCH -n 1
##SBATCH --nodelist=MarkDupEBS-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=56G

#sudo chmod 777 /data/

###sbatch --export=sample=SampleName  *.sbatch
echo Pipeline partI Start time is `date +%Y/%m/%d--%H:%M`
work_start_time=$(date +%s)
hostname

sample=surban
dataset=NyuWa
WORKDIR=/data
threads=$(nproc)
batch=
#mem=64
fastq=
#'tos://nw8k/056/RawData/h1-H8_L2_392392.R1.clean.fastq.gz tos://nw8k/056/RawData/h1-H8_L2_392392.R2.clean.fastq.gz tos://nw8k/056/RawData/h1-H8_L3_392392.R1.clean.fastq.gz tos://nw8k/056/RawData/h1-H8_L3_392392.R2.clean.fastq.gz'

rm -r ${WORKDIR}/01_CleanData/
rm -r ${WORKDIR}/03_BwaOut/
rm -r ${WORKDIR}/04_GatkOut/

IMGDIR="/shared/data/images"
TMPDIR="/data/tmp"

[ -e ${TMPDIR} ] && : || mkdir -p ${TMPDIR}

# genome and annotations
ANNODIR=/data/wgs/anno/hs1/v0
REFERENCE=${ANNODIR}/hs1/hs1.fasta
bwa_genome=${ANNODIR}/hs1/hs1.fasta
path2adaptor=${ANNODIR}/hs1/TruSeq3-PE-2.fa

# Log files
C_LOGS=/shared/data/nw8k/log/$batch$sample.complete.log
U_LOGS=/shared/data/nw8k/log/$batch$sample.uncomplete.log

####++++++++++++++++ Download Reference +++++++++++++++++

TOSUTIL=/shared/home/chenlin/500-samples-sbatch/scripts/tosutil
TOSCP="$TOSUTIL cp -p 2 -j 5 -conf /shared/home/zhpn/.tosutilconfig "

[ -e ${ANNODIR} ] && : || mkdir -p ${ANNODIR}
echo Download Start Time:`date +%Y/%m/%d--%H:%M`
if [ ! -e /data/wgs/download.hs1.finished ];then
    ${TOSCP} -r tos://bioinfo-databases/nw8k/hs1/  ${ANNODIR}  > /tmp/download.wgs.log 2>&1 && touch /data/wgs/download.hs1.finished
fi

####++++++++++++++++ Trim Adaptor +++++++++++++++++
cleandir=${WORKDIR}/01_CleanData/${sample}
[ -e ${cleandir} ] && : || mkdir -p ${cleandir}

out_paired_R1=${cleandir}/${sample}_R1_clean_paired.fastq.gz
out_unpaired_R1=${cleandir}/${sample}_R1_clean_unpaired.fastq.gz
out_paired_R2=${cleandir}/${sample}_R2_clean_paired.fastq.gz
out_unpaired_R2=${cleandir}/${sample}_R2_clean_unpaired.fastq.gz

if [ ! -e "${out_paired_R1}" ]; then

if [ ! -e ${WORKDIR}/01_CleanData/download.${sample}.finished ];then
  for i in ${fastq}; do
    echo $i
    bn=`basename $i`
    if [ ! -e ${cleandir}/$bn ]; then
      ${TOSCP}  ${i} ${cleandir}/ >> /tmp/fastq.download.log 2>&1 && touch ${WORKDIR}/01_CleanData/download.${sample}.finished
    else
      ${TOSCP}  ${i} ${cleandir}/$bn.1 >> /tmp/fastq.download.log 2>&1 && touch ${WORKDIR}/01_CleanData/download.${sample}.finished
    fi
  done    
fi

echo Trimmomatic Start Time:`date +%Y/%m/%d--%H:%M`
start_time=$(date +%s)

input_R1=`ls ${cleandir}/*[._]R1*.f*q.gz*`
if [ $? -ne 0 ]; then
  input_R1=`ls ${cleandir}/*[._]1*.f*q.gz*`
  input_R2=`ls ${cleandir}/*[._]2*.f*q.gz*`
else
  input_R2=`ls ${cleandir}/*[._]R2*.f*q.gz*`
fi

echo "input_R1: $input_R1"
#input_R2=`ls ${cleandir}/*[._]R2*.f*q.gz* ${cleandir}/*[._]2*.f*q.gz*`
echo "input_R2: $input_R2"

c=`ls $input_R1 | wc -l`
if [ $c -ne 1 ]; then
  cat $input_R1 >${cleandir}/${sample}_combine_1.fastq.gz
  cat $input_R2 >${cleandir}/${sample}_combine_2.fastq.gz

  ls -lh ${cleandir}
  rm $input_R1 $input_R2

  input_R1=${cleandir}/${sample}_combine_1.fastq.gz
  input_R2=${cleandir}/${sample}_combine_2.fastq.gz
fi

#if [ ! -e "${out_paired_R1}" ]; then
  singularity exec -B $WORKDIR:$WORKDIR  ${IMGDIR}/biotools_ibp_snp_v2.sif \
    bash -c "_JAVA_OPTIONS=' -Xmx24g' trimmomatic PE -threads $threads \
    -phred33 ${input_R1} ${input_R2} \
    ${out_paired_R1} ${out_unpaired_R1} \
    ${out_paired_R2} ${out_unpaired_R2} \
    ILLUMINACLIP:${path2adaptor}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70"
  if [ $? -ne 0 ]; then
    echo "Trimmomatic failed"
    echo "${sample} Trimmomatic UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "Uncomplete ${sample} BWA Trimmomatic `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    rm $input_R1 $input_R2
  fi
fi

#echo "## start upload clean fastq"
#${TOSCP} ${out_paired_R1} tos://nw8k-results/fastq/$sample/ >> /tmp/upload.log 2>&1
#${TOSCP} ${out_paired_R2} tos://nw8k-results/fastq/$sample/ >> /tmp/upload.log 2>&1
#echo "## clean fastq upload finished"

end_time=$(date +%s)
echo Trimmomatic End Time:`date +%Y/%m/%d--%H:%M`
((elapsed_time = $end_time - $start_time))
echo Trimmomatic Elapsed Time:$elapsed_time sec


#++++++++++++++++ Mapping +++++++++++++++++
mappingout=${WORKDIR}/03_BwaOut/${sample}
[ -e ${mappingout} ] && : || mkdir -p ${mappingout}
bwaout=${mappingout}/${sample}.bam

#input_R1=$out_paired_R1
#input_R2=$out_paired_R2

echo "#######-- BWA Start Time:`date +%Y/%m/%d--%H:%M`  --#######"
start_time=$(date +%s)

if [ ! -e "${bwaout}" ]; then
  singularity exec -B $WORKDIR:$WORKDIR  \
   ${IMGDIR}/biotools_ibp_snp_v2.sif sh -c "bwa  mem -M -K 100000000 -Y -t ${threads} \
   -R '@RG\tID:$sample\tSM:$sample\tLB:$sample\tPU:$sample\tPL:ILLUMINA' \
   ${bwa_genome} $out_paired_R1 $out_paired_R2  \
   | samtools view -bS - > ${mappingout}/${sample}.0.bam"
if [ $? -ne 0 ]; then
  echo "Mapping failed"
  echo "${sample} BWA UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  echo "Uncomplete ${sample} BWA `date +%Y/%m/%d--%H:%M`"
  exit -1;
fi
mv ${mappingout}/${sample}.0.bam ${bwaout}
fi
#mv ${mappingout}/${sample}.0.bam ${bwaout}

${TOSCP} ${bwaout} tos://nw8k-results/bams/$dataset/$sample/ >> /tmp/upload.bwa.bam.log 2>&1

echo "${sample} GATKpart1_1 Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS

# remove finished fastq.gz
find ${cleandir} -name \*.fastq.gz -exec rm -f {} \;


####+++++++++++++++ ExpansionHunter Denovo job ++++++++++++++++++

echo "#######-- ExpansionHunter Denovo Start time:`date +%Y/%m/%d--%H:%M` --#######"
[ -e ${WORKDIR}/expansionHunterDenovo/${sample} ] && : || mkdir -p ${WORKDIR}/expansionHunterDenovo/${sample}

singularity exec -B $WORKDIR:$WORKDIR  \
  $IMGDIR/ExpansionHunterDenovo.sif  \
  ExpansionHunterDenovo profile --reads ${bwaout} --reference $REFERENCE --output-prefix ${WORKDIR}/expansionHunterDenovo/${sample}/${sample}

if [ $? -ne 0 ]; then
  echo "ExpansionHunterDenovo failed"
  echo "${sample} ExpansionHunterDenovo UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  echo "Uncomplete ${sample} ExpansionHunterDenovo `date +%Y/%m/%d--%H:%M`"
  exit -1;
fi


echo "#######--ExpansionHunter Denovo Finish time:`date +%Y/%m/%d--%H:%M` --#######"

echo "#######--Download Start time:`date +%Y/%m/%d--%H:%M` --#######"

####################++++++++++++++++ Download output files +++++++++++++++++
${TOSCP} -r ${WORKDIR}/expansionHunterDenovo/${sample} tos://nw8k-results/expansionHunterDenovo/$dataset/ >> /tmp/upload.log 2>&1

echo "#######--Download Finish time:`date +%Y/%m/%d--%H:%M` --#######"
echo "${sample} ExpansionHunterDenovo Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS



wait
rm ${bwaout}
rm -rf ${WORKDIR}/expansionHunterDenovo/${sample}
end_time=$(date +%s)
echo "#######-- BWA End Time:`date +%Y/%m/%d--%H:%M`  --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- BWA Elapsed Time:$elapsed_time sec  --#######"
echo "Complete ${sample} BWA"

