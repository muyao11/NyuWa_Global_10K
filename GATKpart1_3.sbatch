#!/bin/bash
### Gatk4 Germline SNPs+Indes Pipeline, part I: Data pre-processing (single sample)
###
### Inherit all current environment variables
#SBATCH --get-user-env
### Job name
#SBATCH --job-name surban.SNP.gatk4-Part1-persample
### Keep Output and Error
#SBATCH -o job-%x-%j.out
### Queue name
#SBATCH -p C16G64
### Specify the number of nodes and thread (ppn) for your job.
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=56G

echo Pipeline partI Start time is `date +%Y/%m/%d--%H:%M`
work_start_time=$(date +%s)



# Sample ID
sample=$sample
WORKDIR=/data
#INPUTDIR=$input_dir
dataset=$dataset
threads=$(nproc)
mem=64


# Input files
## get input reads from command line parameters
markbam=$markbam

TOSUTIL=/shared/home/chenlin/500-samples-sbatch/scripts/tosutil
TOSCP="$TOSUTIL cp -p 2 -j 5 -conf /shared/home/chenlin/500-samples-sbatch/scripts/.tosutilconfig "
IMGDIR="/shared/data/images"
TMPDIR="/data/tmp"




[ -e $TMPDIR ] && : || mkdir -p $TMPDIR


# genome and annotations
ANNODIR=/wgs/anno/hs1/v0
REFERENCE=${ANNODIR}/hs1/hs1.fasta
bwa_genome=${ANNODIR}/hs1/hs1.fasta
path2adaptor=${ANNODIR}/hs1/TruSeq3-PE-2.fa


[ -e ${ANNODIR} ] && : || mkdir -p ${ANNODIR}
echo Download Start Time:`date +%Y/%m/%d--%H:%M`
if [ ! -e /wgs/download.hs1.finished ];then
    ${TOSCP} -r tos://bioinfo-databases/nw8k/hs1/  ${ANNODIR}  > /tmp/download.wgs.log 2>&1 && touch /wgs/download.hs1.finished
fi




# Log files
[ -e $WORKDIR/output ] && : || mkdir -p $WORKDIR/output
C_LOGS=$WORKDIR/output/log_GATK4_PartI_complete.lst
U_LOGS=$WORKDIR/output/log_GATK4_PartI_uncomplete.lst

####++++++++++++++++ Download Reads +++++++++++++++++ 
#
#TOSUTIL=/shared/home/chenlin/500-samples-sbatch/scripts/tosutil 
#TOSCP="$TOSUTIL cp -p 2 -j 5 -conf /shared/home/chenlin/500-samples-sbatch/scripts/.tosutilconfig "

gatkoutdir=${WORKDIR}/04_GatkOut/${sample}
[ -e ${gatkoutdir} ] && : || mkdir -p ${gatkoutdir}
gatkout=${gatkoutdir}/${sample}
echo "${TOSCP} $markbam ${gatkoutdir}/ >  /tmp/markbamdownload.log 2>&1"
${TOSCP} -r ${markbam} ${gatkoutdir}/ >  /tmp/markbamdownload.log 2>&1


if [ ! -e ${gatkout}.Mdup.bam.bai ];then
  singularity exec -B $input_dir:$input_dir -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/tmp \
    ${IMGDIR}/biotools_ibp_snp_v2.sif \
    sh -c "samtools index ${gatkout}.Mdup.bam"
fi

#++++++++++++++++ InDel Realigner by gatk3 +++++++++++++++++
echo "##########-- RealignerTargetCreator Start Time:`date +%Y/%m/%d--%H:%M` --#########"
start_time=$(date +%s)
if [ ! -e ${gatkout}.realignertargetcreator.intervals ];then
  singularity exec -B $input_dir:$input_dir -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/tmp \
    ${IMGDIR}/gatk3_3.8-0.sif \
    sh -c "java -Djava.io.tmpdir=/tmp -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
    -nt ${threads} -R ${REFERENCE} -o ${gatkout}.realignertargetcreator.intervals \
    -I ${gatkout}.Mdup.bam -known ${ANNODIR}/hs1/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.hs1.vcf.gz "
  end_time=$(date +%s)
  if [ $? -ne 0 ]; then
    echo "RealignerTargetCreator failed"
    echo "${sample} UnComplete RealignerTargetCreator `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    exit -1;
  fi
fi


echo "##########-- RealignerTargetCreator End Time:`date +%Y/%m/%d--%H:%M` --#########"
((elapsed_time = $end_time - $start_time))
echo "##########-- RealignerTargetCreator Elapsed Time:$elapsed_time sec --#########"

echo "##########-- IndelRealigner Start Time:`date +%Y/%m/%d--%H:%M` --#########"
start_time=$(date +%s)
singularity exec  -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/tmp \
    ${IMGDIR}/gatk3_3.8-0.sif \
    sh -c "java -Djava.io.tmpdir=/tmp -jar /usr/GenomeAnalysisTK.jar  -T IndelRealigner -I ${gatkout}.Mdup.bam \
    -R ${REFERENCE} -targetIntervals ${gatkout}.realignertargetcreator.intervals \
    -o ${gatkout}.Mdup.realigned.bam -known ${ANNODIR}/hs1/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.hs1.vcf.gz "
if [ $? -ne 0 ]; then
  echo "IndelRealigner failed"
  echo "${sample} UnComplete IndelRealigner `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  exit -1;
fi

## remove 
if [ -f ${gatkout}.Mdup.bam ];then
    rm -f ${gatkout}.Mdup.bam 
fi
end_time=$(date +%s)
echo "##########-- IndelRealigner End Time:`date +%Y/%m/%d--%H:%M` --#########"
((elapsed_time = $end_time - $start_time))
echo "##########-- IndelRealigner Elapsed Time:$elapsed_time sec --#########"


#++++++++++++++++ BaseRecalibrator +++++++++++++++++
echo "#######-- BaseRecalibrator Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)
singularity exec  -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/tmp \
    ${IMGDIR}/gatk4_4.5.0.sif \
    gatk  BaseRecalibrator --java-options "-Xms40G -Xmx40G -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=/tmp/" \
    -R ${REFERENCE} --known-sites ${ANNODIR}/hs1/hs1_dbSNPv155.vcf.gz \
    --known-sites ${ANNODIR}/hs1/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.hs1.vcf.gz \
    --known-sites ${ANNODIR}/hs1/resources-broad-hg38-v0-Homo_sapiens_assembly38.known_indels.hs1.vcf.gz \
    -O ${gatkout}.recalibration.table -I ${gatkout}.Mdup.realigned.bam \
    --tmp-dir /tmp/ --seconds-between-progress-updates 300
if [ $? -ne 0 ]; then
  echo "BaseRecalibrator failed"
  echo "${sample} UnComplete BaseRecalibrator `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  exit -1;
fi

end_time=$(date +%s)
echo "#######-- BaseRecalibrator End Time:`date +%Y/%m/%d--%H:%M`  --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- BaseRecalibrator Elapsed Time:$elapsed_time sec --#######"

#++++++++++++++++ ApplyBQSR +++++++++++++++++
echo "#######-- ApplyBQSR Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)
singularity exec  -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/tmp \
    ${IMGDIR}/gatk4_4.5.0.sif \
    gatk ApplyBQSR --java-options "-Xms40G -Xmx40G -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=/tmp/"  \
    -R ${REFERENCE} -I ${gatkout}.Mdup.realigned.bam -bqsr ${gatkout}.recalibration.table \
    -O ${gatkout}.Mdup.realigned.recal.bam --tmp-dir /tmp/ --static-quantized-quals 10 \
    --static-quantized-quals 20 --static-quantized-quals 30 --seconds-between-progress-updates 300
if [ $? -ne 0 ]; then
  echo "ApplyBQSR failed"
  echo "${sample} UnComplete ApplyBQSR `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  exit -1;
fi

if [ -f ${gatkout}.Mdup.realigned.bam ];then
    rm -f ${gatkout}.Mdup.realigned.bam
fi 

end_time=$(date +%s)
echo "#######-- ApplyBQSR End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- ApplyBQSR Elapsed Time:$elapsed_time sec --#######"


echo "#######-- !!!!!!!!!!!!!!!!!!! Mdup.realigned.recal.bam have done !!!!!!!!!!!!!! --#######"
#++++++++++++++++ Compress Bam +++++++++++++++++
echo "#######-- Compress Bam Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

if [ ! -f "${gatkout}.final.bam" ]; then
singularity exec -B $input_dir:$input_dir -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/tmp \
    ${IMGDIR}/biotools_ibp_snp_v2.sif \
    sh -c "samtools merge -l 5 --threads ${threads} ${gatkout}.final.bam ${gatkout}.Mdup.realigned.recal.bam"
if [ $? -ne 0 ]; then
  echo "Compress Bam failed"
  echo "${sample} UnComplete CompressBam `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  echo "Uncomplete ${sample} CompressBam `date +%Y/%m/%d--%H:%M`"
  exit -1;
fi
if [ ! -e ${gatkout}.final.bam.bai ];then
  singularity exec -B $input_dir:$input_dir -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/tmp \
    ${IMGDIR}/biotools_ibp_snp_v2.sif \
    sh -c "samtools index ${gatkout}.final.bam"
if [ $? -ne 0 ]; then
  echo "Compress Bam Index failed"
  echo "${sample} UnComplete CompressIndexBam `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  echo "Uncomplete ${sample} CompressIndexBam `date +%Y/%m/%d--%H:%M`"
  exit -1;
fi
if [ -f ${gatkout}.Mdup.realigned.bam ];then
    rm -f ${gatkout}.Mdup.realigned.bam
fi 
if [ -f ${gatkout}.Mdup.realigned.recal.bam ];then
    rm -f ${gatkout}.Mdup.realigned.recal.bam
fi 
#if [ -f ${WORKDIR}/03_BwaOut/${sample}/${sample}.bam ];then
#    rm -f ${WORKDIR}/03_BwaOut/${sample}/${sample}.bam
#fi


${TOSCP} ${gatkout}.final.bam tos://nw8k-results/finalbam/$sample/ > /tmp/upload.log 2>&1
${TOSCP} ${gatkout}.final.bam.bai tos://nw8k-results/finalbam/$sample/ > /tmp/upload.log 2>&1

echo "${sample} Complete GATKpart1_3 `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
end_time=$(date +%s)
echo "#######-- Compress Bam End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- Compress Bam Elapsed Time:$elapsed_time sec --#######"

#++++++++++++++++ HaplotypeCaller +++++++++++++++++
#echo "#######-- HaplotypeCaller Start Time:`date +%Y/%m/%d--%H:%M` --#######"
#start_time=$(date +%s)
#singularity exec  -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/tmp \
#    ${IMGDIR}/gatk4_4.5.0.sif \
#    gatk HaplotypeCaller --java-options "-Xms40G -Xmx40G -XX:ParallelGCThreads=${threads} -Djava.io.tmpdir=/tmp/" \
#    -R ${REFERENCE} -ERC GVCF -I ${gatkout}.final.bam \
#    -O ${gatkout}.raw.realign.snps.indels.g.vcf.gz  --tmp-dir /tmp/ \
#    --native-pair-hmm-threads ${threads} --seconds-between-progress-updates 300
#if [ $? -ne 0 ]; then
#  echo "HaplotypeCaller failed"
#  echo "${sample} UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
#  exit -1;
#fi

#${TOSCP} ${gatkout}.raw.realign.snps.indels.g.vcf.gz  tos://nw8k-results/gvcf/$dataset/ > /tmp/upload.log 2>&1
#${TOSCP} ${gatkout}.raw.realign.snps.indels.g.vcf.gz.tbi  tos://nw8k-results/gvcf/$dataset/ > /tmp/upload.log 2>&1

#end_time=$(date +%s)
#echo "#######-- HaplotypeCaller End Time:`date +%Y/%m/%d--%H:%M` --#######"
#((elapsed_time = $end_time - $start_time))
#echo "#######-- HaplotypeCaller Elapsed Time:$elapsed_time sec --#######"

#++++++++++++++++ MappingResult QC +++++++++++++++++
unset DISPLAY
bamQCout=${WORKDIR}/GATKprocess/02_QC/QC_map/${sample}
[ -e ${bamQCout} ] && : || mkdir -p ${bamQCout}
echo "#######-- QC_map Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

if [ ! -f "${bamQCout}/genome_results.txt" ]; then
singularity exec -B $WORKDIR:$WORKDIR -B /wgs:/wgs -B ${TMPDIR}:/datatmp \
  /shared/data/images/qctools_v1.1.sif /bin/bash -c "
  qualimap bamqc -bam ${gatkout}.final.bam -c --java-mem-size=40G -nt ${threads} -outdir ${bamQCout} --skip-duplicated
  "
fi
${TOSCP} -r ${bamQCout}/  tos://nw8k-results/bamQC/$sample/ > /tmp/upload.log 2>&1

end_time=$(date +%s)
echo "#######-- QC_map End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- QC_map Elapsed Time:$elapsed_time sec --#######"



if [ -f $WORKDIR ];then
  rm -rf $WORKDIR/*
fi

work_end_time=$(date +%s)
echo "#######-- GATK4 Pipeline partI End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $work_end_time - $work_start_time))
echo "#######-- $sample GATK4 Cohort PartI Elapsed Time:$elapsed_time sec --#######"
