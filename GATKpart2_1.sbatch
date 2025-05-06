#!/bin/bash
### Gatk Germline SNPs+Indes Pipeline, part 2: GenomicsDBImport ~~ MakeSitesOnlyVcf
### run separately per chromosomesregion
#SBATCH --get-user-env
### Job name
#SBATCH --job-name 4.surban
### Keep Output and Error
#SBATCH -o job-%x-%j.out
### Queue name
#SBATCH -p cpu6
### Specify the number of nodes and thread (ppn) for your job.
#SBATCH -n 1
##SBATCH --nodelist=MarkDupEBS-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=115G

hostname
### need to edit the next 1 line
# vcfs and samplemap list in inputdir
# input_dir=$input_dir
#sample_map.lst:  sample_name--tab--path_to_sample_vcf.gz
#samplemap=$samplemap
mem=115
threads=16
#$(nproc)
dataset=$dataset
batch=$dataset
### 所有样本的GVCF文件会比较大， 需要通过vepfs支持
### 手动将 gvcf文件下载到 vepfs 路径中， 生成 sample_map.lst
### 这里测试使用 nas 下的 /shared/data/nw8k/1KGP/
### 实际大批量跑的时候需要根据文件路径修改
#VCFDIR=/HGDP
#/shared/data/nw8k/1KGP/
#SiteOnlyDIR=/shared/data/nw8k/SiteOnlyvcf
#vfilterDIR=/shared/data/nw8k/vfilter

[ -e ${VCFDIR} ] && : || mkdir -p ${VCFDIR}

### 通过 sbatch export 参数制定 chr start end 区间
CHR=surban
sample=$CHR
#${chr}_${START}-${End}
region=${chr}:${START}-${End}
Pro=$dataset
#"LDSfam2"
WORKDIR=/HGDP
samplemap=$WORKDIR/sample_map_${Pro}.txt

#
IMGDIR="/shared/data/images"
TMPDIR="$WORKDIR/tmp"
[ -e ${TMPDIR} ] && : || mkdir -p ${TMPDIR}

# genome and annotations
ANNODIR=$WORKDIR/anno
REFERENCE=$ANNODIR/hs1/hs1.fasta

{ while true; do echo `date`  `free -g|grep Mem|awk '{print $3}'`G; sleep 300; done } &
###########################################################################################
####################++++++++++++++++ Download sample files +++++++++++++++++###############
###########################################################################################
TOSUTIL=/shared/home/chenlin/500-samples-sbatch/scripts/tosutil
TOSCP="$TOSUTIL cp -p 2 -j 5 -conf /shared/home/chenlin/500-samples-sbatch/scripts/.tosutilconfig "

#echo Download Reference Start Time:`date +%Y/%m/%d--%H:%M`
#EFDIR="/data/wgs"
[ -e ${ANNODIR} ] && : || mkdir -p ${ANNODIR}
if [ ! -e $ANNODIR/download.finished ];then
  echo Download Reference Start Time:`date +%Y/%m/%d--%H:%M`
  ${TOSCP} -r tos://bioinfo-databases/nw8k/hs1/  ${ANNODIR}  > /tmp/download.wgs.log 2>&1 # && touch /wgs/download.finished
  ${TOSCP} -r tos://bioinfo-databases/nw8k/dbSNP155 ${ANNODIR}  > /tmp/download.wgs.log 2>&1 && touch $ANNODIR/download.finished
fi

## 这里就不需要自动处理
# ${TOSCP} -r tos://nw8k-results/gvcf/$dataset/  ${VCFDIR} > /tmp/vcfs-download.log 2>&1

#echo Download Reference Finished Time:`date +%Y/%m/%d--%H:%M`

C_LOGS=/shared/data/nw8k/log/$batch/$sample.complete.log
U_LOGS=/shared/data/nw8k/log/$batch/$sample.uncomplete.log


echo "#######--  Pipeline partII Start time is `date +%Y/%m/%d--%H:%M` --#######"
work_start_time=$(date +%s)

#++++++++++++++++ GenomicsDBImport +++++++++++++++++

chrdatabasedir=${WORKDIR}/GATKprocess/05_${Pro}_DB
chrdatabase=${chrdatabasedir}/${CHR}_DB
### 测试用
[ -e ${chrdatabasedir} ] && : || rm -rf ${chrdatabasedir}
###

[ -e ${chrdatabasedir} ] && : || mkdir -p ${chrdatabasedir}
echo "#######-- chromsome_${CHR} GenomicsDBImport Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

export TILEDB_DISABLE_FILE_LOCKING=1

echo "${CHR}"

if true; then
  rm -r ${chrdatabase}
  singularity exec -B $WORKDIR:$WORKDIR  $IMGDIR/gatk4_4.5.0.sif  \
    gatk --java-options "-Xmx${mem}g -XX:ParallelGCThreads=$threads -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport --sample-name-map ${samplemap} \
    --genomicsdb-workspace-path ${chrdatabase} -L  ${region} --batch-size 100 --tmp-dir /tmp \
    --reader-threads 4 --genomicsdb-shared-posixfs-optimizations

  if [ $? -ne 0 ]; then
    echo "GenomicsDBImport failed"
    echo "${sample} GenomicsDBImport UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "Uncomplete ${sample} GenomicsDBImport `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    echo "${sample} GenomicsDBImport Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${sample} GenomicsDBImport Complete `date +%Y/%m/%d--%H:%M`"
  fi
fi

end_time=$(date +%s)
echo "#######-- chromsome_${CHR} GenomicsDBImport End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- chromsome_${CHR} GenomicsDBImport Elapsed Time:$elapsed_time sec --#######"

#++++++++++++++++ GenotypeGVCFs +++++++++++++++++
cohortout=${WORKDIR}/GATKprocess/06_${Pro}_Cohort
[ -e ${cohortout} ] && : || mkdir -p ${cohortout}
echo "#######-- ${CHR} GenotypeGVCFs Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

if true; then
  singularity exec -B $WORKDIR:$WORKDIR $IMGDIR/gatk4_4.5.0.sif \
    sh -c " gatk --java-options '-Xmx${mem}g -XX:ParallelGCThreads=$threads -Djava.io.tmpdir=/tmp' GenotypeGVCFs \
    -R ${REFERENCE} -D ${ANNODIR}/dbSNP155/hs1_dbSNPv155.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation \
    -V gendb://${chrdatabase} -O ${cohortout}/${Pro}_${CHR}.vcf.gz --tmp-dir /tmp \
    --genomicsdb-shared-posixfs-optimizations --seconds-between-progress-updates 300"

  if [ $? -ne 0 ]; then
    echo "GenotypeGVCFs failed"
    echo "${CHR} GenotypeGVCFs UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${CHR} GenotypeGVCFs UnComplete `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    echo "${CHR} GenotypeGVCFs Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} GenotypeGVCFs Complete `date +%Y/%m/%d--%H:%M`"
  fi

fi

#${TOSCP} -r ${cohortout}/${Pro}_${CHR}.vcf.gz tos://1kg-lius/GATK4-JointCalling/ >> /tmp/upload.log 2>&1
#${TOSCP} -r ${cohortout}/${Pro}_${CHR}.vcf.gz.tbi tos://1kg-lius/GATK4-JointCalling/ >> /tmp/upload.log 2>&1

end_time=$(date +%s)
echo "#######-- ${CHR} GenotypeGVCFs End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- ${CHR} GenotypeGVCFs Elapsed Time:$elapsed_time sec --#######"


#++++++++++++++++ VariantFiltration +++++++++++++++++
echo "#######-- ${CHR} VariantFiltration Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

if true; then
  singularity exec -B $WORKDIR:$WORKDIR $IMGDIR/gatk4_4.5.0.sif \
    sh -c " gatk VariantFiltration --java-options '-XX:ParallelGCThreads=$threads -Djava.io.tmpdir=/tmp' \
    --filter-expression 'ExcessHet > 54.69' --filter-name ExcessHet \
    -V ${cohortout}/${Pro}_${CHR}.vcf.gz -O ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz"

  if [ $? -ne 0 ]; then
    echo "VariantFiltration failed"
    echo "${CHR} VariantFiltration UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${CHR} VariantFiltration UnComplete `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    echo "${CHR} VariantFiltration Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} VariantFiltration Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP}  ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz tos://nw8k-results/vfilter/$dataset/ >> /tmp/upload.log 2>&1

  if [ $? -ne 0 ]; then
    echo "vfiltered Upload failed"
    echo "${sample} vfiltered_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${sample} vfiltered_Upload UnComplete `date +%Y/%m/%d--%H:%M`"
    #exit -1;
  else
    #echo "${CHR} vfiltered_Upload Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} vfiltered_Upload Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP}  ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz.tbi tos://nw8k-results/vfilter/$dataset/ >> /tmp/upload.log 2>&1

fi
#[ -e ${vfilterDIR} ] && : || mkdir -p ${vfilterDIR}/
#cp ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz ${vfilterDIR}/
#cp ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz.tbi ${vfilterDIR}/

end_time=$(date +%s)
echo "#######-- ${CHR} VariantFiltration End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- ${CHR} VariantFiltration Elapsed Time:$elapsed_time sec --#######"


#++++++++++++++++ MakeSitesOnlyVcf +++++++++++++++++
echo "#######-- ${CHR} MakeSitesOnlyVcf Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

if true; then
  singularity exec -B $WORKDIR:$WORKDIR $IMGDIR/gatk4_4.5.0.sif \
    sh -c " gatk MakeSitesOnlyVcf --java-options '-XX:ParallelGCThreads=$threads -Djava.io.tmpdir=/tmp' \
    -I ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz -O ${cohortout}/${Pro}_${CHR}_SiteOnly.vcf.gz "

  if [ $? -ne 0 ]; then
    echo "MakeSitesOnlyVcf failed"
    echo "${CHR} MakeSitesOnlyVcf UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${CHR} MakeSitesOnlyVcf UnComplete `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    echo "${CHR} MakeSitesOnlyVcf Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} MakeSitesOnlyVcf Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP}  ${cohortout}/${Pro}_${CHR}_SiteOnly.vcf.gz tos://nw8k-results/SiteOnly/$dataset/ >> /tmp/upload.log 2>&1

  if [ $? -ne 0 ]; then
    echo "SiteOnly Upload failed"
    echo "${sample} SiteOnly_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${sample} SiteOnly_Upload UnComplete `date +%Y/%m/%d--%H:%M`"
    #exit -1;
  else
    #echo "${CHR} SiteOnly_Upload Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} SiteOnly_Upload Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP}  ${cohortout}/${Pro}_${CHR}_SiteOnly.vcf.gz.tbi tos://nw8k-results/SiteOnly/$dataset/ >> /tmp/upload.log 2>&1

fi

#[ -e ${SiteOnlyDIR} ] && : || mkdir -p ${SiteOnlyDIR}/
#cp ${cohortout}/${Pro}_${CHR}_SiteOnly.vcf.gz ${SiteOnlyDIR}/
#cp ${cohortout}/${Pro}_${CHR}_SiteOnly.vcf.gz.tbi ${SiteOnlyDIR}/


end_time=$(date +%s)
echo "#######-- ${CHR} MakeSitesOnlyVcf End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- ${CHR} MakeSitesOnlyVcf Elapsed Time:$elapsed_time sec --#######"

#++++++++++++++++ Ending +++++++++++++++++
work_end_time=$(date +%s)
echo "#######-- Gatk4 Pipeline partII End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $work_end_time - $work_start_time))
echo "#######-- $Pro GATK4 Cohort PartII Elapsed Time:$elapsed_time sec --#######"
echo "GenotypeGVCFs ${CHR} Complete"

#rm -rf ${chrdatabasedir}
#rm -rf ${cohortout}
