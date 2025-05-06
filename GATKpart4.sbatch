#!/bin/bash
### Gatk4 Germline SNPs+Indes Pipeline, part 4: VQSR ApplyVQSR
### run separately per chromosomes
#SBATCH --get-user-env
### Job name
#SBATCH --job-name 7.surban
### Keep Output and Error
#SBATCH -o job-%x-%j.out
### Queue name
#SBATCH -p cpu6
### Specify the number of nodes and thread (ppn) for your job.
#SBATCH -n 1
##SBATCH --nodelist=MarkDupEBS-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G


### need to edit the next 1 line

CHR=surban
sample=$CHR
WORKDIR=/data
dataset=$dataset
batch=$dataset
threads=$(nproc)
Pro=$dataset


IMGDIR="/shared/data/images"
TMPDIR="/data/tmp"
[ -e ${TMPDIR} ] && : || mkdir -p ${TMPDIR}

cohortout=${WORKDIR}/GATKprocess/06_${Pro}_Cohort
[ -e ${cohortout} ] && : || mkdir -p ${cohortout}

# genome and annotations
ANNODIR=$WORKDIR/anno
REFERENCE=$ANNODIR/hs1/hs1.fasta


C_LOGS=/shared/data/nw8k/log/$batch/$sample.complete.log
U_LOGS=/shared/data/nw8k/log/$batch/$sample.uncomplete.log


###########################################################################################
####################++++++++++++++++ Download sample files +++++++++++++++++############### 
###########################################################################################
TOSUTIL=/shared/home/chenlin/500-samples-sbatch/scripts/tosutil 
TOSCP="$TOSUTIL cp -p 2 -j 5 -conf /shared/home/chenlin/500-samples-sbatch/scripts/.tosutilconfig "

#echo Download Reference Start Time:`date +%Y/%m/%d--%H:%M`
#[ -e ${ANNODIR} ] && : || mkdir -p ${ANNODIR}
#if [ ! -e /wgs/download.finished ];then
#  ${TOSCP} -r tos://bioinfo-databases/nw8k/hs1/  ${ANNODIR}  > /tmp/download.wgs.log 2>&1 && touch /wgs/download.finished
#fi

############################ Download SiteOnly.vcf.gz  ############################
#${TOSCP} -r tos://nw8k-results/GATK4-P3/$dataset/06_Cohort/  $(dirname ${cohortout})/  > /tmp/download.vcf.log 2>&1

#mv ${WORKDIR}/$dataset/*  ${cohortout}
#mv ${cohortout}/06_Cohort/* ${cohortout}/
#mv ${cohortout}/$chr/* ${cohortout}

#${TOSCP} -r tos://1kg-lius/GATK4-P3/ ${WORKDIR}/  > /tmp/download.vcf.log 2>&1
#mv ${WORKDIR}/GATK4-P3/*  ${cohortout}

#echo Download Reference Finished Time:`date +%Y/%m/%d--%H:%M`



echo "#######--  Pipeline partIV Start time is `date +%Y/%m/%d--%H:%M` --#######"
work_start_time=$(date +%s)


#++++++++++++++++ ApplyVQSR for snp +++++++++++++++++
echo "#######-- chromsome_${CHR} ApplyVQSR SNP Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

if true; then
  singularity exec -B $WORKDIR:$WORKDIR  ${IMGDIR}/gatk4_4.5.0.sif \
    gatk  --java-options "-Xmx100g -XX:ParallelGCThreads=$threads" ApplyVQSR \
    -mode SNP -R ${REFERENCE} --truth-sensitivity-filter-level 99.8  \
    --tranches-file ${cohortout}/${Pro}_recalibrate_SNP.tranches  \
    --recal-file ${cohortout}/${Pro}_recalibrate_SNP.recal \
    -V ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz  \
    -O ${cohortout}/${Pro}_${CHR}.recalibrate_SNP.vcf.gz

  if [ $? -ne 0 ]; then
    echo "ApplyVQSR SNP failed"
    echo "${CHR} ApplyVQSR_SNP UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${CHR} ApplyVQSR_SNP UnComplete `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    #echo "${CHR} GenotypeGVCFs Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} ApplyVQSR_SNP Complete `date +%Y/%m/%d--%H:%M`"
  fi

fi
end_time=$(date +%s)
echo "#######-- chromsome_${CHR} ApplyVQSR SNP End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- chromsome_${CHR} ApplyVQSR SNP Elapsed Time:$elapsed_time sec --#######"


#++++++++++++++++ ApplyVQSR for INDEL +++++++++++++++++
echo "#######-- chromsome_${CHR} ApplyVQSR INDEL Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

#20200309- run with truth-sensitivity-filter-levle = 99.0
if true; then
  singularity exec -B $WORKDIR:$WORKDIR ${IMGDIR}/gatk4_4.5.0.sif \
    gatk  --java-options "-Xmx100g -XX:ParallelGCThreads=$threads" ApplyVQSR -mode INDEL \
    -R ${REFERENCE} --truth-sensitivity-filter-level 99.0  \
    --tranches-file ${cohortout}/${Pro}_recalibrate_INDEL.tranches  \
    --recal-file ${cohortout}/${Pro}_recalibrate_INDEL.recal \
    -V ${cohortout}/${Pro}_${CHR}.recalibrate_SNP.vcf.gz \
    -O ${cohortout}/${Pro}_${CHR}.recalibrate_SNP_INDEL.vcf.gz

  if [ $? -ne 0 ]; then
    echo "ApplyVQSR INDEL failed"
    echo "${CHR} ApplyVQSR_INDEL UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${CHR} ApplyVQSR_INDEL UnComplete `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    #echo "${CHR} GenotypeGVCFs Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} ApplyVQSR_INDEL Complete `date +%Y/%m/%d--%H:%M`"
  fi

 {
  ${TOSCP}  ${cohortout}/${Pro}_${CHR}.recalibrate_SNP_INDEL.vcf.gz tos://nw8k-results/final_recal/$dataset/ >> /tmp/upload.log 2>&1

  if [ $? -ne 0 ]; then
    echo "VQSR Upload failed"
    echo "${sample} VQSR_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${sample} VQSR_Upload UnComplete `date +%Y/%m/%d--%H:%M`"
    #exit -1;
  else
    #echo "${CHR} vfiltered_Upload Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} VQSR_Upload Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP}  ${cohortout}/${Pro}_${CHR}.recalibrate_SNP_INDEL.vcf.gz.tbi tos://nw8k-results/final_recal/$dataset/ >> /tmp/upload.log 2>&1
 } &

fi
end_time=$(date +%s)
echo "#######-- chromsome_${CHR} ApplyVQSR INDEL End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- chromsome_${CHR} ApplyVQSR INDEL Elapsed Time:$elapsed_time sec --#######"

#++++++++++++++++ MakeSitesOnlyVcf +++++++++++++++++

if true; then
 {
  echo "#######-- chr${CHR} MakeSitesOnlyVcf Start Time:`date +%Y/%m/%d--%H:%M` --#######"
  start_time=$(date +%s)
  singularity exec -B $WORKDIR:$WORKDIR ${IMGDIR}/gatk4_4.5.0.sif gatk MakeSitesOnlyVcf  \
    -I ${cohortout}/${Pro}_${CHR}.recalibrate_SNP_INDEL.vcf.gz \
    -O ${cohortout}/${Pro}_${CHR}.recalibrate_SNP_INDEL_SiteOnly.vcf.gz

  if [ $? -ne 0 ]; then
    echo "VQSR SitesOnly Upload failed"
    echo "${sample} VQSR_SitesOnly UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${sample} VQSR_SitesOnly UnComplete `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    #echo "${CHR} vfiltered_Upload Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} VQSR_SitesOnly Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP}  ${cohortout}/${Pro}_${CHR}.recalibrate_SNP_INDEL_SiteOnly.vcf.gz tos://nw8k-results/final_recal/$dataset/ >> /tmp/upload.log 2>&1

  if [ $? -ne 0 ]; then
    echo "VQSR SitesOnly Upload failed"
    echo "${sample} VQSR_SitesOnly_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${sample} VQSR_SitesOnly_Upload UnComplete `date +%Y/%m/%d--%H:%M`"
    #exit -1;
  else
    echo "${CHR} VQSR Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} VQSR_SitesOnly_Upload Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP}  ${cohortout}/${Pro}_${CHR}.recalibrate_SNP_INDEL_SiteOnly.vcf.gz.tbi tos://nw8k-results/final_recal/$dataset/ >> /tmp/upload.log 2>&1

  end_time=$(date +%s)
  echo "#######-- chr${CHR} MakeSitesOnlyVcf End Time:`date +%Y/%m/%d--%H:%M` --#######"
  ((elapsed_time = $end_time - $start_time))
  echo "#######-- chr${CHR} MakeSitesOnlyVcf Elapsed Time:$elapsed_time sec --#######"
 } &
fi

if true; then

  echo "#######-- chr${CHR} MakePassVcf Start Time:`date +%Y/%m/%d--%H:%M` --#######"
  singularity exec -B $WORKDIR:$WORKDIR  ${IMGDIR}/biotools_ibp_v2.sif \
  sh -c " bcftools view -f PASS --threads 5 -Oz -o ${cohortout}/${Pro}_${CHR}.PASS.vcf.gz ${cohortout}/${Pro}_${CHR}.recalibrate_SNP_INDEL.vcf.gz  " 
  #&& tabix ${cohortout}/${Pro}_${CHR}.PASS.vcf.gz "

  if [ $? -ne 0 ]; then
    echo "Filter_PASS failed"
    echo "${CHR} Filter_PASS UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${CHR} Filter_PASS UnComplete `date +%Y/%m/%d--%H:%M`"
    exit -1;
  fi

 {
  ${TOSCP}  ${cohortout}/${Pro}_${CHR}.PASS.vcf.gz tos://nw8k-results/final_pass/$dataset/ >> /tmp/upload.log 2>&1

  if [ $? -ne 0 ]; then
    echo "Filter_PASS Upload failed"
    echo "${sample} Filter_PASS_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${sample} Filter_PASS_Upload UnComplete `date +%Y/%m/%d--%H:%M`"
    #exit -1;
  else
    echo "${CHR} Filter_PASS Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} Filter_PASS_Upload Complete `date +%Y/%m/%d--%H:%M`"
  fi
 } &

  singularity exec -B $WORKDIR:$WORKDIR  ${IMGDIR}/biotools_ibp_v2.sif \
    tabix ${cohortout}/${Pro}_${CHR}.PASS.vcf.gz

  ${TOSCP}  ${cohortout}/${Pro}_${CHR}.PASS.vcf.gz.tbi tos://nw8k-results/final_pass/$dataset/ >> /tmp/upload.log 2>&1

fi

wait

#rsync -avvR ${cohortout}/./${Pro}_${CHR}.recalibrate_SNP* ${OUTDIR}/$dataset/ 
#${TOSCP} -r ${OUTDIR}/ tos://nw8k-results/ > /tmp/gatk4-p4-upload.log   2>&1 


#++++++++++++++++ Annotation +++++++++++++++++
#echo "#######-- Annotation Start Time:`date +%Y/%m/%d--%H:%M` --#######"
#start_time=$(date +%s)

#$annovar/table_annovar.pl --vcfinput ${cohortout}/${Pro1}_chr${CHR}.recalibrate_SNP_INDEL_SiteOnly.vcf.gz  $annovar/humandb_hg38/ -buildver hg38 -out ${cohortout}/${Pro1}_chr${CHR}.recalibrate_SNP_INDEL.annovar.recal -remove -protocol refGene,ensGene,NONCODEv5,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,1000g2015aug_amr,1000g2015aug_afr,gnomad211_genome,gnomad30_genome,clinvar_20210123,exac03,avsnp150,dbnsfp41a,icgc28 -operation g,g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -polish --otherinfo

#$annovar/table_annovar.pl --vcfinput ${cohortout}/${Pro}_chr${CHR}.recalibrate_INDEL.vcf.gz  $annovar/humandb_hg38/ -buildver hg38 -out ${cohortout}/${Pro}_chr${CHR}.recalibrate_INDEL.annovar.recal -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,gnomad_genome,clinvar_20170130,exac03,avsnp147,dbnsfp30a -operation g,f,f,f,f,f,f,f -nastring . -polish --otherinfo
#bgzip ${cohortout}/${Pro1}_chr${CHR}.recalibrate_SNP_INDEL.annovar.recal.hg38_multianno.vcf
#tabix ${cohortout}/${Pro1}_chr${CHR}.recalibrate_SNP_INDEL.annovar.recal.hg38_multianno.vcf.gz

#end_time=$(date +%s)
#echo "#######-- Annotation End Time:`date +%Y/%m/%d--%H:%M` --#######"
#((elapsed_time = $end_time - $start_time))
#echo "#######-- Annotation Elapsed Time:$elapsed_time sec --#######"


#++++++++++++++++ Ending +++++++++++++++++
work_end_time=$(date +%s)
echo "#######-- Gatk4 Pipeline partIV End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $work_end_time - $work_start_time))
echo "#######-- $Pro GATK4 Cohort PartIV Elapsed Time:$elapsed_time sec --#######"
echo "chr_${CHR} Complete `date +%Y/%m/%d--%H:%M`" 
