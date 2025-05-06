#!/bin/bash
### Gatk4 Germline SNPs+Indes Pipeline, part 3: VQSR VariantRecalibrator snp
### run separately per chromosomes
#SBATCH --get-user-env
### Job name
#SBATCH --job-name 6.surban
### Keep Output and Error
#SBATCH -o job-%x-%j.out
### Queue name
#SBATCH -p cpu1
### Specify the number of nodes and thread (ppn) for your job.
#SBATCH -n 1
##SBATCH --nodelist=
#SBATCH --cpus-per-task=16
#SBATCH --mem=230G


WORKDIR=/HGDP
dataset=$dataset
batch=$dataset
threads=$(nproc)
Pro=$dataset
sample=$dataset

### 
IMGDIR="/shared/data/images"
TMPDIR="/data/tmp"
[ -e ${TMPDIR} ] && : || mkdir -p ${TMPDIR}


### need to edit the next 1 line
cohortout=${WORKDIR}/GATKprocess/06_${Pro}_Cohort
[ -e ${cohortout} ] && : || mkdir -p ${cohortout}

# genome and annotations
ANNODIR=$WORKDIR/anno
REFERENCE=$ANNODIR/hs1/hs1.fasta

axiomPoly_resource_vcf=${ANNODIR}/hs1/resources-broad-hg38-v0-Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.hs1.vcf.gz

U_LOGS=/shared/data/nw8k/log/$batch/$sample.uncomplete.log
C_LOGS=/shared/data/nw8k/log/$batch/$sample.complete.log


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
#${TOSCP} -r tos://nw8k-results/final_vcf/$dataset/ ${WORKDIR}/  > /tmp/download.vcf.log 2>&1
#${TOSCP} -r tos://nw8k-results/final_vcf/$dataset/ ${WORKDIR}/  > /tmp/download.vcf.log 2>&1
#mv ${WORKDIR}/$dataset/*  ${cohortout}

#singularity exec -B $cohortout:$cohortout -B /data/wgs:/wgs -B ${TMPDIR}:/datatmp \
# /shared/data/images/biotools_ibp_snp_v2.sif sh -c "
#for f in `ls ${cohortout}/*.vcf.gz`
#do
#    tabix $f
#    ${TOSCP} $f.tbi tos://1kg-lius/GATK4-JointCalling/ > /tmp/upload.vcf.log 2>&1
#done
#"

#echo Download Reference Finished Time:`date +%Y/%m/%d--%H:%M`

echo "#######--  Pipeline partIII snp Start time is `date +%Y/%m/%d--%H:%M` --#######"
work_start_time=$(date +%s)


#++++++++++++++++ VariantRecalibrator for snp +++++++++++++++++
echo "#######-- ${Pro} VariantRecalibrator SNP Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)

#for CHR in `ls ${cohortout}/*_SiteOnly.vcf.gz`
#for CHR in `find ${cohortout}/ -type f -name \*_SiteOnly.vcf.gz`
for CHR in {1..22} X Y
do
    inputchr=" ${inputchr} -V ${cohortout}/${Pro}_chr${CHR}_SiteOnly.vcf.gz "
done

## SNP
if true;then
{
  singularity exec -B $WORKDIR:$WORKDIR  $IMGDIR/gatk4_4.5.0.sif  sh -c "gatk \
    --java-options '-Xmx150g -XX:ParallelGCThreads=16' VariantRecalibrator -mode SNP -R ${REFERENCE} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ANNODIR/hs1/resources-broad-hg38-v0-hapmap_3.3.hg38.hs1.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 $ANNODIR/hs1/resources-broad-hg38-v0-1000G_omni2.5.hg38.hs1.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $ANNODIR/hs1/resources-broad-hg38-v0-1000G_phase1.snps.high_confidence.hg38.hs1.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7 $ANNODIR/hs1/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.hs1.vcf.gz \
    -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 \
    -tranche 99.5 -tranche 99.4 -tranche 99.3  -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    --trust-all-polymorphic --max-gaussians 6  ${inputchr} -O ${cohortout}/${Pro}_recalibrate_SNP.recal \
    --tranches-file ${cohortout}/${Pro}_recalibrate_SNP.tranches \
    --rscript-file ${cohortout}/${Pro}_recalibrate_SNP.plots.R" 

  if [ $? -ne 0 ]; then
    echo "VariantRecalibrator_SNP failed"
    echo "${sample} VariantRecalibrator_SNP UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1
    #comp_vqsr=0
  else
    echo "${sample} VariantRecalibrator_SNP Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
  fi

  ${TOSCP} ${cohortout}/${Pro}_recalibrate_SNP.recal tos://nw8k-results/vqsr/$dataset/ > /tmp/tos-gatk-p3-upload.log 2>&1

  if [ $? -ne 0 ]; then
    echo "VariantRecalibrator_SNP Upload failed"
    echo "${sample} VariantRecalibrator_SNP_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${sample} VariantRecalibrator_SNP_Upload UnComplete `date +%Y/%m/%d--%H:%M`"
    #exit -1;
  else
    #echo "${CHR} SiteOnly_Upload Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} VariantRecalibrator_SNP_Upload Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP} ${cohortout}/${Pro}_recalibrate_SNP.tranches tos://nw8k-results/vqsr/$dataset/ > /tmp/tos-gatk-p3-upload.log 2>&1
  ${TOSCP} ${cohortout}/${Pro}_recalibrate_SNP.plots.R tos://nw8k-results/vqsr/$dataset/ > /tmp/tos-gatk-p3-upload.log 2>&1
}&

fi

## INDEL
if true;then
{
  singularity exec -B $WORKDIR:$WORKDIR  $IMGDIR/gatk4_4.5.0.sif sh -c "gatk \
    --java-options '-Xmx120g -XX:ParallelGCThreads=16' VariantRecalibrator \
    -mode INDEL -R ${REFERENCE} --max-gaussians 4 \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 $ANNODIR/hs1/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.hs1.vcf.gz \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $ANNODIR/hs1/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.hs1.vcf.gz  \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP  -tranche 100.0 \
    -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
    -tranche 95.0 -tranche 94.0  -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0   \
    --trust-all-polymorphic ${inputchr} -O ${cohortout}/${Pro}_recalibrate_INDEL.recal \
    --tranches-file ${cohortout}/${Pro}_recalibrate_INDEL.tranches \
    --rscript-file ${cohortout}/${Pro}_recalibrate_INDEL.plots.R "

  if [ $? -ne 0 ]; then
    echo "VariantRecalibrator_INDEL failed"
    echo "${sample} VariantRecalibrator_INDEL UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1
    #comp_vqsr=0
  else
    echo "${sample} VariantRecalibrator_INDEL Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
  fi

  ${TOSCP} ${cohortout}/${Pro}_recalibrate_INDEL.recal tos://nw8k-results/vqsr/$dataset/ > /tmp/tos-gatk-p3-upload.log 2>&1
  if [ $? -ne 0 ]; then
    echo "VariantRecalibrator_INDEL Upload failed"
    echo "${sample} VariantRecalibrator_INDEL_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${sample} VariantRecalibrator_INDEL_Upload UnComplete `date +%Y/%m/%d--%H:%M`"
    #exit -1;
  else
    #echo "${CHR} SiteOnly_Upload Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} VariantRecalibrator_INDEL_Upload Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP} ${cohortout}/${Pro}_recalibrate_INDEL.tranches tos://nw8k-results/vqsr/$dataset/ > /tmp/tos-gatk-p3-upload.log 2>&1
  ${TOSCP} ${cohortout}/${Pro}_recalibrate_INDEL.plots.R tos://nw8k-results/vqsr/$dataset/ > /tmp/tos-gatk-p3-upload.log 2>&1


}
fi
wait


end_time=$(date +%s)
echo "#######-- ${Pro} VariantRecalibrator SNP End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- ${Pro} VariantRecalibrator SNP Elapsed Time:$elapsed_time sec --#######"
#${TOSCP} -r ${cohortout} tos://nw8k-results/GATK4-P3/$dataset/ > /tmp/tos-gatk-p3-upload.log 2>&1

#++++++++++++++++ Ending +++++++++++++++++
work_end_time=$(date +%s)
echo "#######-- Gatk4 Pipeline partIII snp End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $work_end_time - $work_start_time))
echo "#######-- $Pro GATK4 Cohort PartIII snp Elapsed Time:$elapsed_time sec --#######"
