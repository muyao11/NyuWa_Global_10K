#!/bin/bash
### Gatk Germline SNPs+Indes Pipeline, part 2: GenomicsDBImport ~~ MakeSitesOnlyVcf 
### run separately per chromosomesregion
#SBATCH --get-user-env
### Job name
#SBATCH --job-name 5.surban
### Keep Output and Error
#SBATCH -o job-%x-%j.out
### Queue name
#SBATCH -p cpu6
### Specify the number of nodes and thread (ppn) for your job.
#SBATCH -n 1
##SBATCH --nodelist=MarkDupEBS-3
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G


### need to edit the next 1 line 

CHR=surban
#${chr}
sample=$CHR
WORKDIR=/HGDP
dataset=$dataset
batch=$dataset
threads=$(nproc)
Pro=$dataset

#SiteOnly_file_list=$SiteOnly_file_list
#vfiltered_file_list=$vfiltered_file_list
firstRegion=
addRegion=
breakpoint_list=$breakpoint_list
# 
IMGDIR="/shared/data/images"
TMPDIR="$WORKDIR/tmp"
[ -e ${TMPDIR} ] && : || mkdir -p ${TMPDIR}

# genome and annotations
ANNODIR=$WORKDIR/anno
REFERENCE=$ANNODIR/hs1/hs1.fasta

###########################################################################################
####################++++++++++++++++ Download sample files +++++++++++++++++############### 
###########################################################################################
TOSUTIL=/shared/home/chenlin/500-samples-sbatch/scripts/tosutil 
TOSCP="$TOSUTIL cp -p 2 -j 5 -conf /shared/home/chenlin/500-samples-sbatch/scripts/.tosutilconfig "

#echo Download Reference Start Time:`date +%Y/%m/%d--%H:%M`
[ -e ${ANNODIR} ] && : || mkdir -p ${ANNODIR}
#if [ ! -e /wgs/download.finished ];then
#  ${TOSCP} -r tos://bioinfo-databases/nw8k/hs1/  ${ANNODIR}  > /tmp/download.wgs.log 2>&1 && touch /wgs/download.finished
#fi

#${TOSCP} -r ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz tos://nw8k-results/vfilter/$dataset/$chr/ ${VCFDIRvfilter} >> /tmp/upload.log 2>&1
#${TOSCP} -r ${cohortout}/${Pro}_${CHR}_vfiltered.vcf.gz.tbi tos://nw8k-results/vfilter/$dataset/$chr/ ${VCFDIRvfilter} >> /tmp/upload.log 2>&1

#${TOSCP} -r ${cohortout}/${Pro}_${CHR}_SiteOnly.vcf.gz tos://nw8k-results/SiteOnly/$dataset/$chr/ ${VCFDIRsiteonly} >> /tmp/upload.log 2>&1
#${TOSCP} -r ${cohortout}/${Pro}_${CHR}_SiteOnly.vcf.gz.tbi tos://nw8k-results/SiteOnly/$dataset/$chr ${VCFDIRsiteonly} >> /tmp/upload.log 2>&1


#echo Download Reference Finished Time:`date +%Y/%m/%d--%H:%M`

C_LOGS=/shared/data/nw8k/log/$batch/$sample.complete.log
U_LOGS=/shared/data/nw8k/log/$batch/$sample.uncomplete.log

echo "#######--  Pipeline partII Start time is `date +%Y/%m/%d--%H:%M` --#######"
hostname
work_start_time=$(date +%s)

#++++++++++++++++ Merge +++++++++++++++++
cohortout=${WORKDIR}/GATKprocess/06_${Pro}_Cohort
[ -e ${cohortout} ] && : || mkdir -p ${cohortout}

echo "#######-- chromsome_${CHR} Merge Start Time:`date +%Y/%m/%d--%H:%M` --#######"
start_time=$(date +%s)
echo "${CHR}"

for t in SiteOnly vfiltered ; do

  addvcfs=""

  for region in $addRegion; do
    addvcfs="$addvcfs,${cohortout}/${Pro}_${region}_$t.vcf.gz"
  done
  addvcfs=`echo $addvcfs |sed -e 's/^,//'`

  singularity exec  -B /shared:/shared -B $WORKDIR:$WORKDIR  ${IMGDIR}/biotools_ibp_snp_v2.sif \
    bash -c " set -o pipefail
    python /shared/data/nw8k/catvcfs.py ${cohortout}/${Pro}_${firstRegion}_$t.vcf.gz ${addvcfs} ${breakpoint_list} | bgzip --threads 5 > $cohortout/${Pro}_${CHR}_$t.vcf.gz && \
    tabix $cohortout/${Pro}_${CHR}_$t.vcf.gz"

  if [ $? -ne 0 ]; then
    echo "Merge_$t failed"
    echo "${CHR} Merge_$t UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${CHR} Merge_$t UnComplete `date +%Y/%m/%d--%H:%M`"
    exit -1;
  else
    echo "${CHR} Merge_$t Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} Merge_$t Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP} -r ${cohortout}/${Pro}_${CHR}_$t.vcf.gz tos://nw8k-results/final_vcf/$dataset/ >> /tmp/upload.log 2>&1

  if [ $? -ne 0 ]; then
    echo "$t Upload failed"
    echo "${CHR} Merge_${t}_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    echo "${CHR} Merge_${t}_Upload UnComplete `date +%Y/%m/%d--%H:%M`"
    #exit -1;
  else
    #echo "${CHR} SiteOnly_Upload Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    echo "${CHR} Merge_${t}_Upload Complete `date +%Y/%m/%d--%H:%M`"
  fi

  ${TOSCP} -r ${cohortout}/${Pro}_${CHR}_$t.vcf.gz.tbi tos://nw8k-results/final_vcf/$dataset/ >> /tmp/upload.log 2>&1

done

end_time=$(date +%s)
echo "#######-- chromsome_${CHR} Merge End Time:`date +%Y/%m/%d--%H:%M` --#######"
((elapsed_time = $end_time - $start_time))
echo "#######-- chromsome_${CHR} Merge Elapsed Time:$elapsed_time sec --#######"

