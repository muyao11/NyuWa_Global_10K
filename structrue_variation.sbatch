#!/bin/bash
### Inherit all current environment variables
#SBATCH --get-user-env
### Keep Output and Error
#SBATCH -o job-%x-%j.out
### Job name
#SBATCH --job-name 3.surban 
### Queue name
#SBATCH -p C16G64
### Specify the number of nodes and thread (ppn) for your job.
#SBATCH -n 1
##SBATCH --nodelist=MarkDupEBS-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=56G


echo Pipeline Start time @ $(date +%Y/%m/%d--%H:%M)
work_start_time=$(date +%s)
hostname

sample=surban
sex=$sex
bam=$bam
bai=$bai
dataset=$dataset
threads=$(nproc)
WORKDIR=/data
REFDIR=$WORKDIR/ref
batch=
#[ -e ${REFDIR} ] && : || mkdir -p ${REFDIR}

U_LOGS=/shared/data/nw8k/log/$batch$sample.uncomplete.log
C_LOGS=/shared/data/nw8k/log/$batch$sample.complete.log

IMGDIR="/shared/data/images"

TMPDIR="$WORKDIR/tmp"
[ -e ${TMPDIR} ] && : || mkdir -p ${TMPDIR}
TOSUTIL=/shared/home/chenlin/500-samples-sbatch/scripts/tosutil
TOSCP="$TOSUTIL cp -p 2 -j 5 -conf /shared/home/zhpn/.tosutilconfig "

# genome and annotations
ANNODIR=$WORKDIR/wgs/anno/hs1/v0
REFERENCE=${ANNODIR}/multihs1/hs1.fasta

rm -r ${WORKDIR}/01_CleanData/
rm -r ${WORKDIR}/03_BwaOut/
rm -r ${WORKDIR}/04_GatkOut/
rm -r ${WORKDIR}/advntr*
rm -r ${WORKDIR}/megane

####++++++++++++++++ Download sample files +++++++++++++++++

echo Download Reference Start Time:`date +%Y/%m/%d--%H:%M`
[ -e ${ANNODIR} ] && : || mkdir -p ${ANNODIR}
echo Download Start Time:`date +%Y/%m/%d--%H:%M`
if [ ! -e $WORKDIR/wgs/download.multihs1.finished ];then
    ${TOSCP} -r tos://bioinfo-databases/nw8k/multihs1  ${ANNODIR}  > /tmp/download.wgs.log 2>&1 && touch $WORKDIR/wgs/download.multihs1.finished
fi

gatkoutdir=${WORKDIR}/04_GatkOut/${sample}
[ -e ${gatkoutdir} ] && : || mkdir -p ${gatkoutdir}
gatkout=${gatkoutdir}/${sample}

input_bam=$gatkout.final.bam

if [ ! -e $input_bam ];then
${TOSCP} tos://nw8k-results/finalbam/$dataset/$sample/$sample.final.bam $gatkout.final.bam >> /tmp/bam-download.log 2>&1

${TOSCP} tos://nw8k-results/finalbam/$dataset/$sample/$sample.final.bam.bai $gatkout.final.bam.bai >> /tmp/bam-download.log 2>&1
fi

if [ ! -e $WORKDIR/mrefdownload.finished ];then
     ${TOSCP} -r tos://bioinfo-databases/ref ${WORKDIR}  > /tmp/download.ref.log 2>&1 && touch $WORKDIR/mrefdownload.finished
fi

touch ${ANNODIR}/multihs1/hs1.fasta.fai

if [ ! -e $WORKDIR/dbrefdownload.finished ];then
     ${TOSCP} -r -flat tos://bioinfo-databases/dbrefv2  ${WORKDIR}/ref/  > /tmp/download.ref.log 2>&1 && touch $WORKDIR/dbrefdownload.finished
fi


echo Download Reference End Time:`date +%Y/%m/%d--%H:%M`

run=true
##########################################################################################
##################################hapcaller###############################################
if $run; then
{

  echo "#######-- HaplotypeCaller Start Time:`date +%Y/%m/%d--%H:%M` --#######"
  start_time=$(date +%s)
  singularity exec  -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp  ${IMGDIR}/gatk4_4.5.0.sif \
    gatk HaplotypeCaller --java-options " -Xmx40G -XX:ParallelGCThreads=5 -Djava.io.tmpdir=/tmp/" \
    -R ${REFERENCE} -ERC GVCF -I ${gatkout}.final.bam \
    -O ${gatkout}.raw.realign.snps.indels.g.vcf.gz  --tmp-dir /tmp/ \
    --native-pair-hmm-threads 5 --seconds-between-progress-updates 300
  if [ $? -ne 0 ]; then
    echo "HaplotypeCaller failed"
    echo "${sample} HaplotypeCaller UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  #exit -1;
  fi

  ${TOSCP} ${gatkout}.raw.realign.snps.indels.g.vcf.gz  tos://nw8k-results/gvcf/$dataset/ > /tmp/upload.log 2>&1

  if [ $? -ne 0 ]; then
    echo "HaplotypeCaller Upload failed"
    echo "${sample} HaplotypeCaller_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
  #exit -1;
  else
    echo "${sample} HaplotypeCaller Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
  fi

  ${TOSCP} ${gatkout}.raw.realign.snps.indels.g.vcf.gz.tbi  tos://nw8k-results/gvcf/$dataset/ > /tmp/upload.log 2>&1

  end_time=$(date +%s)
  echo "#######-- HaplotypeCaller End Time:`date +%Y/%m/%d--%H:%M` --#######"
  ((elapsed_time = $end_time - $start_time))
  echo "#######-- HaplotypeCaller Elapsed Time:$elapsed_time sec --#######"
} &
fi

##########################################################################################
##################################  adVNTR  ##############################################
comp_ad=1
if $run; then

  for i in {0..10};do
  {
    echo "######-- adVNTR $i start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)

    OUTDIR=${WORKDIR}/advntr_$i
    [ -e ${OUTDIR}/$sample ] && : || mkdir -p ${OUTDIR}/$sample

    echo "advntrExtracting Allosome Bam Start time is `date +%Y/%m/%d--%H:%M`"
    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/ \
        ${IMGDIR}/advntr-v1.4.1.sif /bin/bash -c "
            source activate advntr
            cd $OUTDIR
            /opt/miniconda3/envs/advntr/bin/advntr genotype -t 1  --alignment_file $input_bam  \
                --working_directory $OUTDIR/ -m $REFDIR/chang_$i.db \
                -of vcf -o $OUTDIR/${sample}_chang_$i.vcf  "

    if [ $? -ne 0 ]; then
        echo "advntr failed"
        echo "${sample} adVNTR_chang_$i UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_ad=0
    fi
    rm $OUTDIR/*.unmapped.fasta

    echo "Extracting Allosome Bam End time is `date +%Y/%m/%d--%H:%M`"

    #**************** Ending Check ********************

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  ${IMGDIR}/biotools_ibp_v2.sif \
      sh -c " bgzip $OUTDIR/${sample}_chang_$i.vcf && tabix $OUTDIR/${sample}_chang_$i.vcf.gz "
    if [ $? -eq 0 ]; then
    #if [ -f $OUTDIR/${sample}_chang_$i.vcf ];then
        echo "${sample} chang $i Complete `date +%Y/%m/%d--%H:%M`"
    else
        echo "${sample} chang $i UnComplete"
        echo "${sample} adVNTR_chang_bgzip_$i UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        comp_ad=0
    fi

    echo ${sample} adVNTR $i complete `date +%Y/%m/%d--%H:%M`

    end_time=$(date +%s)
    ((elapsed_time = $end_time - $start_time))
    echo "#######-- adVNTR $i Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
    echo "#######-- adVNTR $i Elapsed Time:$elapsed_time sec --#######"
  } &
  sleep 300
  done

fi

wait

if $run;then

  for i in {11..22};do 
  {
    echo "######-- adVNTR $i start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)

    OUTDIR=${WORKDIR}/advntr_$i
    [ -e ${OUTDIR}/$sample ] && : || mkdir -p ${OUTDIR}/$sample

    echo "advntrExtracting Allosome Bam Start time is `date +%Y/%m/%d--%H:%M`"
    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  \
        ${IMGDIR}/advntr-v1.4.1.sif /bin/bash -c "
            source activate advntr
            cd $OUTDIR
            /opt/miniconda3/envs/advntr/bin/advntr genotype -t 1  --alignment_file $input_bam  \
                --working_directory $OUTDIR/ -m $REFDIR/chang_$i.db \
                -of vcf -o $OUTDIR/${sample}_chang_$i.vcf  "

    if [ $? -ne 0 ]; then
        echo "advntr failed"
        echo "${sample} adVNTR_chang_$i UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_ad=0
    fi
    rm $OUTDIR/*.unmapped.fasta

    echo "Extracting Allosome Bam End time is `date +%Y/%m/%d--%H:%M`"

    #**************** Ending Check ********************
    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  ${IMGDIR}/biotools_ibp_v2.sif \
      sh -c " bgzip $OUTDIR/${sample}_chang_$i.vcf && tabix $OUTDIR/${sample}_chang_$i.vcf.gz "
    if [ $? -eq 0 ]; then
    #if [ -f $OUTDIR/${sample}_chang_$i.vcf ];then
        echo "${sample} chang $i Complete `date +%Y/%m/%d--%H:%M`"
    else
        echo "${sample} chang $i UnComplete"
        echo "${sample} adVNTR_chang_bgzip_$i UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        comp_ad=0
    fi

    echo ${sample} adVNTR complete `date +%Y/%m/%d--%H:%M`

    end_time=$(date +%s)
    ((elapsed_time = $end_time - $start_time))
    echo "#######-- adVNTR  Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
    echo "#######-- adVNTR  Elapsed Time:$elapsed_time sec --#######"

  } &
  sleep 300
  done

 if [ $sex -eq 0 ]; then

  echo "########## Skip adVNTR X ###########"
  while false ; do
  #for i in {0..1}; do
  {
    echo "######-- adVNTR X $i start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)

    OUTDIR=${WORKDIR}/advntr_X_${i}
    [ -e ${OUTDIR}/$sample ] && : || mkdir -p ${OUTDIR}/$sample

    echo "advntrExtracting Allosome Bam Start time is `date +%Y/%m/%d--%H:%M`"
    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/ \
        ${IMGDIR}/advntr-v1.4.1.sif /bin/bash -c "
            source activate advntr
            cd $OUTDIR
            /opt/miniconda3/envs/advntr/bin/advntr genotype -t 1 --haploid  \
                --alignment_file $input_bam  \
                --working_directory $OUTDIR/ -m $REFDIR/X_${i}.db -of vcf \
                -o $OUTDIR/${sample}_X_${i}.vcf "
    if [ $? -ne 0 ]; then
        echo "advntr failed"
        echo "${sample} adVNTR_X_${i} UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_ad=0
    fi

    echo "Extracting Allosome Bam End time is `date +%Y/%m/%d--%H:%M`"

    #**************** Ending Check ********************
    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  ${IMGDIR}/biotools_ibp_v2.sif \
      sh -c " bgzip $OUTDIR/${sample}_X_$i.vcf && tabix $OUTDIR/${sample}_X_$i.vcf.gz "
    if [ $? -eq 0 ]; then
    #if [ -f $OUTDIR/${sample}_X_$i.vcf ];then
        echo "${sample} X Complete `date +%Y/%m/%d--%H:%M`"
    else
        echo "${sample} X UnComplete"
        echo "${sample} adVNTR_X_$i UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_ad=0
    fi

    echo ${sample} adVNTR X complete `date +%Y/%m/%d--%H:%M`

    end_time=$(date +%s)
    ((elapsed_time = $end_time - $start_time))
    echo "#######-- adVNTR  Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
    echo "#######-- adVNTR  Elapsed Time:$elapsed_time sec --#######"

  } &
  sleep 300
  done

 else
  for i in {0..3};do
  {
    echo "######-- adVNTR XY $i start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)

    OUTDIR=${WORKDIR}/advntr_XY_$i
    [ -e ${OUTDIR}/$sample ] && : || mkdir -p ${OUTDIR}/$sample

    echo "advntrExtracting Allosome Bam Start time is `date +%Y/%m/%d--%H:%M`"
    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/ \
        ${IMGDIR}/advntr-v1.4.1.sif /bin/bash -c "
            source activate advntr
            cd $OUTDIR
            /opt/miniconda3/envs/advntr/bin/advntr genotype -t 1 --haploid  \
                --alignment_file $input_bam  \
                --working_directory $OUTDIR/ -m $REFDIR/XY_$i.db -of vcf \
                -o $OUTDIR/${sample}_XY_$i.vcf "
    if [ $? -ne 0 ]; then
        echo "advntr failed"
        echo "${sample} advntr_XY_$i UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_ad=0
    fi

    echo "Extracting Allosome Bam End time is `date +%Y/%m/%d--%H:%M`"

    #**************** Ending Check ********************
    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  ${IMGDIR}/biotools_ibp_v2.sif \
      sh -c " bgzip $OUTDIR/${sample}_XY_$i.vcf && tabix $OUTDIR/${sample}_XY_$i.vcf.gz "
    if [ $? -eq 0 ]; then
    #if [ -f $OUTDIR/${sample}_XY_$i.vcf ];then
        echo "${sample} XY$i Complete `date +%Y/%m/%d--%H:%M`"
    else
        echo   "${sample} XY$i UnComplete "
        echo "${sample} adVNTR_XY_${i} UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_ad=0
    fi

  } &
  sleep 300
  done
 
 fi
fi

wait

if $run; then

for i in {0..22};do
  ${TOSCP} ${WORKDIR}/advntr_$i/${sample}_chang_$i.vcf.gz tos://nw8k-results/advntr/$dataset/$sample/ >> /tmp/${sample}.advntr.upload.log 2>&1
  if [ $? -ne 0 ]; then
    echo "${sample} adVNTR_chang_${i}_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    comp_ad=0
  fi
  ${TOSCP} ${WORKDIR}/advntr_$i/${sample}_chang_$i.vcf.gz.tbi tos://nw8k-results/advntr/$dataset/$sample/ >> /tmp/${sample}.advntr.upload.log 2>&1
done
if [ $sex -eq 0 ]; then
  while false ; do
  #for i in 0 1;do
    ${TOSCP} -r ${WORKDIR}/advntr_X_$i/${sample}_X_$i.vcf.gz tos://nw8k-results/advntr/$dataset/$sample/ >> /tmp/${sample}.advntr.upload.log 2>&1
    if [ $? -ne 0 ]; then
      echo "${sample} adVNTR_X_${i}_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
      comp_ad=0
    fi
    ${TOSCP} -r ${WORKDIR}/advntr_X_$i/${sample}_X_$i.vcf.gz.tbi tos://nw8k-results/advntr/$dataset/$sample/ >> /tmp/${sample}.advntr.upload.log 2>&1
  done
else
  for i in  0 1 2 3; do
    ${TOSCP} -r ${WORKDIR}/advntr_XY_$i/${sample}_XY_$i.vcf.gz tos://nw8k-results/advntr/$dataset/$sample/ >> /tmp/${sample}.advntr.upload.log 2>&1
    if [ $? -ne 0 ]; then
      echo "${sample} adVNTR_XY_${i}_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
      comp_ad=0
    fi
    ${TOSCP} -r ${WORKDIR}/advntr_XY_$i/${sample}_XY_$i.vcf.gz.tbi tos://nw8k-results/advntr/$dataset/$sample/ >> /tmp/${sample}.advntr.upload.log 2>&1
  done
fi
wait

if [ $comp_ad -eq 1 ]; then
  echo "${sample} adVNTR Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
fi

rm -rf ${WORKDIR}/advntr*

fi

#####################advntr end######################
#######################SV############################

##########################################################################################
####################################  manta pipeline  ####################################
##########################################################################################
comp_manta=1
if $run; then
{
    echo "#######-- Manta start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)

    ##################  tools
    path2manta=/opt/miniconda3/envs/python27/share/manta-1.6.0-2/bin/configManta.py
    #path2convertInversion=/opt/miniconda3/envs/python27/bin/convertInversion.py
    path2convertInversion=/shared/data/Parastor300s_G30S/zhangjj/software/manta.convertInversion.change.py
    path2samtools=/opt/miniconda3/bin/samtools
    path2bgzip=/opt/miniconda3/envs/python27/bin/bgzip
    path2tabix=/opt/miniconda3/envs/python27/bin/tabix

    OUTDIR=${WORKDIR}/manta
    [ -e ${OUTDIR}/$sample ] && : || mkdir -p ${OUTDIR}/$sample

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/ \
        ${IMGDIR}/biotools_ibp_v2.sif /bin/bash -c "
            mamba init
            . ~/.bashrc
            mamba activate python27
            echo 'start manta'
            $path2manta --bam $input_bam --referenceFasta $REFERENCE --runDir $OUTDIR/$sample;
            $OUTDIR/$sample/runWorkflow.py -j 8
        "

    if [ $? -ne 0 ]; then
        echo "manta workflow failed"
        echo "${sample} manta UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_manta=0
    fi

    echo "#######-- manta workflow End time:`date +%Y/%m/%d--%H:%M`  --#######"

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/ -B /shared/data:/shared/data \
        ${IMGDIR}/biotools_ibp_v2.sif /bin/bash -c "
            mamba init
            . ~/.bashrc
            mamba activate python27
            $path2convertInversion $path2samtools $REFERENCE $OUTDIR/$sample/results/variants/diploidSV.vcf.gz | \
                $path2bgzip -c >$OUTDIR/$sample/results/variants/diploidSV.convertInversion.vcf.gz
            $path2tabix -p vcf $OUTDIR/$sample/results/variants/diploidSV.convertInversion.vcf.gz
        "

    echo "#######-- manta convert inversion End time:`date +%Y/%m/%d--%H:%M`  --#######"

    #**************** Ending Check ********************
    if [ -e ${OUTDIR}/$sample/results/variants/diploidSV.convertInversion.vcf.gz.tbi ];then
        echo "${sample} manta convertInversion Complete `date +%Y/%m/%d--%H:%M`"
    else
        echo "${sample} manta convertInversion UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        comp_manta=0
    fi

    ${TOSCP} -r $OUTDIR/$sample/results/variants tos://nw8k-results/manta/$dataset/$sample/ >> /tmp/${sample}.manta.upload.log 2>&1
    if [ $? -ne 0 ]; then
      echo "${sample} manta_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
      comp_manta=0
    fi

    if [ $comp_manta -eq 1 ]; then
      echo ${sample} manta complete `date +%Y/%m/%d--%H:%M`
      echo "${sample} Manta Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    else
      echo "${sample} Manta UnComplete `date +%Y/%m/%d--%H:%M`"
    fi

    end_time=$(date +%s)
    ((elapsed_time = $end_time - $start_time))
    echo "#######-- manta All Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
    echo "#######-- manta All Elapsed Time:$elapsed_time sec --#######"

} &
fi

##########################################################################################
#################################  CNVnator pipeline  ####################################
##########################################################################################
comp_cnv=1
if $run;then
{
    echo "#######-- cnvnator start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)

    bin_size=100

    path2cnvnator=/usr/src/CNVnator/cnvnator
    path2cnvnator2VCF=/usr/src/CNVnator/cnvnator2VCF.pl

    OUTDIR=${WORKDIR}/cnvnator/$sample
    [ -e ${OUTDIR} ] && : || mkdir -p ${OUTDIR}

    #**************** cnvnator pipeline ********************

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  \
        ${IMGDIR}/cnvnator.sif /bin/bash -c " \
        cd $OUTDIR
        export OMP_NUM_THREADS=6
        echo '## step1: extract mapping information'
        $path2cnvnator -root $OUTDIR/$sample.root -tree $input_bam
        echo '## step2: generate histogram'
        $path2cnvnator -root $OUTDIR/$sample.root -genome hs1 -his ${bin_size} -d $REFERENCE
        echo '## step3: generate the statistic'
        $path2cnvnator -root $OUTDIR/$sample.root -stat ${bin_size}
        echo '## step4: RD information partition'
        $path2cnvnator -root $OUTDIR/$sample.root -partition ${bin_size}
        echo '## step5: call cnv'
        $path2cnvnator -root $OUTDIR/$sample.root -call ${bin_size} >$OUTDIR/${sample}.cnv
        echo '## step6: cnv2vcf'
        $path2cnvnator2VCF $OUTDIR/${sample}.cnv > $OUTDIR/${sample}.vcf
    "

    if [ $? -ne 0 ]; then
        echo "cnvnator failed"
        echo "${sample} cnvnator UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_cnv=0
    fi
    end_time=$(date +%s)
    echo "#######-- cnvnator End time:`date +%Y/%m/%d--%H:%M`  --#######"

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  ${IMGDIR}/biotools_ibp_v2.sif \
      sh -c "bgzip -c $OUTDIR/${sample}.vcf >$OUTDIR/${sample}.vcf.gz && tabix $OUTDIR/${sample}.vcf.gz "

    if [ $? -ne 0 ]; then
        echo "cnvnator bgzip failed"
        #echo "${sample} cnvnator_bgzip UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        #comp_cnv=0
    fi

    #$singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  ${IMGDIR}/biotools_ibp_v2.sif \
    #  sh -c "bgzip $OUTDIR/${sample}.cnv "

    ${TOSCP} $OUTDIR/${sample}.vcf.gz tos://nw8k-results/cnv/$dataset/$sample/  >> /tmp/${sample}.cnvnator.upload.log 2>&1
    if [ $? -ne 0 ]; then
        echo "cnvnator Upload failed"
        echo "${sample} cnvnator_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_cnv=0
    fi
    ${TOSCP} $OUTDIR/${sample}.vcf.gz.tbi tos://nw8k-results/cnv/$dataset/$sample/  >> /tmp/${sample}.cnvnator.upload.log 2>&1
    ${TOSCP} $OUTDIR/${sample}.cnv tos://nw8k-results/cnv/$dataset/$sample/  >> /tmp/${sample}.cnvnator.upload.log 2>&1

    #**************** Ending Check ********************

    if [ $comp_cnv -eq 1 ]; then
      echo ${sample} cnvnator complete `date +%Y/%m/%d--%H:%M`
      echo "${sample} CNVnator Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    else
        echo "${sample} CNVnator UnComplete"
    fi
    end_time=$(date +%s)
    ((elapsed_time = $end_time - $start_time))
    echo "#######-- cnvnator  Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
    echo "#######-- cnvnator  Elapsed Time:$elapsed_time sec --#######"
} &
fi

wait

##########################################################################################
#################################  whamg  pipeline    ####################################
##########################################################################################
comp_wh=1
if $run;then
{
    echo "######-- whamg start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)

    OUTDIR=${WORKDIR}/whamg/$sample
    [ -e ${OUTDIR} ] && : || mkdir -p ${OUTDIR}

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/ \
        ${IMGDIR}/wham_latest.sif /bin/bash -c "
            export EXCLUDE=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM
            whamg -x 8 -c \$EXCLUDE -a $REFERENCE -f $input_bam > $OUTDIR/${sample}.vcf 2> $OUTDIR/${sample}.err
    "
    if [ $? -ne 0 ]; then
        echo "whamg failed"
        echo "${sample} whamg UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_wh=0
    fi

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  ${IMGDIR}/biotools_ibp_v2.sif \
      sh -c "bgzip -c $OUTDIR/${sample}.vcf >$OUTDIR/${sample}.vcf.gz && tabix $OUTDIR/${sample}.vcf.gz "
    if [ $? -ne 0 ]; then
        echo "whamg bgzip failed"
        echo "${sample} whamg_bgzip UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        comp_wh=0
    fi

    ${TOSCP} $OUTDIR/${sample}.vcf.gz tos://nw8k-results/whamg/$dataset/$sample/ >> /tmp/${sample}.whamg.upload.log 2>&1
    if [ $? -ne 0 ]; then
        echo "whamg Upload failed"
        echo "${sample} whamg_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        comp_wh=0
    fi

    ${TOSCP} $OUTDIR/${sample}.vcf.gz.tbi tos://nw8k-results/whamg/$dataset/$sample/ >> /tmp/${sample}.whamg.upload.log 2>&1

    end_time=$(date +%s)
    echo "#######-- whamg End time:`date +%Y/%m/%d--%H:%M`  --#######"

    #**************** Ending Check ********************

    if [ $comp_wh -eq 1 ]; then
      echo "${sample} whamg complete `date +%Y/%m/%d--%H:%M`"
      echo "${sample} whamg Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    else
      echo "${sample} UnComplete `date +%Y/%m/%d--%H:%M`"
    fi
    end_time=$(date +%s)
    ((elapsed_time = $end_time - $start_time))
    echo "#######-- whamg  Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
    echo "#######-- whamg  Elapsed Time:$elapsed_time sec --#######"
} &
fi

##########################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 等待后台进程结束， 勿删除该语句，否则会导致流程步骤结果不完整
echo "wait for jobs to finish"
wait

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##########################################################
#MOPLINE
comp_mop=1
if $run; then
  echo "######-- MOPline start time is `date +%Y/%m/%d--%H:%M`"
  start_time=$(date +%s)
  OUTDIR=${WORKDIR}/MOPline
  [ -e ${OUTDIR} ] && : || mkdir -p ${OUTDIR}
  [ -e ${OUTDIR}/$sample ] && : || mkdir -p ${OUTDIR}/$sample

  echo "${sample}" > $OUTDIR/sample.lst
  echo "${WORKDIR}/04_GatkOut/${sample} ${sample}.final.bam" > $OUTDIR/bam.lst

  mkdir -p ${OUTDIR}/$sample/MELT/
  touch ${OUTDIR}/$sample/MELT/MELT.${sample}.vcf

  mkdir -p ${OUTDIR}/$sample/MATCHCLIP/
  touch ${OUTDIR}/$sample/MATCHCLIP/MATCHCLIP.${sample}.vcf

  mkdir -p ${OUTDIR}/$sample/GRIDSS/
  touch ${OUTDIR}/$sample/GRIDSS/GRIDSS.${sample}.vcf

  mkdir -p ${OUTDIR}/$sample/inGAP/
  touch ${OUTDIR}/$sample/inGAP/inGAP.${sample}.vcf
  mkdir -p ${OUTDIR}/$sample/Manta/
  mkdir -p ${OUTDIR}/$sample/CNVnator/
  mkdir -p ${OUTDIR}/$sample/Wham/
  mopdir=/opt/local/tools
  singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/ \
    /shared/home/zhangjj/dockerfiles/mopline.sif \
    sh -c "$mopdir/scripts/run_SVcallers/convert_Manta_vcf.pl ${WORKDIR}/manta/${sample}/results/variants/diploidSV.vcf.gz > ${OUTDIR}/$sample/Manta/Manta.${sample}.vcf;\
    $mopdir/scripts/run_SVcallers/convert_CNVnator_vcf.pl ${WORKDIR}/cnvnator/${sample}/${sample}.cnv $mopdir/MOPline/Data/gap.bed >${OUTDIR}/$sample/CNVnator/CNVnator.${sample}.vcf;\
    $mopdir/scripts/run_SVcallers/convert_Wham_vcf.pl ${WORKDIR}/whamg/${sample}/${sample}.vcf  >${OUTDIR}/$sample/Wham/Wham.${sample}.vcf;\
    cd ${OUTDIR}/$sample;\
    echo $threads;\
    $mopdir/mopline create_cov -b $OUTDIR/bam.lst -r $REFERENCE -rl 150 -n $threads;\
    $mopdir/mopline merge_7tools --sample ${sample} -rl 150;\
    cd ${OUTDIR}
    $mopdir/mopline add_cov --sample_list $OUTDIR/sample.lst --toolset 7tools --build T2T --vcf_dname Merge_7tools -n 4;\
  "
    if [ $? -ne 0 ]; then
        echo "mopline failed"
        echo "${sample} mopline UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_mop=0
    fi

  ${TOSCP} -r $OUTDIR/$sample tos://nw8k-results/MOPline/$dataset/ >> /tmp/${sample}.mopline.upload.log 2>&1
  if [ $? -ne 0 ]; then
    echo "mopline Upload failed"
    echo "${sample} mopline_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    comp_mop=0
  fi

  if [ $comp_wh -eq 1 ]; then
    echo ${sample} mopline complete `date +%Y/%m/%d--%H:%M`
    echo "${sample} mopline Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    rm -rf $WORKDIR/MOPline
    rm -rf $WORKDIR/manta
    rm -rf $WORKDIR/cnvnator
    rm -rf $WORKDIR/Whamg
    #rm -rf $OUTDIR/bam.lst
    #rm -rf $OUTDIR/sample.lst
  fi
end_time=$(date +%s)
((elapsed_time = $end_time - $start_time))
echo "#######-- mopline  Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
echo "#######-- mopline  Elapsed Time:$elapsed_time sec --#######"

fi

##########################################################################################
#################################  gangSTR  pipeline   ###################################
##########################################################################################
comp_gan=1
if $run;then
{
    echo "######-- gangSTR start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)

    OUTDIR=${WORKDIR}/gangSTR/$sample
    [ -e ${OUTDIR} ] && : || mkdir -p ${OUTDIR}

    ################# run gangstr autosome
    echo gangstr autosome start time is `date +%Y/%m/%d--%H:%M`

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  \
        ${IMGDIR}/biotools_ibp_v2.sif /bin/bash -c "
            GangSTR --ref $REFERENCE --regions $REFDIR/hs1_gangstr.bed  \
              --bam $input_bam  --out $OUTDIR/$sample --include-ggl

            echo ${sample} GangSTR complete `date +%Y/%m/%d--%H:%M`

            echo sort zip index start time is `date +%Y/%m/%d--%H:%M`
            cd $OUTDIR
            bcftools sort $sample.vcf > $sample-sort.vcf
            bgzip $sample-sort.vcf
            tabix -p vcf -f $sample-sort.vcf.gz
            echo sort zip index finish time is `date +%Y/%m/%d--%H:%M`

            echo gangstr autosome finish time is `date +%Y/%m/%d--%H:%M`
            ################# run gangstr allosome
            echo gangstr allosome start time is `date +%Y/%m/%d--%H:%M`
            # chrX-haploid
            if [ $sex -eq 1 ];then
            GangSTR  --bam  $input_bam  --ref $REFERENCE  \
                  --regions $REFDIR/hs1_gangstr.bed  \
                  --bam-samps $sample  --samp-sex  M --chrom chrX   \
                  --out $OUTDIR/$sample.chrX  --ploidy 1

            bcftools sort  $OUTDIR/$sample.chrX.vcf  >  $OUTDIR/$sample.chrX.sort.vcf
            bgzip   $OUTDIR/$sample.chrX.sort.vcf
            tabix  -p  vcf -f  $OUTDIR/$sample.chrX.sort.vcf.gz
            #mv $OUTDIR/$sample.chrX.vcf  $OUTDIR/trash
            # chrY-haploid
            GangSTR  --bam  $input_bam  --ref $REFERENCE  \
                --regions $REFDIR/hs1_gangstr.bed   \
                --bam-samps $sample  --samp-sex  M --chrom chrY   \
                --out $OUTDIR/$sample.chrY  --ploidy 1

            bcftools sort  $OUTDIR/$sample.chrY.vcf  >  $OUTDIR/$sample.chrY.sort.vcf
            bgzip   $OUTDIR/$sample.chrY.sort.vcf
            tabix  -p  vcf -f  $OUTDIR/$sample.chrY.sort.vcf.gz
            #mv  $OUTDIR/$sample.chrY.vcf  $OUTDIR/trash

            # chrX-diploid
            else
            GangSTR  --bam  $input_bam  --ref $REFERENCE  \
                --regions $REFDIR/hs1_gangstr.bed  \
                --bam-samps $sample  --samp-sex  F --chrom chrX   \
                --out $OUTDIR/$sample.chrXdiploid  --ploidy 2

            bcftools sort  $OUTDIR/$sample.chrXdiploid.vcf  >  $OUTDIR/$sample.chrXdiploid.sort.vcf
            bgzip   $OUTDIR/$sample.chrXdiploid.sort.vcf
            tabix  -p  vcf  -f $OUTDIR/$sample.chrXdiploid.sort.vcf.gz
            #mv  $OUTDIR/$sample.chrXdiploid.vcf  $OUTDIR/trash
            fi
        "

    if [ ! -e $OUTDIR/$sample-sort.vcf.gz.tbi ]; then
       echo "GangSTR failed"
       echo "${sample} GangSTR UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
       #exit -1;
       comp_gan=0
    fi

    rm $OUTDIR/*.vcf

    ${TOSCP} -r $OUTDIR tos://nw8k-results/gangstr/$dataset/ >> /tmp/${sample}.gangstr.upload.log 2>&1
    if [ $? -ne 0 ]; then
       echo "GangSTR Upload failed"
       echo "${sample} GangSTR_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
       comp_gan=0
    fi

    if [ $comp_gan -eq 1 ]; then
      echo "GangSTR Complete"
      echo "${sample} GangSTR Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    fi

    end_time=$(date +%s)
    ((elapsed_time = $end_time - $start_time))
    echo "#######-- gangstr allosome  Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
    echo "#######-- gangstr allosome  Elapsed Time:$elapsed_time sec --#######"
rm -rf ${WORKDIR}/gangSTR
} &
fi

##########################################################################################
#################################  MEGAnE  pipeline   ####################################
##########################################################################################
comp_meg=1
if $run;then
{
if [ ! -e $WORKDIR/wgs/meirefdownload.finished ];then
  ${TOSCP} -r tos://bioinfo-databases/meiref/    $WORKDIR/wgs    > /tmp/download.wgs.log 2>&1 && touch $WORKDIR/wgs/meirefdownload.finished
fi

    echo "######-- megane start time is `date +%Y/%m/%d--%H:%M`"
    start_time=$(date +%s)


    megane_mk=$WORKDIR/wgs/meiref/megane_kmer_set
    megane_lib=$WORKDIR/wgs/meiref/megane_repeat_lib
    path2megane=/usr/local/bin/call_genotype_38

    OUTDIR=${WORKDIR}/megane/$sample
    [ -e ${OUTDIR} ] && : || mkdir -p ${OUTDIR}

    singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  \
        ${IMGDIR}/MEGAnE.sif /bin/bash -c "
            $path2megane -i $input_bam -sample_name $sample -outdir $OUTDIR \
                -fa $REFERENCE -mk $megane_mk/hs1.fasta.mk \
                -fadb $megane_lib/T2T_blastdb -rep $megane_lib/Dfam_custom.rep \
                -repout $megane_lib/hs1.fasta.out \
                -repremove $megane_lib/non_ME_rep.txt -pA_ME $megane_lib/ME_with_pA.txt \
                -mainchr $megane_lib/main_chrs.txt -p 2
        "

    if [ $? -ne 0 ]; then
        echo "megane failed"
        echo "${sample} megane UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
        #exit -1
        comp_meg=0
    fi
    end_time=$(date +%s)
    echo "#######-- megane End time:`date +%Y/%m/%d--%H:%M`  --#######"

    ${TOSCP} -r $OUTDIR tos://nw8k-results/megane/$dataset/ >> /tmp/${sample}.megane.upload.log 2>&1
    if [ $? -ne 0 ]; then
      echo "megane Upload failed"
      echo "${sample} megane_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
      comp_meg=0
    fi

    #**************** Ending Check ********************
    cd $OUTDIR/
    if [ -s MEI_final_percentile_genotyped.vcf ];then
        echo "${sample} MEGAnE Complete `date +%Y/%m/%d--%H:%M`"
    else
        echo "${sample} MEGAnE UnComplete"
        comp_meg=0
    fi

    if [ $comp_meg -eq 1 ]; then
      echo ${sample} megane complete `date +%Y/%m/%d--%H:%M`
      echo "${sample} MEGAnE Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
      #rm -rf ${OUTDIR}
    fi

    rm -rf ${OUTDIR}
    end_time=$(date +%s)
    ((elapsed_time = $end_time - $start_time))
    echo "#######-- MEGAnE  Finish Time:`date +%Y/%m/%d--%H:%M` --#######"
    echo "#######-- MEGAnE  Elapsed Time:$elapsed_time sec --#######"
#rm -rf ${WORKDIR}/megane
} &
fi


#################################  HapCUT2  pipeline   ####################################
comp_hap=1
if $run;then
{

  OUTDIR=${WORKDIR}/hapcut/$sample
  [ -e ${OUTDIR} ] && : || mkdir -p ${OUTDIR}
  hcout=${OUTDIR}/$sample
  echo "#######-- HapCUT2 Start Time:`date +%Y/%m/%d--%H:%M` --#######"
  start_time=$(date +%s)

  if [ ! -e ${gatkout}.raw.realign.snps.indels.g.vcf.gz ]; then

    ${TOSCP} tos://nw8k-results/gvcf/$dataset/$sample.raw.realign.snps.indels.g.vcf.gz ${gatkout}.raw.realign.snps.indels.g.vcf.gz  > /tmp/upload.log 2>&1
    if [ $? -ne 0 ]; then
      echo "HapCUT2_gvcf Download failed"
      echo "${sample} HapCUT2_gvcf_Download UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
      #exit -1;
      comp_hap=0
    fi

    ${TOSCP} tos://nw8k-results/gvcf/$dataset/$sample.raw.realign.snps.indels.g.vcf.gz.tbi ${gatkout}.raw.realign.snps.indels.g.vcf.gz.tbi  > /tmp/upload.log 2>&1

  fi

  regions=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX
  singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp $IMGDIR/bcftools.sif \
    sh -c "bcftools convert -f $REFERENCE --gvcf2vcf  -e 'ALT=\"<NON_REF>\"' -r $regions ${gatkout}.raw.realign.snps.indels.g.vcf.gz | \
    bcftools view -e 'QUAL=\".\" || QUAL<30' |sed  's/,<NON_REF>//' |bcftools annotate -x QUAL,INFO,FORMAT | \
    grep -v '\*' |grep -v '\s.*,.*,' > $hcout.final.vcf "

  singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp --env LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/HapCUT2/htslib $IMGDIR/hapcut2.sif \
    extractHAIRS --bam $gatkout.final.bam --VCF $hcout.final.vcf --out $hcout.fragment --maxfragments 2000000 --triallelic 1 --indels 1 --realign_variants 1 --ref $REFERENCE

  if [ $? -ne 0 ]; then
    echo "extractHAIRS failed"
    echo "${sample} extractHAIRS UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_hap=0
  fi

  singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp $IMGDIR/hapcut2.sif \
    HAPCUT2 --fragments $hcout.fragment --VCF $hcout.final.vcf  --output $hcout.hapcut --outvcf 1

  if [ $? -ne 0 ]; then
    echo "HapCUT2 failed"
    echo "${sample} HapCUT2 UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_hap=0
  fi

  singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp $IMGDIR/gatk4_4.5.0.sif \
    sh -c "bgzip $hcout.hapcut.phased.VCF
    tabix $hcout.hapcut.phased.VCF.gz"

  ${TOSCP} ${hcout}.hapcut.phased.VCF.gz tos://nw8k-results/hapcut/$dataset/$sample/ >> /tmp/upload.log 2>&1
  if [ $? -ne 0 ]; then
    echo "HapCUT2 Upload failed"
    echo "${sample} HapCUT2_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_hap=0
  fi

  ${TOSCP} ${hcout}.hapcut.phased.VCF.gz.tbi tos://nw8k-results/hapcut/$dataset/$sample/ >> /tmp/upload.log 2>&1

  if [ $comp_hap -eq 1 ]; then
    echo ${sample} HapCUT2 complete `date +%Y/%m/%d--%H:%M`
    echo "${sample} HapCUT2 Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    rm -r $OUTDIR
  fi

  end_time=$(date +%s)
  echo "#######-- HapCUT2 End time:`date +%Y/%m/%d--%H:%M`  --#######"
  ((elapsed_time = $end_time - $start_time))
  echo "#######-- HapCUT2 Elapsed Time:$elapsed_time sec --#######"

}
fi

#################################  telomerecat pipeline   ####################################
comp_tel=1
if $run;then
{
  OUTDIR=${WORKDIR}/telomerecat/$sample
  [ -e ${OUTDIR} ] && : || mkdir -p ${OUTDIR}
  tmout=${OUTDIR}/$sample

  echo "#######--telomerelength Start time:`date +%Y/%m/%d--%H:%M` --#######"

  singularity exec  -B $WORKDIR:$WORKDIR --env LC_ALL=C.UTF-8 --env LANG=C.UTF-8 \
    $IMGDIR/telomerecat_d_4.0.2.sif \
    telomerecat bam2length -v2 $gatkout.final.bam --output $tmout.csv

  if [ $? -ne 0 ]; then
    echo "telomerecat failed"
    echo "${sample} telomerecat UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_tel=0
  fi

  ${TOSCP}  -r ${OUTDIR} tos://nw8k-results/telomerecat/$dataset/ >> /tmp/upload.log 2>&1
  if [ $? -ne 0 ]; then
    echo "telomerecat Upload failed"
    echo "${sample} telomerecat_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_tel=0
  fi

 if [ $comp_tel -eq 1 ]; then
    echo ${sample} telomerecat complete `date +%Y/%m/%d--%H:%M`
    echo "${sample} telomerecat Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
  fi

}
fi

####+++++++++++++++ strling extract job ++++++++++++++++++
comp_stl=1
if $run; then
{
  echo "#######--strling extract Start time:`date +%Y/%m/%d--%H:%M` --#######"
  mkdir -p ${WORKDIR}/strling/$sample
  singularity exec  -B $WORKDIR:$WORKDIR  $IMGDIR/strling.sif  \
    strling extract -f $REFERENCE -g $REFDIR/hs1.str -v $input_bam ${WORKDIR}/strling/$sample/$sample.str.bin

  if [ $? -ne 0 ]; then
    echo "strling failed"
    echo "${sample} strling UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_stl=0
  fi

  echo "#######--strling extract Finish time:`date +%Y/%m/%d--%H:%M` --#######"

####++++++++++++++++ upload output files +++++++++++++++++

  echo "#######--upload Start time:`date +%Y/%m/%d--%H:%M` --#######"

  ${TOSCP} -r ${WORKDIR}/strling/$sample tos://nw8k-results/strling/$dataset/ >> /tmp/upload.log 2>&1
  if [ $? -ne 0 ]; then
    echo "strling Upload failed"
    echo "${sample} strling_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_stl=0
  fi

  if [ $comp_stl -eq 1 ]; then
    echo ${sample} strling complete `date +%Y/%m/%d--%H:%M`
    echo "${sample} strling Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    rm -rf ${WORKDIR}/strling/$sample
  fi


  echo "#######--upload Finish time:`date +%Y/%m/%d--%H:%M` --#######"
  #rm -rf ${WORKDIR}/strling/$sample
} &
fi

wait

######################################################
comp_eph=1
if false; then
#{
  echo "#######--ExpansionHunter  Start Time:`date +%Y/%m/%d--%H:%M` --#######"
  mkdir -p $WORKDIR/EH/$sample
  singularity exec  -B $WORKDIR:$WORKDIR  $IMGDIR/ExpansionHunter.sif \
    ExpansionHunter --reads $input_bam --reference $REFERENCE  --variant-catalog ${REFDIR}/filtered_hs1.json  --threads 15   --output-prefix $WORKDIR/EH/$sample/$sample

  if [ $? -ne 0 ]; then
    echo "ExpansionHunter failed"
    echo "${sample} ExpansionHunter UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_eph=0
  fi

  echo "#######--ExpansionHunter Finish Time:`date +%Y/%m/%d--%H:%M` --#######"

####################++++++++++++++++ Download output files +++++++++++++++++###############

  echo "#######--ExpansionHunter Start Time:`date +%Y/%m/%d--%H:%M` --#######"

  singularity exec -B $WORKDIR:$WORKDIR -B ${TMPDIR}:/tmp/  ${IMGDIR}/biotools_ibp_v2.sif \
    sh -c " bgzip $WORKDIR/EH/$sample/${sample}.vcf && tabix $WORKDIR/EH/$sample/${sample}.vcf.gz "
  if [ $? -ne 0 ]; then
    echo "${sample} ExpansionHunter bgzip UnComplete "
    echo "${sample} ExpansionHunter_bgzip UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    comp_eph=0
  fi
  gzip $WORKDIR/EH/$sample/$sample.json
  rm $WORKDIR/EH/$sample/${sample}_realigned.bam

  ${TOSCP} -r ${WORKDIR}/EH/$sample tos://nw8k-results/EH/$dataset/ >> /tmp/upload.log 2>&1
  if [ $? -ne 0 ]; then
    echo "ExpansionHunter Upload failed"
    echo "${sample} ExpansionHunter_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    #exit -1;
    comp_eph=0
  fi
  if [ $comp_eph -eq 1 ]; then
    echo ${sample} ExpansionHunter complete `date +%Y/%m/%d--%H:%M`
    echo "${sample} ExpansionHunter Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    rm -rf ${WORKDIR}/EH/$sample
  fi

  echo "#######--ExpansionHunter Finish Time:`date +%Y/%m/%d--%H:%M` --#######"

#} &
fi

comp_dbtk=1
if $run; then
#{
#rm -rf ${WORKDIR}/EH/$sample/

  echo danbingtk 1.3.2 Pipeline Start time is `date +%Y/%m/%d--%H:%M`

  echo "#######-- danbingtk align Start Time:`date +%Y/%m/%d--%H:%M`  --#######"
  danbingout=$WORKDIR/danbing/${sample}
  mkdir -p $danbingout
  input_fasta=$danbingout/${sample}.fasta
  if [ ! -s "$input_fasta" ]; then
    singularity exec -B $WORKDIR:$WORKDIR -B $WORKDIR/tmp:/tmp  ${IMGDIR}/biotools_ibp_snp_v2.sif \
      sh -c "samtools fasta -n ${input_bam} > ${input_fasta}"
    if [ $? -ne 0 ]; then
      echo "danbingtk fasta failed"
      echo "${sample} danbingtk_Fasta UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
      exit -1;
      comp_dbtk=0
    fi
  fi
  rm ${input_bam}*

  if [ ! -f "${sample}.tr.kmers" ]; then
    singularity exec -B $WORKDIR:$WORKDIR -B $WORKDIR/tmp:/tmp  ${IMGDIR}/danbing-tk_v1.3.2.sif \
      sh -c "danbing-tk -gc 85 -ae -kf 4 1 -cth 45 -on $danbingout/${sample} -k 21 -qs $REFDIR/pan -fa $input_fasta -p 8|gzip > $danbingout/${sample}.aln.gz"
      #sh -c "danbing-tk -gc 85 -ae -kf 4 1 -cth 45 -o $WORKDIR/danbing/$sample -k 21 -qs /data/ref/pan -fa $WORKDIR/${sample}.fasta -p 8"
    if [ $? -ne 0 ]; then
      echo "danbingtk failed"
      echo "${sample} danbingtk UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
      exit -1;
      comp_dbtk=0
    fi
  fi
  rm $input_fasta
  gzip $danbingout/$sample.tr.kmers

  ${TOSCP}  $danbingout/$sample.tr.kmers.gz tos://nw8k-results/danbing-tk/$dataset/$sample/ >> /tmp/download.log 2>&1
  if [ $? -ne 0 ]; then
    echo "danbingtk Upload failed"
    echo "${sample} danbingtk_Upload UnComplete `date +%Y/%m/%d--%H:%M`" >>$U_LOGS
    exit -1;
    comp_dbtk=0
  fi

  if [ $comp_dbtk -eq 1 ]; then
    echo ${sample} danbingtk complete `date +%Y/%m/%d--%H:%M`
    echo "${sample} danbingtk Complete `date +%Y/%m/%d--%H:%M`" >>$C_LOGS
    #rm -f $input_fasta
    rm -rf $danbingout
  fi


    #rm -f input_fasta=$WORKDIR/${sample}.fasta
    #rm -rf $danbingout
  echo "#######-- danbingtk align End Time:`date +%Y/%m/%d--%H:%M`  --#######"

#} &
fi

rm ${input_bam}*




