Scripts for genome alignment and variant calling.
Before use, it is necessary to check the sample and path information in the script and modify them according to the actual running situation.
Function:
GATKpart1_1.sbatch:fastq processing, initial genome alignment, ExpansionHunter Denovo running
GATKpart1_2.sbatch:markduplicate
GATKpart1_3.sbatch:subsequent processing of genome Bam files, SNP/InDel analysis of individual samples
GATKpart2_1.sbatch:combining SNP/InDel single sample mutation results
GATKpart2_2.sbatch:merge the results from the previous step
GATKpart3.sbatch:variantRecalibrator for SNP/InDel
GATKpart4.sbatch:ApplyVQSR for SNP/InDel
structrue_variation.sbatch:variant calling for TRV, MEV and SV
