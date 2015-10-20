#!/bin/bash -l

# This script is an example test of the ataqc package

#export module
module add bedtools java picard-tools preseq python_anaconda r samtools

# Directories and prefixes
WORKDIR='/mnt/lab_data/kundaje/projects/skin/data/fastqs/atac/2015-07-24_freshfrozen/merged/KeraFreshDay01minA'
OUTDIR='/srv/scratch/dskim89/tmp/ataqc'
OUTPREFIX='KeraFreshDay01minA_1.trim'
INPREFIX='KeraFreshDay01minA_1.trim'

# Annotation files
ANNOTDIR='/users/dskim89/git/ataqc/annotations'
DNASE_BED="${ANNOTDIR}/reg2map_honeybadger2_dnase_all_p2.bed.gz"
BLACKLIST_BED="${ANNOTDIR}/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz"
TSS_BED='/users/dskim89/git/pipelines/atac/tss/parsed_hg19_RefSeq.merged.ANS.bed'
REF_FASTA='/srv/scratch/dskim89/tmp/encodeHg19Male.fa'

python run_ataqc.py --workdir $WORKDIR --outdir $OUTDIR --outprefix $OUTPREFIX --ref $REF_FASTA --tss $TSS_BED --dnase $DNASE_BED --blacklist $BLACKLIST_BED --pipeline kundajelab --inprefix $INPREFIX
