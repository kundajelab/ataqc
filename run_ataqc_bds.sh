#!/bin/bash -l

# Example script with all arguments defined to run in BDS

# Don't forget to run `export -f module` first
module add bedtools java picard-tools preseq python_anaconda r samtools ucsc_tools

# Be sure that environmental variables are set. For example, $PICARD should
# point to the directory with the picard.jar file.

# Directories and prefixes
WORKDIR=$1
OUTDIR=$2
OUTPREFIX=$3
INPREFIX=$4
GENOME='hg19' # This is the only genome that currently works

# Annotation files
ANNOTDIR="/mnt/lab_data/kundaje/users/dskim89/ataqc/annotations"
DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
BLACKLIST_BED="${ANNOTDIR}/${GENOME}/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz"
TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"

python ~/git/ataqc/run_ataqc.py \
    --workdir $WORKDIR \
    --outdir $OUTDIR \
    --outprefix $OUTPREFIX \
    --genome $GENOME \
    --ref $REF_FASTA \
    --tss $TSS_BED \
    --dnase $DNASE_BED \
    --blacklist $BLACKLIST_BED \
    --prom $PROM \
    --enh $ENH \
    --reg2map $REG2MAP \
    --meta $ROADMAP_META \
    --fastq1 $5 \
    --fastq2 $6 \
    --alignedbam $7 \
    --alignmentlog $8 \
    --coordsortbam $9 \
    --duplog ${10} \
    --finalbam ${11} \
    --finalbed ${12} \
    --bigwig ${13} \
    --peaks ${14}
