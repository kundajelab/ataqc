# Daniel Kim, CS Foo
# 2015-10-15
# Script to run ataqc, all parts

import matplotlib
matplotlib.use('Agg')

import os, sys
import pysam
import pybedtools
import metaseq
import multiprocessing
import numpy as np
from matplotlib import pyplot as plt

#from make_report import *


def getBowtieStats(bowtieAlignmentLog):
    '''
    From Bowtie log, get relevant stats and return
    num raw reads, percent aligned, etc
    '''
    print 'bowtie alignment log...'
    line_num = 1
    bowtieStats = []
    with open(bowtieAlignmentLog, 'rb') as fp:
        for line in fp:
            print line.strip()
            bowtieStats.append(line.strip())
    return bowtieStats

def getChrM(sortedBamFile):
    '''
    Get fraction of reads that are Chrom M
    '''
    print 'mitochondrial chromosome...'
    chromList = pysam.idxstats(sortedBamFile)
    tot_reads = 0
    for chrom in chromList:
        chrom_stats = chrom.split('\t')
        if chrom_stats[0] == 'chrM':
            chrM_reads = int(chrom_stats[2])
        tot_reads += int(chrom_stats[2])
    return float(chrM_reads) / tot_reads

def getGC(qsortedBamFile, referenceFasta, prefix):
    '''
    Use picard tools. Note that the reference MUST be the same fasta file that generated the
    bowtie indices. Assumes picard was already loaded into space
    '''
    print 'GC bias...'
    getGcMetrics = 'java -Xmx4G -jar /software/picard-tools/1.129/picard.jar CollectGcBiasMetrics R={0} I={1} O={2}_GcBias.txt CHART={2}_GcBiasPlot.pdf S={2}_GcSummary.txt'.format(referenceFasta, qsortedBamFile, prefix)
    print getGcMetrics
    os.system(getGcMetrics)
    return None

def runPreseq(sortedDupsBam, OUTDIR): # PUT THE RESULTS INTO OUTPUT FOLDER
    '''
    Runs preseq. Look at preseq data output to get PBC/NRF?
    '''
    preseqData = '{0}.preseq.dat'.format(sortedDupsBam)
    preseqLog = '{0}.preseq.log'.format(sortedDupsBam)
    preseq = "preseq lc_extrap -P -B -o {0} {1} -v 2> {2}".format(preseqData, sortedDupsBam, preseqLog)
    print preseq
    os.system(preseq)
    return None

def makeVPlot(bamFile, prefix, bins=100, processes=8):
    '''
    Take bootstraps, generate V-plots, and get a mean and standard deviation on the plot
    '''
    print 'vplot...'

    print 'tss...'
    tss = pybedtools.BedTool('/users/dskim89/git/pipelines/atac/tss/parsed_hg19_RefSeq.merged.ANS.bed')
    tss_1kb = tss.slop(b=1000, genome='hg19')

    print 'bams...'
    bam = metaseq.genomic_signal(bamFile, 'bam')
    bam_array = bam.array(tss_1kb, bins=bins,processes=processes)
    bam_array /= bam.mapped_read_count() / 1e6

    print 'making plot...'
    # Generate a line plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-1000, 1000, 100)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')

    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig('{0}_Vplot.pdf'.format(prefix))

    # Print a more complicated plot with lots of info
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    fig = metaseq.plotutils.imshow(bam_array, x=x, figsize=(3,7), vmin=5, vmax=99, percentile=True, line_kwargs=dict(color='k', label='All'), fill_kwargs=dict(color='k',alpha=0.3), sort_by=bam_array.mean(axis=1))

    fig.savefig('{0}_biggerplot.pdf'.format(prefix))

    return None

def getPicardDupStats(picardDupFile):
    '''
    Parse MarkDuplicates metrics file
    '''
    print 'picard markduplicates...'
    mark = 0
    dupStats = {}
    with open(picardDupFile) as fp:
        for line in fp:
            if '##' in line:
                if 'METRICS CLASS' in line:
                    mark = 1
                continue

            if mark == 2:
                line_elems = line.strip().split('\t')
                dupStats['PERCENT_DUPLICATION'] = line_elems[7]
                return dupStats

            if mark > 0:
                mark += 1

    return None

def getSamtoolsFlagstat(bamFile):
    '''
    samtools flagstat output
    '''
    print 'samtools flagstat...'
    results = pysam.flagstat(bamFile)
    for line in results:
        print line.strip()
    return None

def getFractMapQ(bamFile, q=30):
    '''
    samtools mapq reads
    '''
    print 'samtools mapq 30...'
    numQReads = int(pysam.view('-q', '30', '-c', bamFile)[0].strip())
    totReads = int(pysam.view('-c', bamFile)[0].strip()) # Change, this counts alignments not reads

    return float(numQReads)/totReads

def getFinalReadCount(firstBam, lastBam):
    '''
    Get final mapped reads compared to initial reads
    '''
    print 'final read counts...'
    numLastBam = int(pysam.view('-c', lastBam)[0].strip())
    numFirstBam = int(pysam.view('-c', firstBam)[0].strip())
    return numLastBam, float(numLastBam)/numFirstBam

def getInsertDistribution(finalBam):
    '''
    Call picard
    '''
    print 'insert size distribution...'
    graphInsertDist = 'java -Xmx4G -jar /software/picard-tools/1.129/picard.jar CollectInsertSizeMetrics INPUT={0} OUTPUT={0}.hist_data.log H={0}.hist_graph.pdf W=1000 STOP_AFTER=5000000'.format(finalBam)
    print graphInsertDist
    os.system(graphInsertDist)
    return None

def getSignalToNoise(finalBed, dnase_regions, blacklist_regions):
    '''
    Given region sets, determine whether reads are falling in or outside these regions
    '''
    print 'signal to noise...'
    # Load files into pybedtools
    dnaseBedTool = pybedtools.BedTool(dnase_regions)
    blacklistBedTool = pybedtools.BedTool(blacklist_regions)
    tn5BedTool = pybedtools.BedTool(finalBed)
    tn5ReadCount = tn5BedTool.count()

    # Dnase regions
    print 'dnase...'
    dnaseReads = dnaseBedTool.intersect(tn5BedTool, c=True)
    dnaseReadCount = 0

    for interval in dnaseReads:
        dnaseReadCount += int(interval[-1])

    # Blacklist regions
    print 'blacklist...'
    blacklistReads = blacklistBedTool.intersect(tn5BedTool, c=True)
    blacklistReadCount = 0

    for interval in blacklistReads:
        blacklistReadCount += int(interval[-1])

    return dnaseReadCount, blacklistReadCount, tn5ReadCount


def getFinalReadStats():
    '''
    Final useable reads and map ratio (final useable reads over all reads)
    '''

    return None

def main():

    # Set up basic variables
    WORKDIR='/mnt/lab_data/kundaje/projects/skin/data/fastqs/atac/2015-07-24_freshfrozen/merged/KeraFreshDay01minA'
    OUTDIR='/srv/scratch/dskim89/tmp/ataqc'
    FILE_PREFIX= 'KeraFreshDay01minA_1.trim'
    PREFIX = '{0}/{1}'.format(WORKDIR, FILE_PREFIX)
    OUTPREFIX = '{0}/{1}'.format(OUTDIR, FILE_PREFIX)
    REF_FASTA = '/srv/scratch/dskim89/tmp/encodeHg19Male.fa'

    # Annotation files
    ANNOTDIR='/users/dskim89/git/ataqc/annotations'
    dnase_regions = '{0}/{1}'.format(ANNOTDIR, 'reg2map_honeybadger2_dnase_all_p2.bed.gz')
    blacklist_regions = '{0}/{1}'.format(ANNOTDIR, 'Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz')

    # Define the different files needed for QC
    justAlignedBam = '{0}.bam'.format(PREFIX)
    justAlignedBamLog = '{0}.align.log'.format(PREFIX)
    justAlignedQnameSortBam = '{0}.sort.bam'.format(PREFIX)
    justAlignedCoordSortBam = '{0}.filt.srt.nodup.bam'.format(PREFIX)
    picardDupMetricsFile = '{0}.filt.srt.dup.qc'.format(PREFIX)
    finalBam = '{0}.filt.srt.nodup.nonchrM.bam'.format(PREFIX)
    finalBed = '{0}.filt.srt.nodup.nonchrM.tn5.bed.gz'.format(PREFIX)

    # Troubleshooting
    runPreseq(justAlignedQnameSortBam, OUTPREFIX)

    # Seqencing metrics from Bowtie1/2
    bowtieStats = getBowtieStats(justAlignedBamLog)

    # Chr M fraction
    fraction_chrM = getChrM(justAlignedQnameSortBam)

    # GC bias
    getGC(justAlignedCoordSortBam, REF_FASTA, OUTPREFIX)

    # Preseq results - CS's stuff, to integrate. Get PBC and NRF
    runPreseq(justAlignedQnameSortBam)

    # V plot for enrichment
    makeVPlot(justAlignedCoordSortBam, OUTPREFIX)

    # Picard output: library complexity, filter metrics, duplicate metrics
    dupStats = getPicardDupStats(picardDupMetricsFile)

    # Samtools flagstat: library complexity, etc
    getSamtoolsFlagstat(finalBam)

    # percent reads above mapq 30
    fractMapQ = getFractMapQ(justAlignedBam)

    # Final read statistics: num final usable reads, map ratio
    finalReadCount, fractReadsLeft = getFinalReadCount(justAlignedBam, finalBam)

    # Insert size distribution
    getInsertDistribution(finalBam) # bug here

    # Signal to noise: reads in DHS regions vs not, reads falling into blacklist regions
    getSignalToNoise(finalBed, dnase_regions, blacklist_regions)

    return None

main()
