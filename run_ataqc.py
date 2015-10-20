# Daniel Kim, CS Foo
# 2015-10-15
# Script to run ataqc, all parts

import matplotlib
matplotlib.use('Agg')

import os
import sys
import pysam
import pybedtools
import metaseq
import multiprocessing
import numpy as np
from matplotlib import pyplot as plt

import argparse
import logging
# logging.logging


def get_bowtie_stats(bowtie_alignment_log):
    '''
    From the Bowtie alignment log, get relevant stats and return
    the file in a list format where each line is an element in
    the list. Can be parsed further if desired.
    '''
    print 'Reading bowtie alignment log...'
    bowtie_stats = []
    with open(bowtie_alignment_log, 'rb') as fp:
        for line in fp:
            print line.strip()
            bowtie_stats.append(line)
    return bowtie_stats


def get_chr_m(sorted_bam_file):
    '''
    Get fraction of reads that are mitochondrial (chr M).
    '''
    print 'Getting mitochondrial chromosome fraction...'
    chrom_list = pysam.idxstats(sorted_bam_file)
    tot_reads = 0
    for chrom in chrom_list:
        chrom_stats = chrom.split('\t')
        if chrom_stats[0] == 'chrM':
            chr_m_reads = int(chrom_stats[2])
        tot_reads += int(chrom_stats[2])
    return float(chr_m_reads) / tot_reads


def get_gc(qsorted_bam_file, reference_fasta, prefix):
    '''
    Uses picard tools (CollectGcBiasMetrics). Note that the reference
    MUST be the same fasta file that generated the bowtie indices.
    Assumes picard was already loaded into space (module add picard-tools)
    '''
    print 'Getting GC bias...'
    output_file = '{0}_gc.txt'.format(prefix)
    plot_file = '{0}_gcPlot.pdf'.format(prefix)
    summary_file = '{0}_gcSummary.txt'.format(prefix)
    get_gc_metrics = ('java -Xmx4G -jar '
                      '/software/picard-tools/1.129/picard.jar ' # environment variable $PICARD instead of hardcoding path
                      'CollectGcBiasMetrics R={0} I={1} O={2} '
                      'CHART={3} S={4}').format(reference_fasta,
                                                qsorted_bam_file,
                                                output_file,
                                                plot_file,
                                                summary_file)
    print get_gc_metrics
    os.system(get_gc_metrics)
    return output_file, plot_file, summary_file

def plot_gc(data_file):
    '''
    Replot the Picard output as png file to put into the html
    '''
    
    return None


def run_preseq(sorted_dups_bam, prefix): # PUT THE RESULTS INTO OUTPUT FOLDER
    '''
    Runs preseq. Look at preseq data output to get PBC/NRF.
    '''
    print 'Running preseq...'
    preseq_data = '{0}.preseq.dat'.format(prefix)
    preseq_log = '{0}.preseq.log'.format(prefix)
    preseq = ('preseq lc_extrap '
              '-P -B -o {0} {1} -v 2> {2}').format(preseq_data,
                                                   sorted_dups_bam,
                                                   preseq_log)
    print preseq
    os.system(preseq)
    return preseq_data, preseq_log


def make_vplot(bam_file, tss, prefix, bins=100, processes=8):
    '''
    Take bootstraps, generate V-plots, and get a mean and
    standard deviation on the plot. Produces 2 plots. One is the
    aggregation plot alone, while the other also shows the signal
    at each TSS ordered by strength.
    '''
    print 'Generating vplot...'
    vplot_file = '{0}_vplot.png'.format(prefix)
    vplot_large_file = '{0}_large_vplot.png'.format(prefix)

    # Load the TSS file
    tss = pybedtools.BedTool(tss)
    tss_1kb = tss.slop(b=1000, genome='hg19')

    # Load the bam file
    bam = metaseq.genomic_signal(bam_file, 'bam')
    bam_array = bam.array(tss_1kb, bins=bins,processes=processes)
    bam_array /= bam.mapped_read_count() / 1e6   # Change, to match Greenleaf lab metric

    # Generate a line plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-1000, 1000, 100)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')

    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig(vplot_file)

    # Print a more complicated plot with lots of info
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    fig = metaseq.plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(3,7),
                                   vmin=5, vmax=99, percentile=True,
                                   line_kwargs=dict(color='k', label='All'),
                                   fill_kwargs=dict(color='k',alpha=0.3),
                                   sort_by=bam_array.mean(axis=1))

    # And save the file
    fig.savefig(vplot_large_file)

    return vplot_file, vplot_large_file


def get_picard_dup_stats(picard_dup_file):
    '''
    Parse Picard's MarkDuplicates metrics file
    '''
    print 'Running Picard MarkDuplicates...'
    mark = 0
    dup_stats = {}
    with open(picard_dup_file) as fp:
        for line in fp:
            if '##' in line:
                if 'METRICS CLASS' in line:
                    mark = 1
                continue

            if mark == 2:
                line_elems = line.strip().split('\t')
                dup_stats['PERCENT_DUPLICATION'] = line_elems[7]
                return line_elems[7]

            if mark > 0:
                mark += 1
    return None


def get_samtools_flagstat(bam_file):
    '''
    Runs samtools flagstat to get read metrics
    '''
    print 'samtools flagstat...'
    results = pysam.flagstat(bam_file)
    for line in results:
        print line.strip()
    return results


def get_fract_mapq(bam_file, q=30):
    '''
    Runs samtools view to get the fraction of reads of a certain
    map quality.
    '''
    print 'samtools mapq 30...'
    num_qreads = int(pysam.view('-q', '30', '-c', bam_file)[0].strip())
    tot_reads = int(pysam.view('-c', bam_file)[0].strip()) # Change, this counts alignments not reads
    return float(num_qreads)/tot_reads


def get_final_read_count(first_bam, last_bam):
    '''
    Get final mapped reads compared to initial reads
    '''
    print 'final read counts...'
    num_reads_last_bam = int(pysam.view('-c', last_bam)[0].strip())
    num_reads_first_bam = int(pysam.view('-c', first_bam)[0].strip())
    return num_reads_last_bam, float(num_reads_last_bam)/num_reads_first_bam


def get_insert_distribution(final_bam, prefix):
    '''
    Calls Picard CollectInsertSizeMetrics
    '''
    print 'insert size distribution...'
    insert_data = '{0}.inserts.hist_data.log'.format(prefix)
    insert_plot = '{0}.inserts.hist_graph.pdf'.format(prefix)
    graph_insert_dist = ('java -Xmx4G -jar '
                         '/software/picard-tools/1.129/picard.jar '
                         'CollectInsertSizeMetrics '
                         'INPUT={0} OUTPUT={1} H={2} '
                         'W=1000 STOP_AFTER=5000000').format(final_bam,
                                                             insert_data,
                                                             insert_plot)
    print graph_insert_dist
    os.system(graph_insert_dist)
    return insert_data, insert_plot


def get_signal_to_noise(final_bed, dnase_regions, blacklist_regions):
    '''
    Given region sets, determine whether reads are
    falling in or outside these regions
    '''
    print 'signal to noise...'
    # Load files into pybedtools
    dnase_bedtool = pybedtools.BedTool(dnase_regions)
    blacklist_bedtool = pybedtools.BedTool(blacklist_regions)
    tn5_bedtool = pybedtools.BedTool(final_bed)
    tn5_read_count = tn5_bedtool.count()

    # Dnase regions
    dnase_reads = dnase_bedtool.intersect(tn5_bedtool, c=True)
    dnase_read_count = 0
    for interval in dnase_reads:
        dnase_read_count += int(interval[-1])

    # Blacklist regions
    blacklist_reads = blacklist_bedtool.intersect(tn5_bedtool, c=True)
    blacklist_read_count = 0
    for interval in blacklist_reads:
        blacklist_read_count += int(interval[-1])

    return float(dnase_read_count)/tn5_read_count, \
           float(blacklist_read_count)/tn5_read_count


# ===========================================================

from base64 import b64encode
from collections import namedtuple
from collections import OrderedDict
from io import BytesIO

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import numpy as np
from scipy.signal import find_peaks_cwt

from jinja2 import Template


def preseq_plot(data_file):
    data = np.loadtxt(data_file, skiprows=1)
    data /= 1e6  # scale to millions of reads

    fig = plt.figure()

    # Plot the average expected yield
    plt.plot(data[:, 0], data[:, 1], 'r-')

    # Plot confidence intervals
    ci_lower, = plt.plot(data[:, 0], data[:, 2], 'b--')
    ci_upper, = plt.plot(data[:, 0], data[:, 3], 'b--')
    plt.legend([ci_lower], ['95% confidence interval'], loc=4)

    plt.title('Preseq estimated yield')
    plt.xlabel('Sequenced fragments [ millions ]')
    plt.ylabel('Expected distinct fragments [ millions ]')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return plot_img.getvalue()


def read_picard_histogram(data_file):
    with open(data_file) as fp:
        for line in fp:
            if line.startswith('## HISTOGRAM'):
                break
        data = np.loadtxt(fp, skiprows=1)

    return data


QCResult = namedtuple('QCResult', ['metric', 'qc_pass', 'message'])
INF = float("inf")


class QCCheck(object):
    def __init__(self, metric):
        self.metric = metric

    def check(self, value):
        return True

    def message(self, value, qc_pass):
        return ('{} - OK'.format(value) if qc_pass
                else '{} - Failed'.format(value))

    def __call__(self, value):
        qc_pass = self.check(value)
        return QCResult(self.metric, qc_pass, self.message(value, qc_pass))


class QCIntervalCheck(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCIntervalCheck, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, value):
        return self.lower <= value <= self.upper

    def message(self, value, qc_pass):
        return ('{} - OK'.format(value) if qc_pass else
                '{} out of range [{}, {}]'.format(value, self.lower,
                                                  self.upper))


class QCLessThanEqualCheck(QCIntervalCheck):
    def __init__(self, metric, upper):
        super(QCLessThanEqualCheck, self).__init__(metric, -INF, upper)


class QCGreaterThanEqualCheck(QCIntervalCheck):
    def __init__(self, metric, lower):
        super(QCGreaterThanEqualCheck, self).__init__(metric, lower, INF)


class QCHasElementInRange(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCHasElementInRange, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, elems):
        return (len([elem for elem in elems
                    if self.lower <= elem <= self.upper]) > 0)

    def message(self, elems, qc_pass):
        return ('OK' if qc_pass else
                'Cannot find element in range [{}, {}]'.format(
                    self.lower, self.upper))


def fragment_length_qc(data):
    results = []

    NFR_UPPER_LIMIT = 150
    MONO_NUC_LOWER_LIMIT = 150
    MONO_NUC_UPPER_LIMIT = 300

    # % of NFR vs res
    percent_nfr = data[:NFR_UPPER_LIMIT].sum() / data.sum()
    results.append(
        QCGreaterThanEqualCheck('Fraction of reads in NFR', 0.4)(percent_nfr))

    # % of NFR vs mononucleosome
    percent_nfr_vs_mono_nuc = (
        data[:NFR_UPPER_LIMIT].sum() /
        data[MONO_NUC_LOWER_LIMIT:MONO_NUC_UPPER_LIMIT + 1].sum())
    results.append(
        QCGreaterThanEqualCheck('NFR / mono-nuc reads', 2.5)(
            percent_nfr_vs_mono_nuc))

    # peak locations
    peaks = find_peaks_cwt(data[:, 1], np.array([25]))
    nuc_range_metrics = [('Presence of NFR peak', 20, 90),
                         ('Presence of Mono-Nuc peak', 120, 250),
                         ('Presence of Di-Nuc peak', 300, 500)]
    for range_metric in nuc_range_metrics:
        results.append(QCHasElementInRange(*range_metric)(peaks))

    return results


def fragment_length_plot(data_file, peaks=None):
    data = read_picard_histogram(data_file)

    fig = plt.figure()
    plt.bar(data[:, 0], data[:, 1])
    plt.xlim((0, 1000))

    if peaks:
        peak_vals = [data[peak_x, 1] for peak_x in peaks]
        plt.plot(peaks, peak_vals, 'ro')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return plot_img.getvalue()


html_template = Template("""
{% macro inline_img(base64_img, img_type='png') -%}
    <img src="data:image/{{ img_type }};base64,{{ base64_img }}">
{%- endmacro %}

<html>

<head>
  <title>{{ sample['name'] }} - ATAqC report</title>
</head>

<body>
  <h2>Basic Information</h2>
  <table>
    <tbody>
      {% for field, value in sample['basic_info'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>

  <h2>Alignment statistics</h2>

  <h3>Bowtie alignment log</h3>
  <table>
    <tbody>
      {% for line in sample['bowtie_stats'] %}
      <tr>
        <td>{{ line }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>

  <h3>Mitochondrial fraction</h3>
  {{ sample['fraction_chr_m'] }}

  <h3>GC bias</h3>
  <embed src="data:application/pdf;base64,{{ sample['gc_bias'] }}" width='600' height='600'>

  <h3>Mapping quality (Fraction > q30)</h3>
  {{ sample['fract_mapq'] }}

  <h3>Percent duplication (Picard MarkDuplicates)</h3>
  {{ sample['percent_dup'] }}

  <h3>Samtools flagstat</h3>
  <table>
    <tbody>
      {% for line in sample['samtools_flagstat'] %}
      <tr>
        <td>{{ line }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>

  <h3>Final read stats</h3>
  <table>
    <tbody>
      {% for field, value in sample['final_stats'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>


  <h2>Fragment length distribution</h2>
  {{ inline_img(sample['fraglen_dist']) }}

  <h2>Enrichment plots</h2>
  <h3>TSS enrichment plot</h3>
  {{ inline_img(sample['enrichment_plots']['tss']) }}

  <h3>Annotation enrichments</h3>
  <table>
    <tbody>
      {% for field, value in sample['annot_enrichments'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>

  <h2>Library complexity</h2>

  <h2>Yield prediction</h2>
  {{ inline_img(sample['yield_prediction']) }}
</body>

</html>
""")


# ===========================================================

def parse_args():
    '''
    Set up the package to be run from the command line
    '''

    parser = argparse.ArgumentParser(description='ATAC-seq QC package')

    # Directories and prefixes
    parser.add_argument('--workdir',help='Working directory')
    parser.add_argument('--outdir', help='Output directory')
    parser.add_argument('--outprefix', help='Output prefix')

    # Annotation files
    parser.add_argument('--ref', help='Reference fasta file')
    parser.add_argument('--tss', help='TSS file')
    parser.add_argument('--dnase', help='Open chromatin region file')
    parser.add_argument('--blacklist', help='Blacklisted region file')

    # Choose which mode
    parser.add_argument('--pipeline',
                        default='kundajelab',
                        help='Specific pipeline was used')

    # Mode 1: Provide an input prefix
    parser.add_argument('--inprefix', help='Input file prefix')

    # Mode 2: Define every possible QC file
    parser.add_argument('--alignedbam', help='BAM file from the aligner')
    parser.add_argument('--alignmentlog', help='Alignment log')
    parser.add_argument('--qsortedbam', help='BAM file sorted by QNAME')
    parser.add_argument('--coordsortbam', help='BAM file sorted by coordinate')
    parser.add_argument('--duplog', help='Picard duplicate metrics file')
    parser.add_argument('--finalbam', help='Final filtered BAM file')
    parser.add_argument('--finalbed',
                        help='Final filtered alignments in BED format')

    args = parser.parse_args()

    # Set up all variables
    INPUT_PREFIX = '{0}/{1}'.format(args.workdir, args.inprefix)
    OUTPUT_PREFIX = '{0}/{1}'.format(args.outdir, args.outprefix)
    NAME = args.outprefix

    # Set up annotations
    REF = args.ref
    TSS = args.tss
    DNASE = args.dnase
    BLACKLIST = args.blacklist

    # If mode 1
    if args.pipeline == 'kundajelab':
        ALIGNED_BAM = '{0}.bam'.format(INPUT_PREFIX)
        ALIGNMENT_LOG = '{0}.align.log'.format(INPUT_PREFIX)
        QSORT_BAM = '{0}.sort.bam'.format(INPUT_PREFIX)
        COORDSORT_BAM = '{0}.filt.srt.nodup.bam'.format(INPUT_PREFIX)
        DUP_LOG = '{0}.filt.srt.dup.qc'.format(INPUT_PREFIX)
        FINAL_BAM = '{0}.filt.srt.nodup.nonchrM.bam'.format(INPUT_PREFIX)
        FINAL_BED = '{0}.filt.srt.nodup.nonchrM.tn5.bed.gz'.format(INPUT_PREFIX)
    else:
        ALIGNED_BAM = args.alignedbam
        ALIGNMENT_LOG = args.alignmentlog
        QSORT_BAM = args.qsortedbam
        COORDSORT_BAM = args.coordsortbam
        DUP_LOG = args.duplog
        FINAL_BAM = args.finalbam
        FINAL_BED = args.finalbed



    return NAME, OUTPUT_PREFIX, REF, TSS, DNASE, BLACKLIST, ALIGNED_BAM, \
           ALIGNMENT_LOG, QSORT_BAM, COORDSORT_BAM, DUP_LOG, FINAL_BAM, \
           FINAL_BED


def main():

    # # Set up basic variables
    # WORKDIR='/mnt/lab_data/kundaje/projects/skin/data/fastqs/atac/2015-07-24_freshfrozen/merged/KeraFreshDay01minA'
    # OUTDIR='/srv/scratch/dskim89/tmp/ataqc'
    # FILE_PREFIX= 'KeraFreshDay01minA_1.trim'
    # PREFIX = '{0}/{1}'.format(WORKDIR, FILE_PREFIX)
    # OUTPREFIX = '{0}/{1}'.format(OUTDIR, FILE_PREFIX)
    # REF_FASTA = '/srv/scratch/dskim89/tmp/encodeHg19Male.fa'
    #
    # # Annotation files
    # ANNOTDIR='/users/dskim89/git/ataqc/annotations'
    # dnase_regions = '{0}/{1}'.format(ANNOTDIR, 'reg2map_honeybadger2_dnase_all_p2.bed.gz')
    # blacklist_regions = '{0}/{1}'.format(ANNOTDIR, 'Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz')
    # tss = '/users/dskim89/git/pipelines/atac/tss/parsed_hg19_RefSeq.merged.ANS.bed'
    #
    # # Define the different files needed for QC
    # justAlignedBam = '{0}.bam'.format(PREFIX)
    # justAlignedBamLog = '{0}.align.log'.format(PREFIX)
    # justAlignedQnameSortBam = '{0}.sort.bam'.format(PREFIX)
    # justAlignedCoordSortBam = '{0}.filt.srt.nodup.bam'.format(PREFIX)
    # picardDupMetricsFile = '{0}.filt.srt.dup.qc'.format(PREFIX)
    # finalBam = '{0}.filt.srt.nodup.nonchrM.bam'.format(PREFIX)
    # finalBed = '{0}.filt.srt.nodup.nonchrM.tn5.bed.gz'.format(PREFIX)

    [ NAME, OUTPUT_PREFIX, REF, TSS, DNASE, BLACKLIST, \
      ALIGNED_BAM, ALIGNMENT_LOG, QSORT_BAM, COORDSORT_BAM, \
      DUP_LOG, FINAL_BAM, FINAL_BED ] = parse_args()

    # Sequencing metrics: Bowtie1/2 alignment log, chrM, GC bias
    BOWTIE_STATS = get_bowtie_stats(ALIGNMENT_LOG)
    fraction_chr_m = get_chr_m(QSORT_BAM)
    gc_out, gc_plot, gc_summary = get_gc(COORDSORT_BAM,
                                         REF,
                                         OUTPUT_PREFIX)

    # Library complexity: Preseq results - CS's stuff, to integrate.
    # Get PBC and NRF?
#    run_preseq(justAlignedQnameSortBam, OUTPREFIX)

    # Filtering metrics: duplicates, map quality
    fract_mapq = get_fract_mapq(ALIGNED_BAM)
    percent_dup = get_picard_dup_stats(DUP_LOG)
    flagstat = get_samtools_flagstat(ALIGNED_BAM)

    # Final read statistics
    final_read_count, fract_reads_left = get_final_read_count(ALIGNED_BAM,
                                                              FINAL_BAM)

    # Insert size distribution
    insert_data, insert_plot = get_insert_distribution(FINAL_BAM,
                                                       OUTPUT_PREFIX) # bug here
#     fragment_length_plot(insert_data)
#     #fragment_length_qc(insert_data)
#
    # Enrichments: V plot for enrichment
    vplot_file, vplot_large_file = make_vplot(COORDSORT_BAM,
                                              TSS, OUTPUT_PREFIX)

    # Signal to noise: reads in DHS regions vs not, reads falling
    # into blacklist regions
    fract_dnase, fract_blacklist = get_signal_to_noise(FINAL_BED,
                                                       DNASE,
                                                       BLACKLIST)

    # Take all this info and render the html file
    BASIC_INFO = OrderedDict([
        ('Filename', NAME),
        ('Genome', 'hg19'),
    ])

    TEST_ENRICHMENT_PLOTS = {
        'tss': b64encode(open(vplot_large_file, 'rb').read())
    }

    FINAL_BAM_STATS = OrderedDict([
        ('Final read count', final_read_count),
        ('Fraction of reads after filtering', fract_reads_left),
    ])

    ANNOT_ENRICHMENTS = OrderedDict([
        ('Fraction of reads in universal DHS regions', fract_dnase),
        ('Fraction of reads in blacklist regions', fract_blacklist),
    ])

    TEST_SAMPLE = {
        'name': NAME,
        'basic_info': BASIC_INFO,
        'bowtie_stats': BOWTIE_STATS,
        'fraction_chr_m': fraction_chr_m,
        'gc_bias': open(gc_plot, 'rb').read().encode('base64'),
        'fract_mapq': fract_mapq,
        'percent_dup': percent_dup,
        'samtools_flagstat': flagstat,
        'final_stats': FINAL_BAM_STATS,
        'annot_enrichments': ANNOT_ENRICHMENTS,
        'enrichment_plots': TEST_ENRICHMENT_PLOTS,
        'yield_prediction': b64encode(preseq_plot('/srv/scratch/dskim89/tmp/ataqc/KeraFreshDay01minA_1.trim.preseq.dat')),
        'fraglen_dist': b64encode(fragment_length_plot(insert_data)),
    }

    test_html = open('test.html', 'w')
    test_html.write(html_template.render(sample=TEST_SAMPLE))
    test_html.close()

    return None

main()
