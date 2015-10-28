# Daniel Kim, CS Foo
# 2015-10-20
# Script to run ataqc, all parts

import matplotlib
matplotlib.use('Agg')

import os
import sys
import pysam
import pybedtools
import metaseq
import subprocess
import multiprocessing
import numpy as np
from matplotlib import pyplot as plt

import argparse
import logging

from base64 import b64encode
from collections import namedtuple
from collections import OrderedDict
from io import BytesIO

from scipy.signal import find_peaks_cwt

from jinja2 import Template

### QC STUFF

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



def get_bowtie_stats(bowtie_alignment_log):
    '''
    From the Bowtie alignment log, get relevant stats and return
    the file in a list format where each line is an element in
    the list. Can be parsed further if desired.
    '''
    logging.info('Reading bowtie alignment log...')
    bowtie_text = ''
    with open(bowtie_alignment_log, 'rb') as fp:
        for line in fp:
            logging.info(line.strip())
            bowtie_text += line
    return bowtie_text


def get_chr_m(sorted_bam_file):
    '''
    Get fraction of reads that are mitochondrial (chr M).
    '''
    logging.info('Getting mitochondrial chromosome fraction...')
    chrom_list = pysam.idxstats(sorted_bam_file)
    tot_reads = 0
    for chrom in chrom_list:
        chrom_stats = chrom.split('\t')
        if chrom_stats[0] == 'chrM':
            chr_m_reads = int(chrom_stats[2])
        tot_reads += int(chrom_stats[2])
    fract_chr_m = float(chr_m_reads) / tot_reads

    # QC check


    return chr_m_reads, fract_chr_m


def get_gc(qsorted_bam_file, reference_fasta, prefix):
    '''
    Uses picard tools (CollectGcBiasMetrics). Note that the reference
    MUST be the same fasta file that generated the bowtie indices.
    Assumes picard was already loaded into space (module add picard-tools)
    '''
    logging.info('Getting GC bias...')
    output_file = '{0}_gc.txt'.format(prefix)
    plot_file = '{0}_gcPlot.pdf'.format(prefix)
    summary_file = '{0}_gcSummary.txt'.format(prefix)
    get_gc_metrics = ('java -Xmx4G -jar '
                      '/software/picard-tools/1.129/picard.jar ' # environment variable $PICARD instead of hardcoding path
                      'CollectGcBiasMetrics R={0} I={1} O={2} '
                      'VERBOSITY=WARNING QUIET=TRUE '
                      'CHART={3} S={4}').format(reference_fasta,
                                                qsorted_bam_file,
                                                output_file,
                                                plot_file,
                                                summary_file)
    logging.info(get_gc_metrics)
    os.system(get_gc_metrics)
    return output_file, plot_file, summary_file


def plot_gc(data_file):
    '''
    Replot the Picard output as png file to put into the html
    '''
    # Load data
    with open(data_file) as fp:
        for line in fp:
            if line.startswith('## METRICS'):
                break
        data = np.loadtxt(fp, skiprows=1)

    # Plot the data
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlim((0, 100))


    lin1 = ax.plot(data[:, 0], data[:, 4],
                   label='Normalized coverage', color='r')
    ax.set_ylabel('Normalized coverage')

    ax2 = ax.twinx()
    lin2 = ax2.plot(data[:, 0], data[:, 3],
                    label='Mean base quality at GC%', color='b')
    ax2.set_ylabel('Mean base quality at GC%')


    ax3 = ax.twinx()
    # new_fixed_axis = ax3.get_grid_helper().new_fixed_axis
    # ax3.axis['right'] = new_fixed_axis(loc='right', axes=ax3, offset=(60,0))
    lin3 = ax3.plot(data[:, 0], data[:, 1]/np.sum(data[:,1]),
                    label='Windows at GC%', color='g')
    ax3.get_yaxis().set_visible(False)

    lns = lin1 + lin2 + lin3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='best')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return plot_img.getvalue()


def run_preseq(sorted_dups_bam, prefix):
    '''
    Runs preseq. Look at preseq data output to get PBC/NRF.
    '''
    logging.info('Running preseq...')
    preseq_data = '{0}.preseq.dat'.format(prefix)
    preseq_log = '{0}.preseq.log'.format(prefix)
    preseq = ('preseq lc_extrap '
              '-P -B -o {0} {1} -v 2> {2}').format(preseq_data,
                                                   sorted_dups_bam,
                                                   preseq_log)
    logging.info(preseq)
    os.system(preseq)
    return preseq_data, preseq_log

def get_encode_complexity_measures(preseq_log):
    '''
    Use info from the preseq log to calculate NRF, PBC1, and PBC2
    '''
    with open(preseq_log, 'rb') as fp:
        for line in fp:
            if line.startswith('TOTAL READS'):
                tot_reads = float(line.strip().split("= ")[1])
            elif  line.startswith('DISTINCT READS'):
                distinct_reads = float(line.strip().split('= ')[1])
            elif line.startswith('1\t'):
                one_pair = float(line.strip().split()[1])
            elif line.startswith('2\t'):
                two_pair = float(line.strip().split()[1])

    NRF = distinct_reads/tot_reads
    PBC1 = one_pair/distinct_reads
    PBC2 = one_pair/two_pair

    # QC check

    return NRF, PBC1, PBC2

def preseq_plot(data_file):
    '''
    Generate a preseq plot
    '''
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


def make_vplot(bam_file, tss, prefix, bins=400, bp_edge=2000, processes=8,
               greenleaf_norm=True):
    '''
    Take bootstraps, generate V-plots, and get a mean and
    standard deviation on the plot. Produces 2 plots. One is the
    aggregation plot alone, while the other also shows the signal
    at each TSS ordered by strength.
    '''
    logging.info('Generating vplot...')
    vplot_file = '{0}_vplot.png'.format(prefix)
    vplot_large_file = '{0}_large_vplot.png'.format(prefix)

    # Load the TSS file
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, genome='hg19')

    # Load the bam file
    bam = metaseq.genomic_signal(bam_file, 'bam')
    bam_array = bam.array(tss_ext, bins=bins,
                          processes=processes, stranded=True)

    # Normalization (Greenleaf style): Find the avg height
    # at the end bins and take fold change over that
    if greenleaf_norm:
        # Use enough bins to cover 100 bp on either end
        num_edge_bins = int(100/(2*bp_edge/bins))
        bin_means = bam_array.mean(axis=0)
        avg_noise = (sum(bin_means[:num_edge_bins]) +
                     sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
        bam_array /= avg_noise
    else:
        bam_array /= bam.mapped_read_count() / 1e6

    # Generate a line plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-bp_edge, bp_edge, bins)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')

    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig(vplot_file)

    # Print a more complicated plot with lots of info
    #plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 8
    fig = metaseq.plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(5,10),
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
    logging.info('Running Picard MarkDuplicates...')
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
                return float(line_elems[7])

            if mark > 0:
                mark += 1
    return None


def get_samtools_flagstat(bam_file):
    '''
    Runs samtools flagstat to get read metrics
    '''
    logging.info('samtools flagstat...')
    results = pysam.flagstat(bam_file)
    flagstat = ''
    for line in results:
        logging.info(line.strip())
        flagstat += line
    return flagstat


def get_fract_mapq(bam_file, q=30):
    '''
    Runs samtools view to get the fraction of reads of a certain
    map quality.
    '''
    # Current bug in pysam.view module...
    logging.info('samtools mapq 30...')

    # There is a bug in pysam.view('-c'), so just use subprocess
    num_qreads = int(subprocess.check_output(['samtools',
                                              'view', '-c',
                                              '-q', str(q), bam_file]).strip())
    tot_reads = int(subprocess.check_output(['samtools',
                                             'view', '-c',
                                             bam_file]).strip())
    fract_good_mapq = float(num_qreads)/tot_reads
    return fract_good_mapq


def get_final_read_count(first_bam, last_bam):
    '''
    Get final mapped reads compared to initial reads
    '''
    logging.info('final read counts...')
    # Bug in pysam.view
    num_reads_last_bam  = int(subprocess.check_output(['samtools',
                                                       'view', '-c',
                                                       last_bam]).strip())
    num_reads_first_bam  = int(subprocess.check_output(['samtools',
                                                        'view', '-c',
                                                        first_bam]).strip())
    fract_reads_left = float(num_reads_last_bam)/num_reads_first_bam

    return num_reads_first_bam, num_reads_last_bam, fract_reads_left


def get_insert_distribution(final_bam, prefix):
    '''
    Calls Picard CollectInsertSizeMetrics
    '''
    logging.info('insert size distribution...')
    insert_data = '{0}.inserts.hist_data.log'.format(prefix)
    insert_plot = '{0}.inserts.hist_graph.pdf'.format(prefix)
    graph_insert_dist = ('java -Xmx4G -jar '
                         '/software/picard-tools/1.129/picard.jar '
                         'CollectInsertSizeMetrics '
                         'INPUT={0} OUTPUT={1} H={2} '
                         'VERBOSITY=WARNING QUIET=TRUE '
                         'W=1000 STOP_AFTER=5000000').format(final_bam,
                                                             insert_data,
                                                             insert_plot)
    logging.info(graph_insert_dist)
    os.system(graph_insert_dist)
    return insert_data, insert_plot

def get_fract_reads_in_regions(reads_bed, regions_bed):
    '''
    Function that takes in bed file of reads and bed file of regions and
    gets fraction of reads sitting in said regions
    '''
    reads_bedtool = pybedtools.BedTool(reads_bed)
    regions_bedtool = pybedtools.BedTool(regions_bed)

    reads = regions_bedtool.intersect(reads_bedtool, c=True)

    read_count = 0
    for interval in reads:
        read_count += int(interval[-1])
    fract_reads = float(read_count)/reads_bedtool.count()

    return fract_reads


def get_signal_to_noise(final_bed, dnase_regions, blacklist_regions, \
                        prom_regions, enh_regions):
    '''
    Given region sets, determine whether reads are
    falling in or outside these regions
    '''
    logging.info('signal to noise...')

    # Dnase regions
    fract_dnase = get_fract_reads_in_regions(final_bed, dnase_regions)

    # Blacklist regions
    fract_blacklist = get_fract_reads_in_regions(final_bed, blacklist_regions)

    # Prom regions
    fract_prom = get_fract_reads_in_regions(final_bed, prom_regions)

    # Enh regions
    fract_enh = get_fract_reads_in_regions(final_bed, enh_regions)

    return fract_dnase, fract_blacklist, fract_prom, fract_enh


def track_reads(reads_list, labels):
    '''
    This function takes in read counts for different stages
    and generates a bar chart to show where reads are lost
    '''
    # Initial bam, filters (q30), dups, chrM
    ind = np.arange(len(reads_list))
    width = 0.35

    # reads_list = []

    # For each file listed, get the read count
    # for bam_file in bam_list:
    #     num_reads = int(subprocess.check_output(['samtools',
    #                                              'view', '-c',
    #                                              bam_file]).strip())
    #     reads_list.append(num_reads)

    fig, ax = plt.subplots()
    ax.bar(ind, reads_list, width, color='b')
    ax.set_ylabel('Read count')
    ax.set_title('Reads at each processing step')
    ax.set_xticks(ind+width)
    ax.set_xticklabels(labels)

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

{% macro qc_table(qc_results) -%}
  <table class='qc_table'>
    <thead>
        <tr>
            <th scope='col'>Metric</th>
            <th scope='col'>Result</th>
        </tr>
    </thead>
    <tbody>
    {% for result in qc_results %}
    <tr>
        <td>{{ result.metric }}</td>
        <td {% if not result.qc_pass %} class='fail' {% endif %}>
          {{ result.message }}
        </td>
    </tr>
    {% endfor %}
    </tbody>
  </table>
{%- endmacro %}

<html>

<head>
  <title>{{ sample['name'] }} - ATAqC report</title>
  <style>
  .qc_table{
      font-family:"Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
      font-size:12px;
      width:480px;
      text-align:left;
      border-collapse:collapse;
      margin:20px;
  }

  .qc_table th{
      font-size:14px;
      font-weight:normal;
      background:#8c1515;
      border-top:4px solid #700000;
      border-bottom:1px solid #fff;
      color:white;
      padding:8px;
  }

  .qc_table td{
      background:#f2f1eb;
      border-bottom:1px solid #fff;
      color:black;
      border-top:1px solid transparent;
      padding:8px;
  }

  .qc_table .fail{
      color:#ff0000;
      font-weight:bold;
  }
  </style>
</head>

<body>
  <h2>Basic Information</h2>
  <table class='qc_table'>
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
  <pre>
{{ sample['bowtie_stats'] }}
  </pre>

  <h3>Samtools flagstat</h3>
  <pre>
{{ sample['samtools_flagstat'] }}
  </pre>

  <h3>Mapping quality (Fraction > q30)</h3>
  {{ '{0:.3f}'.format(sample['fract_mapq']) }}

  <h3>Percent duplication (Picard MarkDuplicates)</h3>
  {{ '{0:.3f}'.format(sample['percent_dup']) }}

  <h3>Mitochondrial fraction</h3>
  {{ '{0:.3f}'.format(sample['fraction_chr_m']) }}

  <h3>Final read stats</h3>
  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['final_stats'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>

  <h3>GC bias</h3>
  {{ inline_img(sample['gc_bias']) }}


  <h2>Fragment length distribution</h2>
  {{ inline_img(sample['fraglen_dist']) }}
  {{ qc_table(sample['nucleosomal']) }}


  <h2>Enrichment plots (TSS)</h2>
  {{ inline_img(sample['enrichment_plots']['tss']) }}
  <pre>
  An average TSS enrichment is above 6-7. A strong
  TSS enrichment is above 10.
  </pre>

  <h2>Annotation enrichments</h2>
  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['annot_enrichments'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ '{0:.3f}'.format(value) }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>

  <h2>Library complexity</h2>
  <h3>ENCODE library complexity metrics</h3>
  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['library_complexity'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ '{0:.3f}'.format(value) }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>

  <h3>Yield prediction</h3>
  {{ inline_img(sample['yield_prediction']) }}

  <h2>Read summary</h2>
  {{ inline_img(sample['read_tracker']) }}

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
    parser.add_argument('--prom', help='Promoter region file')
    parser.add_argument('--enh', help='Enhancer region file')

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
    INPUT_PREFIX = os.path.join(args.workdir, args.inprefix)
    OUTPUT_PREFIX = os.path.join(args.outdir, args.outprefix)
    NAME = args.outprefix

    # Set up annotations
    REF = args.ref
    TSS = args.tss
    DNASE = args.dnase
    BLACKLIST = args.blacklist
    PROM = args.prom
    ENH = args.enh

    # If mode 1 - TO BE DEPRECATED. In this case, the module is run with 
    # Jin's pipeline
    if args.pipeline == 'kundajelab':
        ALIGNED_BAM = '{0}.bam'.format(INPUT_PREFIX)
        ALIGNMENT_LOG = '{0}.align.log'.format(INPUT_PREFIX)
        QSORT_BAM = '{0}.sort.bam'.format(INPUT_PREFIX)
        COORDSORT_BAM = '{0}.filt.srt.nodup.bam'.format(INPUT_PREFIX)
        DUP_LOG = '{0}.filt.srt.dup.qc'.format(INPUT_PREFIX)
        FINAL_BAM = '{0}.filt.srt.nodup.nonchrM.bam'.format(INPUT_PREFIX)
        FINAL_BED = '{0}.filt.srt.nodup.nonchrM.tn5.bed.gz'.format(INPUT_PREFIX)
    else: # mode 2
        ALIGNED_BAM = args.alignedbam
        ALIGNMENT_LOG = args.alignmentlog
        QSORT_BAM = args.qsortedbam
        COORDSORT_BAM = args.coordsortbam
        DUP_LOG = args.duplog
        FINAL_BAM = args.finalbam
        FINAL_BED = args.finalbed

    return NAME, OUTPUT_PREFIX, REF, TSS, DNASE, BLACKLIST, PROM, ENH, \
           ALIGNED_BAM, ALIGNMENT_LOG, QSORT_BAM, COORDSORT_BAM, DUP_LOG, \
           FINAL_BAM, FINAL_BED


def main():

    # Parse args
    [ NAME, OUTPUT_PREFIX, REF, TSS, DNASE, BLACKLIST, PROM, ENH, \
      ALIGNED_BAM, ALIGNMENT_LOG, QSORT_BAM, COORDSORT_BAM, \
      DUP_LOG, FINAL_BAM, FINAL_BED ] = parse_args()

    # Set up the log file
    logging.basicConfig(filename='test.log',level=logging.DEBUG)

    # Sequencing metrics: Bowtie1/2 alignment log, chrM, GC bias
    BOWTIE_STATS = get_bowtie_stats(ALIGNMENT_LOG)
    chr_m_reads, fraction_chr_m = get_chr_m(QSORT_BAM)
    gc_out, gc_plot, gc_summary = get_gc(COORDSORT_BAM,
                                         REF,
                                         OUTPUT_PREFIX)

    # Library complexity: Preseq results, NRF, PBC1, PBC2
    # preseq_data, preseq_log = run_preseq(justAlignedQnameSortBam, OUTPREFIX)
    preseq_log = '/srv/scratch/dskim89/tmp/ataqc/KeraFreshDay01minA_1.trim.preseq.log'
    preseq_data = '/srv/scratch/dskim89/tmp/ataqc/KeraFreshDay01minA_1.trim.preseq.dat'
    NRF, PBC1, PBC2 = get_encode_complexity_measures(preseq_log)

    # Filtering metrics: duplicates, map quality
    fract_mapq = get_fract_mapq(ALIGNED_BAM)
    percent_dup = get_picard_dup_stats(DUP_LOG)
    flagstat = get_samtools_flagstat(ALIGNED_BAM)

    # Final read statistics
    first_read_count, final_read_count, \
    fract_reads_left = get_final_read_count(ALIGNED_BAM,
                                            FINAL_BAM)

    # Insert size distribution
    insert_data, insert_plot = get_insert_distribution(FINAL_BAM,
                                                       OUTPUT_PREFIX)

    # Enrichments: V plot for enrichment
    vplot_file, vplot_large_file = make_vplot(COORDSORT_BAM,
                                              TSS, OUTPUT_PREFIX)

    # Signal to noise: reads in DHS regions vs not, reads falling
    # into blacklist regions
    fract_dnase, fract_blacklist, \
    fract_prom, fract_enh = get_signal_to_noise(FINAL_BED,
                                                DNASE,
                                                BLACKLIST,
                                                PROM,
                                                ENH)

    # Also need to run n-nucleosome estimation
    nucleosomal_qc = fragment_length_qc(read_picard_histogram(insert_data))

    # Finally output the bar chart of reads
    read_count_data = [first_read_count, first_read_count*fract_mapq,
                 first_read_count*fract_mapq*(1-float(percent_dup)),
                 final_read_count]
    read_count_labels = ['Start', 'q>30', 'dups removed', 'chrM removed (final)']
    read_tracker_plot = track_reads(read_count_data, read_count_labels)


    # Take all this info and render the html file
    BASIC_INFO = OrderedDict([
        ('Filename', NAME),
        ('Genome', 'hg19'),
    ])

    TEST_ENRICHMENT_PLOTS = {
        'tss': b64encode(open(vplot_large_file, 'rb').read())
    }

    FINAL_BAM_STATS = OrderedDict([
        ('Final read count', '{:,}'.format(final_read_count)),
        ('Fraction of reads after filtering',
         '{0:.3f}'.format(fract_reads_left)),
    ])

    ANNOT_ENRICHMENTS = OrderedDict([
        ('Fraction of reads in universal DHS regions', fract_dnase),
        ('Fraction of reads in blacklist regions', fract_blacklist),
        ('Fraction of reads in promoter regions', fract_prom),
        ('Fraction of reads in enhancer regions', fract_enh),
    ])

    LIBRARY_COMPLEXITY = OrderedDict([
        ('NRF', NRF),
        ('PBC1', PBC1),
        ('PBC2', PBC2),
    ])

    SAMPLE = {
        'name': NAME,
        'basic_info': BASIC_INFO,
        'bowtie_stats': BOWTIE_STATS,
        'fraction_chr_m': fraction_chr_m,
        'gc_bias': b64encode(plot_gc(gc_out)),
        'fract_mapq': fract_mapq,
        'percent_dup': percent_dup,
        'samtools_flagstat': flagstat,
        'final_stats': FINAL_BAM_STATS,
        'annot_enrichments': ANNOT_ENRICHMENTS,
        'enrichment_plots': TEST_ENRICHMENT_PLOTS,
        'library_complexity': LIBRARY_COMPLEXITY,
        'yield_prediction': b64encode(preseq_plot(preseq_data)),
        'fraglen_dist': b64encode(fragment_length_plot(insert_data)),
        'nucleosomal': nucleosomal_qc,
        'read_tracker': b64encode(read_tracker_plot),
    }

    results = open('test.html', 'w')
    results.write(html_template.render(sample=SAMPLE))
    results.close()

    return None

main()
