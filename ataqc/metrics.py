import matplotlib
matplotlib.use('Agg')

import abc
from base64 import b64encode
from collections import OrderedDict
import gzip
from io import BytesIO
import logging
import os
import subprocess

import numpy as np
import pandas as pd
from scipy.signal import find_peaks_cwt
import scipy.stats
from matplotlib import pyplot as plt
from matplotlib import mlab

from jinja2 import Template
import metaseq
import pybedtools
import pysam
import six


def standard_metric_html_table(header, description, metric_dict):
    template = Template("""
  {% if header is not none %}
  <h2>{{ header }}</h2>
  {% endif %}
  <table class='qc_table'>
    <thead>
        <tr>
            <th scope='col'>Metric</th>
            <th scope='col'>Value</th>
        </tr>
    </thead>
    <tbody>
      {% for field, value in metric_dict.iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
  {% if description is not none %}
  <pre>{{ description }}</pre>
  {% endif %}
    """)
    return template.render(header=header,
                           description=description,
                           metric_dict=metric_dict)


def checked_metric_html_table(header, description, qc_dict):
    template = Template("""
  {% if header is not none %}
  <h2>{{ header }}</h2>
  {% endif %}
  <table class='qc_table'>
    <thead>
        <tr>
            <th scope='col'>Metric</th>
            <th scope='col'>Value</th>
            <th scope='col'>QC Result</th>
        </tr>
    </thead>
    <tbody>
    {% for result in qc_dict %}
    <tr>
        <td>{{ result.metric }}</td>
        <td>{{ result.value }}</td>
        <td {% if not result.qc_pass %} class='fail' {% endif %}>
          {{ result.message }}
        </td>
    </tr>
    {% endfor %}
    </tbody>
  </table>""")
    return template.render(header=header,
                           description=description,
                           qc_dict=qc_dict)


def inline_image_html(base64_img, img_type='png'):
    template = Template("""
    {% if base64_img == '' %}
        <pre>Metric failed.</pre>
    {% else %}
        <img src="data:image/{{ img_type }};base64,{{ base64_img }}">
    {% endif %}""")

    return template.render(base64_img=base64_img, img_type=img_type)


@six.add_metaclass(abc.ABCMeta)
class MetricGroup():
    def __init__(self, **kwargs):
        var_dict = {
            field: kwargs.get(field, None)
            for field in self._init_args
        }
        self.__dict__.update(var_dict)

    @abc.abstractmethod
    def _init_args(self):
        """Returns list of init args that are required and extracted"""
        # Implement as a class level variable
        pass

    @abc.abstractmethod
    def name(self):
        """Returns name; implement as property (or var) in subclass"""
        pass

    @abc.abstractmethod
    def description(self):
        """Returns description; implement as property (or var) in subclass"""
        pass

    @abc.abstractmethod
    def compute(self):
        """Compute the metric(s) and store the values internally.

        Storing metrics in self._dict will enable use of default to_dict()
        """
        pass

    def to_dict(self):
        """Returns the metric(s) as a dict"""
        return self._dict.copy()

    def to_html(self):
        """Returns html for display in the report"""
        return standard_metric_html_table(self.name, self.description,
                                          self.to_dict(), self._image_fields)


class SampleInfo(MetricGroup):
    name = 'Sample Information'
    description = None

    _init_args = ['alignedbam', 'fastq1']

    def compute(self):
        self._dict = {
            'Paired/Single-ended': determine_paired(self.alignedbam),
            'Read length': get_read_length(self.fastq1),
        }


class SummaryReadStats(MetricGroup):
    name = 'Summary Read Statistics'
    description = '''
This bar chart shows the filtering process and where the reads were lost
over the process. Note that each step is sequential - as such, there may
have been more mitochondrial reads which were already filtered because of
high duplication or low mapping quality. Note that all these read counts are
determined using 'samtools view' - as such, these are all reads found in the
file, whether one end of a pair or a single end read. In other words, if your
file is paired end, then you should divide these counts by two.'''

    _init_args = ['fastq1', 'alignedbam', 'finalbam']

    def compute(self):
        self._dict = {
            'Reads in FASTQ': count_reads_fastq(self.fastq1),
            'Aligned reads': count_reads(self.alignedbam),
            'Reads after filtering for mapping quality': ,
            'Reads after removing duplicates':,
            'Reads after removing mitochondrial reads (final)':
                count_reads(self.finalbam),
        }
        pass

    def to_html(self):
        # Add comma formatting to metrics
        # Add image
        formatted_metrics = {
            name: '{0:,}'.format(value)
            for name, value in self.to_dict()
        }
        return standard_metric_html_table(self.name, self.description,
                                          formatted_metrics, self.image_fields)


class AligmentStats(MetricGroup):
    name = 'Alignment Statistics'
    description = ''

    def compute(self):
        pass


class SamtoolsFlagstat(MetricGroup):
    name = 'Samtools Flagstats'
    description = ''

    def compute(self):
        pass


class FilteringStats(MetricGroup):
    name = 'Filtering Statistics'
    description = '''
Mapping quality refers to the quality of the read being aligned to that
particular location in the genome. A standard quality score is > 30.
Duplications are often due to PCR duplication rather than two unique reads
mapping to the same location. High duplication is an indication of poor
libraries. Mitochondrial reads are often high in chromatin accessibility
assays because the mitochondrial genome is very open. A high mitochondrial
fraction is an indication of poor libraries. Based on prior experience, a
final read fraction above 0.70 is a good library.'''

    def compute(self):
        pass

    def to_html(self):
        # TODO: Format metric so number and frac appear next to each other
        #<td>{{ '{0:,}'.format(value[0]) }}</td>
        #<td>{{ '{0:.3f}'.format(value[1]) }}</td>
        formatted_metrics = {
            name: '{0:,}'.format(value)
            for name, value in self.to_dict()
        }
        return standard_metric_html_table(self.name, self.description,
                                          formatted_metrics, self.image_fields)


class LibraryComplexityENCODEStats(MetricGroup):
    name = 'Library Complexity Metrics'
    description = {
        'encode': '''
The non-redundant fraction (NRF) is the fraction of non-redundant mapped reads
in a dataset; it is the ratio between the number of positions in the genome
that uniquely mapped reads map to and the total number of uniquely mappable
reads. The NRF should be > 0.8. The PBC1 is the ratio of genomic locations
with EXACTLY one read pair over the genomic locations with AT LEAST one read
pair. PBC1 is the primary measure, and the PBC1 should be close to 1.
Provisionally 0-0.5 is severe bottlenecking, 0.5-0.8 is moderate bottlenecking,
0.8-0.9 is mild bottlenecking, and 0.9-1.0 is no bottlenecking. The PBC2 is
the ratio of genomic locations with EXACTLY one read pair over the genomic
locations with EXACTLY two read pairs. The PBC2 should be significantly
greater than 1.''',
        'preseq': '''
Preseq performs a yield prediction by subsampling the reads, calculating the
number of distinct reads, and then extrapolating out to see where the
expected number of distinct reads no longer increases. The confidence interval
gives a gauge as to the validity of the yield predictions.'''
    }

    _image_fields = ['preseq_yield_prediction']

    def compute(self):
        pass

    def to_html(self):
        # Need to separate out the ENCODE, Picard and Preseq metrics
        pass


class FragmentLengthStats(MetricGroup):
    name = 'Fragment Length Statistics'
    description = '''
Open chromatin assays show distinct fragment length enrichments, as the cut
sites are only in open chromatin and not in nucleosomes. As such, peaks
representing different n-nucleosomal (ex mono-nucleosomal, di-nucleosomal)
fragment lengths will arise. Good libraries will show these peaks in a
fragment length distribution and will show specific peak ratios.'''

    _image_fields = ['fraglen_dist']

    def compute(self):
        pass

    def to_html(self):
        pass


class PeakStats(MetricGroup):
    name = 'Peak Statistics'
    description = '''
For a good ATAC-seq experiment in human, you expect to get 100k-200k peaks
for a specific cell type'''

    def compute(self):
        pass

    def to_html(self):
        pass


class GCBias(MetricGroup):
    name = 'GC Bias'
    description = '''
Open chromatin assays are known to have significant GC bias. Please take this
into consideration as necessary.
    '''

    def compute(self):
        pass

    def to_html(self):
        pass


class AnnotationStats(MetricGroup):
    name = 'Annotation-based Quality Metrics'
    description = {
        'tss': '''
Open chromatin assays should show enrichment in open chromatin sites, such as
TSS's. An average TSS enrichment is above 6-7. A strong TSS enrichment is
above 10.''',
        'enrichments': '''
Signal to noise can be assessed by considering whether reads are falling into
known open regions (such as DHS regions) or not. A high fraction of reads
should fall into the universal (across cell type) DHS set. A small fraction
should fall into the blacklist regions. A high set (though not all) should
fall into the promoter regions. A high set (though not all) should fall into
the enhancer regions. The promoter regions should not take up all reads, as
it is known that there is a bias for promoters in open chromatin assays.''',
        'roadmap': '''
This bar chart shows the correlation between the Roadmap DNase samples to
your sample, when the signal in the universal DNase peak region sets are
compared. The closer the sample is in signal distribution in the regions
to your sample, the higher the correlation.''',
    }

    def compute(self):
        pass

    def to_html(self):
        # TODO: Format metric so number and frac appear next to each other
        #<td>{{ '{0:,}'.format(value[0]) }}</td>
        #<td>{{ '{0:.3f}'.format(value[1]) }}</td>
        pass


def get_file_handle(filename, mode="r"):
    if filename.endswith('gz') or filename.endswith('.gzip'):
        if (mode == "r"):
            mode = "rb"
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


# CS: should reimplement to only look at first 100-1000 reads
def determine_paired(bam_file):
    '''
    Quick function to determine if the BAM file is paired end or single end
    '''
    num_paired_reads = int(subprocess.check_output(['samtools',
                                                    'view', '-f', '0x1',
                                                    '-c', bam_file]).strip())
    if num_paired_reads > 1:
        return "Paired-ended"
    else:
        return "Single-ended"


# CS: could get away with prob 10k or even 1k reads
def get_read_length(fastq_file):
    '''
    Get read length out of fastq file
    '''
    total_reads_to_consider = 1000000
    line_num = 0
    total_reads_considered = 0
    max_length = 0
    with get_file_handle(fastq_file, 'rb') as fp:
        for line in fp:
            if line_num % 4 == 1:
                if len(line.strip()) > max_length:
                    max_length = len(line.strip())
                total_reads_considered += 1
            if total_reads_considered >= total_reads_to_consider:
                break
            line_num += 1

    return int(max_length)


# CS: parse this file to get fine grained stats
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
    chrom_list = pysam.idxstats(sorted_bam_file, split_lines=True)
    tot_reads = 0
    for chrom in chrom_list:
        chrom_stats = chrom.split('\t')
        if chrom_stats[0] == 'chrM':
            chr_m_reads = int(chrom_stats[2])
        tot_reads += int(chrom_stats[2])
    fract_chr_m = float(chr_m_reads) / tot_reads

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
                      '{5}/picard.jar '
                      'CollectGcBiasMetrics R={0} I={1} O={2} '
                      'VERBOSITY=ERROR QUIET=TRUE '
                      'ASSUME_SORTED=FALSE '
                      'CHART={3} S={4}').format(reference_fasta,
                                                qsorted_bam_file,
                                                output_file,
                                                plot_file,
                                                summary_file,
                                                os.environ['PICARDROOT'])
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
    lin3 = ax3.plot(data[:, 0], data[:, 1]/np.sum(data[:, 1]),
                    label='Windows at GC%', color='g')
    ax3.get_yaxis().set_visible(False)

    lns = lin1 + lin2 + lin3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='best')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return b64encode(plot_img.getvalue())


def run_preseq(bam_w_dups, prefix):
    '''
    Runs preseq. Look at preseq data output to get PBC/NRF.
    '''
    # First sort because this file no longer exists...
    sort_bam = 'samtools sort -o {1}.sorted.bam -T {1} -@ 2 {0}'.format(
        bam_w_dups, prefix)
    os.system(sort_bam)

    logging.info('Running preseq...')
    preseq_data = '{0}.preseq.dat'.format(prefix)
    preseq_log = '{0}.preseq.log'.format(prefix)
    preseq = ('preseq lc_extrap '
              '-P -B -o {0} {1}.sorted.bam -v 2> {2}').format(preseq_data,
                                                              prefix,
                                                              preseq_log)
    logging.info(preseq)
    os.system(preseq)
    os.system('rm {0}.sorted.bam'.format(prefix))
    return preseq_data, preseq_log


def get_encode_complexity_measures(pbc_output):
    '''
    Gets the unique count statistics from the filtered bam file,
    which is consistent with ENCODE metrics for ChIP-seq
    '''
    with open(pbc_output, 'rb') as fp:
        for line in fp:
            l_list = line.strip().split('\t')
            NRF = float(l_list[4])
            PBC1 = float(l_list[5])
            PBC2 = float(l_list[6])
            break

    # QC check
    results = []
    results.append(QCGreaterThanEqualCheck('NRF', 0.8)(NRF))
    results.append(QCGreaterThanEqualCheck('PBC1', 0.8)(PBC1))
    results.append(QCGreaterThanEqualCheck('PBC2', 1.0)(PBC2))

    return results


def get_picard_complexity_metrics(aligned_bam, prefix):
    '''
    Picard EsimateLibraryComplexity
    '''
    out_file = '{0}.picardcomplexity.qc'.format(prefix)
    get_gc_metrics = ('java -Xmx4G -jar '
                      '{2}/picard.jar '
                      'EstimateLibraryComplexity INPUT={0} OUTPUT={1} '
                      'VERBOSITY=ERROR '
                      'QUIET=TRUE').format(aligned_bam,
                                           out_file,
                                           os.environ['PICARDROOT'])
    os.system(get_gc_metrics)

    # Extract the actual estimated library size
    header_seen = False
    est_library_size = 0
    with open(out_file, 'rb') as fp:
        for line in fp:
            if header_seen:
                est_library_size = int(line.strip().split()[-1])
                break
            if 'ESTIMATED_LIBRARY_SIZE' in line:
                header_seen = True

    return est_library_size


def preseq_plot(data_file):
    '''
    Generate a preseq plot
    '''
    try:
        data = np.loadtxt(data_file, skiprows=1)
    except IOError:
        return ''
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

    return b64encode(plot_img.getvalue())


def make_vplot(bam_file, tss, prefix, genome, read_len, bins=400, bp_edge=2000,
               processes=8, greenleaf_norm=True):
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
    tss_ext = tss.slop(b=bp_edge, genome=genome)

    # Load the bam file
    bam = metaseq.genomic_signal(bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
    bam_array = bam.array(tss_ext, bins=bins, shift_width = -read_len/2, # Shift to center the read on the cut site
                          processes=processes, stranded=True)

    # Actually first build an "ends" file
    #get_ends = '''zcat {0} | awk -F '\t' 'BEGIN {{OFS="\t"}} {{if ($6 == "-") {{$2=$3-1; print}} else {{$3=$2+1; print}} }}' | gzip -c > {1}_ends.bed.gz'''.format(bed_file, prefix)
    #print(get_ends)
    #os.system(get_ends)

    #bed_reads = metaseq.genomic_signal('{0}_ends.bed.gz'.format(prefix), 'bed')
    #bam_array = bed_reads.array(tss_ext, bins=bins,
    #                      processes=processes, stranded=True)

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

    # Note the middle high point (TSS)
    tss_point_val = max(bam_array.mean(axis=0))

    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig(vplot_file)

    # Print a more complicated plot with lots of info

    # Find a safe upper percentile - we can't use X if the Xth percentile is 0
    upper_prct = 99
    if mlab.prctile(bam_array.ravel(), upper_prct) == 0.0:
        upper_prct = 100.0

    plt.rcParams['font.size'] = 8
    fig = metaseq.plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(5, 10),
                                   vmin=5, vmax=upper_prct, percentile=True,
                                   line_kwargs=dict(color='k', label='All'),
                                   fill_kwargs=dict(color='k', alpha=0.3),
                                   sort_by=bam_array.mean(axis=1))

    # And save the file
    fig.savefig(vplot_large_file)

    return vplot_file, vplot_large_file, tss_point_val


def get_picard_dup_stats(picard_dup_file, paired_status):
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
                dup_stats['READ_PAIR_DUPLICATES'] = line_elems[5]
                dup_stats['READ_PAIRS_EXAMINED'] = line_elems[2]
                if paired_status == 'Paired-ended':
                    return float(line_elems[5]), float(line_elems[7])
                else:
                    return float(line_elems[4]), float(line_elems[7])

            if mark > 0:
                mark += 1
    return None


def get_sambamba_dup_stats(sambamba_dup_file, paired_status):
    '''
    Parse sambamba markdup's metrics file
    '''
    logging.info('Running sambamba markdup...')
    with open(sambamba_dup_file, 'r') as fp:
        lines = fp.readlines()

    end_pairs = int(lines[1].strip().split()[1])
    single_ends = int(lines[2].strip().split()[1])
    ends_marked_dup = int(lines[4].strip().split()[1])
    if paired_status == 'Paired-ended':
        pairs_marked_dup = 0.5 * float(ends_marked_dup)
        prct_dup = pairs_marked_dup / float(end_pairs)
        return pairs_marked_dup, prct_dup
    else:
        prct_dup = float(ends_marked_dup) / float(single_ends)
        return ends_marked_dup, prct_dup


def get_mito_dups(sorted_bam, prefix, endedness='Paired-ended', use_sambamba=False):
    '''
    Marks duplicates in the original aligned bam file and then determines
    how many reads are duplicates AND from chrM

    To use sambamba markdup instead of picard MarkDuplicates, set
    use_sambamba to True (default False).
    '''

    out_file = '{0}.dupmark.ataqc.bam'.format(prefix)
    metrics_file = '{0}.dup.ataqc'.format(prefix)

    # Filter bam on the flag 0x002
    tmp_filtered_bam = '{0}.filt.bam'.format(prefix)
    tmp_filtered_bam_prefix = tmp_filtered_bam.replace('.bam', '')
    if endedness == 'Paired-ended':
        filter_bam = ('samtools view -F 1804 -f 2 -u {0} | '
                      'samtools sort - {1}'.format(sorted_bam, tmp_filtered_bam_prefix))
    else:
        filter_bam = ('samtools view -F 1804 -u {0} | '
                      'samtools sort - {1}'.format(sorted_bam, tmp_filtered_bam_prefix))
    os.system(filter_bam)

    # Run Picard MarkDuplicates
    mark_duplicates = ('java -Xmx4G -jar '
                       '{0}/picard.jar '
                       'MarkDuplicates INPUT={1} OUTPUT={2} '
                       'METRICS_FILE={3} '
                       'VALIDATION_STRINGENCY=LENIENT '
                       'ASSUME_SORTED=TRUE '
                       'REMOVE_DUPLICATES=FALSE '
                       'VERBOSITY=ERROR '
                       'QUIET=TRUE').format(os.environ['PICARDROOT'],
                                            tmp_filtered_bam,
                                            out_file,
                                            metrics_file)
    if use_sambamba:
        mark_duplicates = ('sambamba markdup -t 8 '
                           '--hash-table-size=17592186044416 '
                           '--overflow-list-size=20000000 '
                           '--io-buffer-size=256 '
                           '{0} '
                           '{1} '
                           '2> {2}').format(tmp_filtered_bam,
                                            out_file,
                                            metrics_file)
    os.system(mark_duplicates)

    # Index the file
    index_file = 'samtools index {0}'.format(out_file)
    os.system(index_file)

    # Get the mitochondrial reads that are marked duplicates
    mito_dups = int(subprocess.check_output(['samtools',
                                             'view', '-f', '1024',
                                             '-c', out_file, 'chrM']).strip())

    total_dups = int(subprocess.check_output(['samtools',
                                              'view', '-f', '1024',
                                              '-c', out_file]).strip())

    # Clean up
    remove_bam = 'rm {0}'.format(out_file)
    os.system(remove_bam)
    remove_metrics_file = 'rm {0}'.format(metrics_file)
    os.system(remove_metrics_file)
    remove_tmp_filtered_bam = 'rm {0}'.format(tmp_filtered_bam)
    os.system(remove_tmp_filtered_bam)

    return mito_dups, float(mito_dups) / total_dups


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
        if "mapped" in line and "mate" not in line:
            mapped_reads = int(line.split('+')[0].strip())
    return flagstat, mapped_reads


# CS: Do we need to compute the fraction of reads here? We already have the
# number of reads (for the aligned bam at least).
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
    return num_qreads, fract_good_mapq


def count_reads(bam):
    cmd_list = ['samtools', 'view', '-c', bam]
    return int(subprocess.check_output(cmd_list).strip())


# CS: deprecate this in favor of count_reads and logic in the class
def get_final_read_count(first_bam, last_bam):
    '''
    Get final mapped reads compared to initial reads
    '''
    logging.info('final read counts...')
    # Bug in pysam.view
    num_reads_last_bam = int(subprocess.check_output(['samtools',
                                                      'view', '-c',
                                                      last_bam]).strip())
    num_reads_first_bam = int(subprocess.check_output(['samtools',
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
                         '{3}/picard.jar '
                         'CollectInsertSizeMetrics '
                         'INPUT={0} OUTPUT={1} H={2} '
                         'VERBOSITY=ERROR QUIET=TRUE '
                         'W=1000 STOP_AFTER=5000000').format(final_bam,
                                                             insert_data,
                                                             insert_plot,
                                                             os.environ['PICARDROOT'])
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

    return read_count, fract_reads


def get_signal_to_noise(final_bed, dnase_regions, blacklist_regions,
                        prom_regions, enh_regions, peaks):
    '''
    Given region sets, determine whether reads are
    falling in or outside these regions
    '''
    logging.info('signal to noise...')

    # Dnase regions
    reads_dnase, fract_dnase = get_fract_reads_in_regions(final_bed,
                                                          dnase_regions)

    # Blacklist regions
    reads_blacklist, \
        fract_blacklist = get_fract_reads_in_regions(final_bed,
                                                     blacklist_regions)

    # Prom regions
    reads_prom, fract_prom = get_fract_reads_in_regions(final_bed,
                                                        prom_regions)

    # Enh regions
    reads_enh, fract_enh = get_fract_reads_in_regions(final_bed, enh_regions)

    # Peak regions
    reads_peaks, fract_peaks = get_fract_reads_in_regions(final_bed, peaks)

    return reads_dnase, fract_dnase, reads_blacklist, fract_blacklist, \
        reads_prom, fract_prom, reads_enh, fract_enh, reads_peaks, \
        fract_peaks


def get_region_size_metrics(peak_file):
    '''
    From the peak file, return a plot of the region size distribution and
    the quartile metrics (summary from R)
    '''

    peak_size_summ = OrderedDict([
        ('Min size', 0),
        ('25 percentile', 0),
        ('50 percentile (median)', 0),
        ('75 percentile', 0),
        ('Max size', 0),
        ('Mean', 0),
    ])

    # If peak file is none, return nothing
    if peak_file == None:
        return peak_size_summ, ''

    # Load peak file. If it fails, return nothing as above
    try:
        peak_df = pd.read_table(peak_file, compression='gzip', header=None)
    except:
        return peak_size_summ, ''

    # Subtract third column from second to get summary
    region_sizes = peak_df.ix[:,2] - peak_df.ix[:,1]

    # Summarize and store in ordered dict
    peak_summary_stats = region_sizes.describe()

    peak_size_summ = OrderedDict([
        ('Min size', peak_summary_stats['min']),
        ('25 percentile', peak_summary_stats['25%']),
        ('50 percentile (median)', peak_summary_stats['50%']),
        ('75 percentile', peak_summary_stats['75%']),
        ('Max size', peak_summary_stats['max']),
        ('Mean', peak_summary_stats['mean']),
    ])

    # Plot density diagram using matplotlib
    fig = plt.figure()
    ax = fig.add_subplot(111)

    y, binEdges = np.histogram(region_sizes, bins=100)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])

    # density = gaussian_kde(y) # from scipy.stats import gaussian_kde
    # density.covariance_factor = lambda : .25
    # density._compute_covariance()

    plt.plot(bincenters, y, '-')
    filename = peak_file.split('/')[-1]
    ax.set_title('Peak width distribution for {0}'.format(filename))
    #ax.set_yscale('log')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return peak_size_summ, b64encode(plot_img.getvalue())


def get_peak_counts(raw_peaks, naive_overlap_peaks=None, idr_peaks=None):
    '''
    Return a table with counts for raw peaks, IDR peaks, and naive
    overlap peaks
    '''

    # Count peaks
    raw_count = sum(1 for line in get_file_handle(raw_peaks))
    if naive_overlap_peaks is not None:
        naive_count = sum(1 for line in get_file_handle(naive_overlap_peaks))
    else:
        naive_count = 0

    if idr_peaks is not None:
        idr_count = sum(1 for line in get_file_handle(idr_peaks))
    else:
        idr_count = 0

    # Literally just throw these into a QC table
    results = []
    results.append(QCGreaterThanEqualCheck('Raw peaks', 10000)(raw_count))
    results.append(QCGreaterThanEqualCheck('Naive overlap peaks',
                                           10000)(naive_count))
    results.append(QCGreaterThanEqualCheck('IDR peaks', 10000)(idr_count))

    return results


def track_reads(reads_list, labels):
    '''
    This function takes in read counts for different stages
    and generates a bar chart to show where reads are lost
    '''
    # Initial bam, filters (q30), dups, chrM
    ind = np.arange(len(reads_list))
    width = 0.35

    fig, ax = plt.subplots()
    ax.bar(ind, reads_list, width, color='b')
    ax.set_ylabel('Read count')
    ax.set_title('Reads at each processing step')
    ax.set_xticks(ind+width)
    ax.set_xticklabels(labels)

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return b64encode(plot_img.getvalue())


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
    try:
        data = read_picard_histogram(data_file)
    except IOError:
        return ''
    except TypeError:
        return ''

    fig = plt.figure()
    plt.bar(data[:, 0], data[:, 1])
    plt.xlim((0, 1000))

    if peaks:
        peak_vals = [data[peak_x, 1] for peak_x in peaks]
        plt.plot(peaks, peak_vals, 'ro')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return b64encode(plot_img.getvalue())


def compare_to_roadmap(bw_file, regions_file, reg2map_file,
                       metadata, output_prefix):
    '''
    Takes a bigwig file and signal file, gets the bwAverageOverBed,
    then compares that signal with the signal in the Roadmap
    regions
    '''

    out_file = '{0}.signal'.format(output_prefix)

    # First get the signal vals for the peak regions
    # remember to use a UCSC formatted bed file for regions
    bw_average_over_bed = 'bigWigAverageOverBed {0} {1} {2}'.format(
                            bw_file, regions_file, out_file)
    logging.info(bw_average_over_bed)
    os.system(bw_average_over_bed)

    # Read the file back in
    sample_data = pd.read_table(out_file, header=None)
    sample_mean0_col = np.array(sample_data.iloc[:, 5])

    # Then, calculate correlations with all other Roadmap samples and rank
    # the correlations
    roadmap_signals = pd.read_table(reg2map_file, compression='gzip')
    (nrow, ncol) = roadmap_signals.shape

    results = pd.DataFrame(columns=('eid', 'corr'))
    for i in range(ncol):
        # Slice, run correlation
        roadmap_i = roadmap_signals.iloc[:, i]
        spearman_corr = scipy.stats.spearmanr(np.array(roadmap_i),
                                              sample_mean0_col)
        results.loc[i] = [roadmap_i.name, spearman_corr[0]]
        logging.info('{0}\t{1}'.format(roadmap_i.name, spearman_corr))

    # Read in metadata to make the chart more understandable
    metadata = pd.read_table(metadata)
    metadata.columns = ['eid', 'mnemonic']

    merged = pd.merge(metadata, results, on='eid')

    sorted_results = merged.sort('corr', ascending=True)

    # Plot results
    pos = np.array(range(ncol)) + 0.5
    fig = plt.figure()
    plt.barh(pos, sorted_results['corr'], align='center', height=1.0)
    plt.yticks(pos, sorted_results['mnemonic'].tolist(), fontsize=7)
    plt.xlabel('Spearmans correlation')
    plt.title('Signal correlation to Roadmap DNase')
    plt.axis('tight')
    ax = plt.axes()
    ax.yaxis.set_ticks_position('none')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png', bbox_inches='tight')
    fig.savefig('test.png', format='png', bbox_inches='tight')

    return b64encode(plot_img.getvalue())


html_template = Template("""
<html>

<head>
  <title>{{ sample_name }} - ATAqC report</title>
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
  <h2>ATAqC</h2>

  {% for group_html in metric_group_html  %}
    {{ group_html }}
  {% endfor %}
</body>

</html>
""")
