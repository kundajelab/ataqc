#!user/bin/env python

from collections import namedtuple
from collections import OrderedDict
import abc
import six
import numpy as np
import pandas as pd
import logging
import os

from matplotlib import pyplot as plt
from io import BytesIO
from base64 import b64encode

# Base class
class QCGroup(object):
    
    # Initialize class attributes
    def __init__(self, metrics, data_files, outprefix):
        self.metrics = metrics
        self.data_files = data_files
        self.outprefix = outprefix
        self.qc = OrderedDict()
        
    def get_name(self):
        pass

    def get_description(self):
        pass

    def get_metrics(self):
        pass

    def run_qc(self):
        pass

    # Returns OrderedDict attribute that contains all the QC metrics
    def get_qc(self):
        return self.qc


class SampleInfo(QCGroup):

    def __init__(self, metrics, data_files, outprefix, sample_name):
        super(SampleInfo, self).__init__(metrics, data_files, outprefix)
        self.sample_name = sample_name

    def get_name(self):
        return "Sample information"
    
    def get_description(self):
        return None
    
    def get_metrics(self):
        table_header = []
        flatten_table = True
        sample_table = OrderedDict([
            ('Sample name', QCNoCheck('Sample name')(self.sample_name)),
            ('Genome', QCNoCheck('Genome')(self.data_files['genome'])),
            ('Paired or single end', QCNoCheck('Paired or single end')(self.metrics['raw_bam']['is_paired'])),
            ('Read length', QCNoCheck('Read length')(self.metrics['raw_bam']['read_length']))
        ])
        
        self.qc['sample_table'] = {'qc': sample_table,
                                   'type': 'table',
                                   'header': 'Sample Information',
                                   'description': None,
                                   'table_header': table_header,
                                   'flatten': flatten_table}
            
            
class SummaryReadStats(QCGroup):

    def get_name(self):
        return "Read filtering summary"

    def get_description(self):
        description = """
        Note that all these read counts are determined using 'samtools view' - as such,
        these are all reads found in the file, whether one end of a pair or a single
        end read. In other words, if your file is paired end, then you should divide
        these counts by two. Each step follows the previous step; for example, the
        duplicate reads were removed after reads were removed for low mapping quality.
        """
        return description

    def get_metrics(self):
        table_header = []
        flatten_table = True
        read_summary_table = OrderedDict([
            ("Read count from sequencer", QCNoCheck('Read count from sequencer')
                                                   (self.metrics['raw_bam']['read_count'])),
            ("Read count successfully aligned", QCNoCheck('Read count successfully aligned')
                                                         (self.metrics['raw_bam']['mapped_count'])),
            ("Read count after filtering for mapping quality", QCNoCheck('Read count after filtering for mapping quality')
                                                                        (self.metrics['raw_bam']['mapq'][0])),
            ("Read count after removing duplicate reads", QCNoCheck('Read count after removing duplicate reads')
                                                                   (int(self.metrics['raw_bam']['mapq'][0] - 
                                                                    self.metrics['raw_bam']['picard_dups'][0]))),
            ("Read count after removing mitochondrial reads (final read count)", QCNoCheck('Read count after removing mitochondrial reads (final read count)')
                                                                                          (self.metrics['integrative']['final_reads'][0]))
        ])

        self.qc['read_summary_table'] = {'qc': read_summary_table,
                                         'type': 'table',
                                         'header': 'Summary',
                                         'description': self.get_description(),
                                         'table_header': table_header,
                                         'flatten': flatten_table}

        # TODO make the bar graph
        def track_reads(reads_list, labels):
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

        read_count_vals = [
            self.metrics['raw_bam']['read_count'],
            self.metrics['raw_bam']['read_count']*self.metrics['raw_bam']['mapq'][1],
            self.metrics['raw_bam']['read_count']*self.metrics['raw_bam']['mapq'][1]*(1 - float(self.metrics['raw_bam']['picard_dups'][1])),
            self.metrics['final_bam']['read_count']
        ]


        read_count_labels = ['Start', 'q>30', 'dups removed',
                         'chrM removed (final)']

        read_summary_plot = track_reads(read_count_vals, read_count_labels)

        plot_description = """
        This bar chart also shows the filtering process and where the reads were lost
        over the process. Note that each step is sequential - as such, there may
        have been more mitochondrial reads which were already filtered because of
        high duplication or low mapping quality. Note that all these read counts are
        determined using 'samtools view' - as such, these are all reads found in
        the file, whether one end of a pair or a single end read. In other words,
        if your file is paired end, then you should divide these counts by two.
        """

        img_header = []

        self.qc['read_summary_plot'] = {'qc': read_summary_plot,
                                        'type': 'plot',
                                        'header': None,
                                        'description': plot_description}


class Fingerprints(QCGroup):

    def get_name(self):
        return "BAM Fingerprints"

    def get_description(self):
        pass

    def run_qc(self):
        def get_fingerprints(self):
            logging.info('Analyzing BAM Fingerprints...')

            count_file = '{0}_fingerprints.tab'.format(self.outprefix)
            plot_file = '{0}_fingerprints.png'.format(self.outprefix)
            bam_files = '{0} {1}'.format(self.data_files['aligned_bam'],
                                         self.data_files['final_bam'])
            labels = '{0} {1}'.format('Raw_Bam',
                                      'Final_Bam')

            plot_fingerprints = ('plotFingerprint '
                                 '-b {0} ' 
                                 '--labels {1} '
                                 '--plotFile {2} '
                                 '--outRawCounts {3}').format(bam_files,
                                                              labels,
                                                              plot_file,
                                                              count_file)

            os.system(plot_fingerprints)
        
        get_fingerprints(self)
        img_src = '{0}_fingerprints.png'.format(self.outprefix)
        fingerprints_plot = b64encode(open(img_src, 'rb').read())

        self.qc['fingerprints_plot'] = {'qc': fingerprints_plot,
                                        'type': 'plot',
                                        'header': 'BAM Fingerprints',
                                        'description': None}


class AlignmentStats(QCGroup):
# TODO: Get bowtie stats and samtool flagstats
    
    def get_name(self):
        return "Alignment Stats"

    def get_description(self):
        pass

    def get_metrics(self):
        alignment_log = self.data_files['alignment_log']

        def get_alignment_text(alignment_log):
            '''
            Get relevant stats from the alignment log and return
            the file in a list format where each line is an element in
            the list. Can be parsed further if desired.
            '''
            logging.info('Reading alignment log...')
            alignment_text = ''
            with open(alignment_log, 'rb') as fp:
                for line in fp:
                    logging.info(line.strip())
                    alignment_text += line
            return alignment_text

        alignment_text = get_alignment_text(alignment_log)
        raw_bam_flagstat = self.metrics['raw_bam']['flagstat']
        final_bam_flagstat = self.metrics['final_bam']['flagstat']

        self.qc['alignment_log'] = {'qc': alignment_text,
                                    'type': 'log',
                                    'header': 'Alignment Log',
                                    'description': None}

        self.qc['raw_bam_flagstat'] = {'qc': raw_bam_flagstat,
                                       'type': 'log',
                                       'header': 'Raw Bam Samtools Flagstats',
                                       'description': None}

        self.qc['final_bam_flagstat'] = {'qc': final_bam_flagstat,
                                         'type': 'log',
                                         'header': 'Final Bam Samtools Flagstats',
                                         'description': None}


class FilteringStats(QCGroup):

    def get_name(self):
        return "Filtering QC"

    def get_description(self):
        description = """
        Mapping quality refers to the quality of the read being aligned to that
        particular location in the genome. A standard quality score is > 30.
        Duplications are often due to PCR duplication rather than two unique reads
        mapping to the same location. High duplication is an indication of poor
        libraries. Mitochondrial reads are often high in chromatin accessibility
        assays because the mitochondrial genome is very open. A high mitochondrial
        fraction is an indication of poor libraries. Based on prior experience, a
        final read fraction above 0.70 is a good library.
        """
        return description

    def get_metrics(self):
        table_header = ['Metric', 'Raw Count', 'Percent of Reads']
        flatten_table = True
        filter_table = OrderedDict([
            ('Mapping quality > q30 (out of total)', QCNoCheck('Mapping quality > q30 (out of total)')
                                                              (self.metrics['raw_bam']['mapq'])),
            ('Duplicates (after filtering)', QCNoCheck('Duplicates (after filtering)')
                                                      (self.metrics['raw_bam']['picard_dups'])),
            ('Duplicates that are mitochondrial (out of all dups)', QCNoCheck('Duplicates that are mitochondrial (out of all dups)')
                                                                             (self.metrics['raw_bam']['chrM_dups'])),
            ('Mitochondrial reads (out of total)', QCNoCheck('Mitochondrial reads (out of total)')
                                                            (self.metrics['raw_bam']['chrM'])),
            ('Final reads (after all filters)', QCNoCheck('Final reads (after all filters)')
                                                         (self.metrics['integrative']['final_reads']))
        ])

        self.qc['filter_table'] = {'qc': filter_table,
                                   'type': 'metric_fraction_table',
                                   'header': 'Filtering Statistics',
                                   'description': self.get_description(),
                                   'table_header': table_header,
                                   'flatten': flatten_table}

#   EDIT
class FragmentLengthStats(QCGroup):

    def get_name(self):
        return "Fragment Length Statistics"

    def get_description(self):
        description = """
        Open chromatin assays show distinct fragment length enrichments, as the cut
        sites are only in open chromatin and not in nucleosomes. As such, peaks
        representing different n-nucleosomal (ex mono-nucleosomal, di-nucleosomal)
        fragment lengths will arise. Good libraries will show these peaks in a
        fragment length distribution and will show specific peak ratios.
        """
        return description

    def run_qc(self):
        nuc_range_metrics = [('Presence of NFR peak', 20, 90),
                         ('Presence of Mono-Nucleosomal peak', 120, 250),
                         ('Presence of Di-Nucleosomal peak', 300, 500)]

        # Perform QC checks and create distribution plot
        table_header_frag_stats = ['Metric', 'Value', 'QC', 'Description']
        flatten_len_stats_table = True
        fragment_len_stats_table = OrderedDict([
            ('Fraction of reads in NFR', QCGreaterThanEqualCheck('Fraction of Reads in NFR', 0.4)
                                                                (self.metrics['final_bam']['fragment_len'][0])),
            ('NFR / Mono-Nucleosomal Reads', QCGreaterThanEqualCheck('NFR / Mono-Nucleosomal Reads', 2.5)
                                                                    (self.metrics['final_bam']['fragment_len'][1]))
            ])

        table_header_peaks = ['Metric', 'Peak Set', 'Selected Peak', 'QC']
        flatten_len_peak_table = False
        fragment_len_peak_table = OrderedDict([
            ('Presence of NFR', QCHasElementInRange(*nuc_range_metrics[0])
                                                    (self.metrics['final_bam']['fragment_len'][2])),
            ('Presence of Mono-Nucleosomal Peak', QCHasElementInRange(*nuc_range_metrics[1])
                                                                    (self.metrics['final_bam']['fragment_len'][2])),
            ('Presence of Di-Nucleosomal Peak', QCHasElementInRange(*nuc_range_metrics[2])
                                                                    (self.metrics['final_bam']['fragment_len'][2]))
            ])

        fragment_len_dist_plot = self.metrics['final_bam']['fragment_len'][3]

        # Setup dictionary for parsing into HTML
        self.qc['fragment_len_dist_plot'] = {'qc': fragment_len_dist_plot,
                                             'type': 'plot',
                                             'header': 'Fragment Length Distribution',
                                             'description': None}
                                             
        self.qc['fragment_len_peak_table'] = {'qc': fragment_len_peak_table,
                                              'type': 'table',
                                              'header': 'Nucleosomal Fragment Peaks',
                                              'description': None,
                                              'table_header': table_header_peaks,
                                              'flatten': flatten_len_peak_table}

        self.qc['fragment_len_stats_table'] = { 'qc': fragment_len_stats_table,
                                                'type': 'table',
                                                'header': 'Fragment Length Statistics',
                                                'description': self.get_description(),
                                                'table_header': table_header_frag_stats,
                                                'flatten': flatten_len_stats_table}


class PeakStats(QCGroup):
    
    def get_name(self):
        return "Peak Statistics"

    def get_description(self):
        description = """
        For a good ATAC-seq experiment in human, you expect to get 100k-200k peaks
        for a specific cell type.
        """

        return description

    def get_metrics(self):
        table_header = []
        flatten_table = True

        raw_peak_table = OrderedDict([
            ("Min Size", QCNoCheck("Min Size")(self.metrics['peaks']['Raw Peaks']['sizes'][0]['Min size'])),
            ("25 Percentile", QCNoCheck("25 Percentile")(self.metrics['peaks']['Raw Peaks']['sizes'][0]['25 percentile'])),
            ("50 Percentile", QCNoCheck("50 Percentile (Median)")(self.metrics['peaks']['Raw Peaks']['sizes'][0]['50 percentile (median)'])),
            ("75 Percentile", QCNoCheck("75 Percentile")(self.metrics['peaks']['Raw Peaks']['sizes'][0]['75 percentile'])),
            ("Max Size", QCNoCheck("Max Size")(self.metrics['peaks']['Raw Peaks']['sizes'][0]['Max size'])),
            ("Mean", QCNoCheck("Mean Size")(self.metrics['peaks']['Raw Peaks']['sizes'][0]['Mean']))
            ])

        raw_peak_plot = self.metrics['peaks']['Raw Peaks']['sizes'][1]


        naive_peak_table = OrderedDict([
            ("Min Size", QCNoCheck("Min Size")(self.metrics['peaks']['Naive Overlap Peaks']['sizes'][0]['Min size'])),
            ("25 Percentile", QCNoCheck("25 Percentile")(self.metrics['peaks']['Naive Overlap Peaks']['sizes'][0]['25 percentile'])),
            ("50 Percentile", QCNoCheck("50 Percentile (Median)")(self.metrics['peaks']['Naive Overlap Peaks']['sizes'][0]['50 percentile (median)'])),
            ("75 Percentile", QCNoCheck("75 Percentile")(self.metrics['peaks']['Naive Overlap Peaks']['sizes'][0]['75 percentile'])),
            ("Max Size", QCNoCheck("Max Size")(self.metrics['peaks']['Naive Overlap Peaks']['sizes'][0]['Max size'])),
            ("Mean", QCNoCheck("Mean Size")(self.metrics['peaks']['Naive Overlap Peaks']['sizes'][0]['Mean']))
            ])

        naive_peak_plot = self.metrics['peaks']['Naive Overlap Peaks']['sizes'][1]


        idr_peak_table = OrderedDict([
            ("Min Size", QCNoCheck("Min Size")(self.metrics['peaks']['IDR Peaks']['sizes'][0]['Min size'])),
            ("25 Percentile", QCNoCheck("25 Percentile")(self.metrics['peaks']['IDR Peaks']['sizes'][0]['25 percentile'])),
            ("50 Percentile", QCNoCheck("50 Percentile (Median)")(self.metrics['peaks']['IDR Peaks']['sizes'][0]['50 percentile (median)'])),
            ("75 Percentile", QCNoCheck("75 Percentile")(self.metrics['peaks']['IDR Peaks']['sizes'][0]['75 percentile'])),
            ("Max Size", QCNoCheck("Max Size")(self.metrics['peaks']['IDR Peaks']['sizes'][0]['Max size'])),
            ("Mean", QCNoCheck("Mean Size")(self.metrics['peaks']['IDR Peaks']['sizes'][0]['Mean']))
            ])

        idr_peak_plot = self.metrics['peaks']['IDR Peaks']['sizes'][1]


        self.qc['raw_peak_table'] = {'qc': raw_peak_table,
                                     'type': 'table',
                                     'header': 'Raw Peak File Statistics',
                                     'description': None,
                                     'table_header': table_header,
                                     'flatten': flatten_table}

        self.qc['raw_peak_plot'] = {'qc': raw_peak_plot,
                                    'type': 'plot',
                                    'header': None,
                                    'description': None}

        self.qc['naive_peak_table'] = {'qc': naive_peak_table,
                                     'type': 'table',
                                     'header': 'Naive Overlap Peak File Statistics',
                                     'description': None,
                                     'table_header': table_header,
                                     'flatten': flatten_table}

        self.qc['naive_peak_plot'] = {'qc': naive_peak_plot,
                                    'type': 'plot',
                                    'header': None,
                                    'description': None}

        self.qc['idr_peak_table'] = {'qc': idr_peak_table,
                                     'type': 'table',
                                     'header': 'IDR Peak File Statistics',
                                     'description': None,
                                     'table_header': table_header,
                                     'flatten': flatten_table}

        self.qc['idr_peak_plot'] = {'qc': idr_peak_plot,
                                    'type': 'plot',
                                    'header': None,
                                    'description': None}

    def run_qc(self):
        table_header = ['Metrics', 'Counts', 'Results', 'QC']
        flatten_table = True
        peak_qc_table = OrderedDict([
            ("Raw Peaks", QCGreaterThanEqualCheck('Raw Peaks', 10000)(self.metrics['peaks']['Raw Peaks']['peak_count'])),
            ("Naive Overlap Peaks", QCGreaterThanEqualCheck('Naive Overlap Peaks', 10000)(self.metrics['peaks']['Naive Overlap Peaks']['peak_count'])),
            ("IDR Peaks", QCGreaterThanEqualCheck('IDR peaks', 10000)(self.metrics['peaks']['IDR Peaks']['peak_count']))
            ])

        self.qc['peak_qc'] = {'qc': peak_qc_table,
                              'type': 'table',
                              'header': 'Peak QC',
                              'description': self.get_description(),
                              'table_header': table_header,
                              'flatten': flatten_table}


class GCBias(QCGroup):

    def get_name(self):
        return "Sequence Quality Metrics"

    def get_description(self):
        description = """
        Open chromatin assays are known to have significant GC bias. Please take this
        into consideration as necessary.
        """

        return description

    def get_metrics(self):

        def plot_gc(img_src):
            '''
            Replot the Picard output as png file to put into the html
            '''
            # Load data
            data = pd.read_table(img_src, comment="#")

            # Plot the data
            fig = plt.figure()
            ax = fig.add_subplot(111)

            plt.xlim((0, 100))

            lin1 = ax.plot(data['GC'], data['NORMALIZED_COVERAGE'],
                           label='Normalized coverage', color='r')
            ax.set_ylabel('Normalized coverage')

            ax2 = ax.twinx()
            lin2 = ax2.plot(data['GC'], data['MEAN_BASE_QUALITY'],
                            label='Mean base quality at GC%', color='b')
            ax2.set_ylabel('Mean base quality at GC%')

            ax3 = ax.twinx()
            lin3 = ax3.plot(data['GC'], data['WINDOWS']/np.sum(data['WINDOWS']),
                            label='Windows at GC%', color='g')
            ax3.get_yaxis().set_visible(False)

            lns = lin1 + lin2 + lin3
            labs = [l.get_label() for l in lns]
            ax.legend(lns, labs, loc='best')

            plot_img = BytesIO()
            fig.savefig(plot_img, format='png')

            return b64encode(plot_img.getvalue())

        img_src = '{0}_gc.txt'.format(self.outprefix)
        gc_bias_plot = plot_gc(img_src)

        self.qc['gc_bias_plot'] = {'qc': gc_bias_plot,
                                   'type': 'plot',
                                   'header': 'GC Bias',
                                   'description': self.get_description()}



class AnnotationQualStats(QCGroup):

    def get_name(self):
        return "Annotation-based Quality Metrics"

    def get_description(self):
        description = """
        Signal to noise can be assessed by considering whether reads are falling into
        known open regions (such as DHS regions) or not. A high fraction of reads
        should fall into the universal (across cell type) DHS set. A small fraction
        should fall into the blacklist regions. A high set (though not all) should
        fall into the promoter regions. A high set (though not all) should fall into
        the enhancer regions. The promoter regions should not take up all reads, as
        it is known that there is a bias for promoters in open chromatin assays.
        """

        return description

    def get_metrics(self):
        img_src = '{0}_large_tss-enrich.png'.format(self.outprefix)
        annot_enrich_plot = b64encode(open(img_src, 'rb').read())

        plot_description = """
        Open chromatin assays should show enrichment in open chromatin sites, such as
        TSS's. An average TSS enrichment is above 6-7. A strong TSS enrichment is
        above 10.
        """

        table_header = ['Fraction of Reads in Region', 'Raw Count of Reads', 'Percent of Reads']
        flatten_table = True
        annot_enrich_table = OrderedDict([
            ("Fraction of reads in universal DHS regions", QCNoCheck("Fraction of reads in universal DHS regions")
                                                                    (self.metrics['integrative']['annotation_enrichments']['DNAse'])),
            ("Fraction of reads in blacklist regions", QCNoCheck("Fraction of reads in blacklist regions")
                                                                (self.metrics['integrative']['annotation_enrichments']['Blacklist'])),
            ("Fraction of reads in promoter regions", QCNoCheck("Fraction of reads in promoter regions")
                                                               (self.metrics['integrative']['annotation_enrichments']['Promoters'])),
            ("Fraction of reads in enhancer regions", QCNoCheck("Fraction of reads in enhancer regions")
                                                               (self.metrics['integrative']['annotation_enrichments']['Enhancers'])),
            ("Fraction of reads in called peak regions", QCNoCheck("Fraction of reads in called peak regions")
                                                                  (self.metrics['integrative']['annotation_enrichments']['Peaks']))
            ])

        self.qc['annot_enrich_plot'] = {'qc': annot_enrich_plot,
                                        'type': 'plot',
                                        'header': 'Enrichment Plot (TSS)',
                                        'description': plot_description}

        self.qc['annot_enrich_table'] = {'qc': annot_enrich_table,
                                         'type': 'table',
                                         'header': 'Annotated Genomic Region Enrichments',
                                         'description': self.get_description(),
                                         'table_header': table_header,
                                         'flatten': flatten_table}


class EncodeStats(QCGroup):

    def __init__(self, metrics, data_files, outprefix, tss_cutoff):
        super(EncodeStats, self).__init__(metrics, data_files, outprefix)
        self.tss_cutoff = tss_cutoff
    
    def get_metrics(self):

        raw_read_count = self.metrics['raw_bam']['read_count']
        aligned_read_count = self.metrics['raw_bam']['mapped_count']
        mapq_filtered_count = self.metrics['raw_bam']['mapq'][0]
        de_duped_count = int(self.metrics['raw_bam']['mapq'][0] - self.metrics['raw_bam']['picard_dups'][0])
        no_mito_count = self.metrics['integrative']['final_reads']

        encode_metrics = self.metrics['raw_bam']['encode_complexity']

        rep_peak_count = self.metrics['peaks']['Raw Peaks']['peak_count']
        idr_peak_count = self.metrics['peaks']['IDR Peaks']['peak_count']

        tss_enrichment_val = self.metrics['final_bam']['tss_enrich'][2]
        tss_cutoff = self.tss_cutoff

        if (self.metrics['raw_bam']['is_paired'] == 'Paired-ended'):
            min_read_count = 50000000
        else:
            min_read_count = 25000000

        nuc_range_metrics = [('Presence of NFR peak', 20, 90),
                         ('Presence of Mono-Nucleosomal peak', 120, 250)]

        # Perform QC checks
        table_header_stats = ['Metric', 'Value', 'QC Passed', 'QC Result']
        flatten_stats_table = True
        stats_table = OrderedDict([
            ('Number of Aligned Reads', QCGreaterThanEqualCheck('Number of Aligned Reads', min_read_count)
                                                               (aligned_read_count)),
            ('Percentage of Mapped Reads (Alignment Rate)', QCGreaterThanEqualCheck('Percentage of Mapped Reads (Alignment Rate)', 0.95)
                                                                                   (aligned_read_count / raw_read_count)),
            ('Non-Redundant Fraction (NRF)', QCGreaterThanEqualCheck('Non-Redundant Fraction (NRF)', 0.9)
                                                    (encode_metrics['NRF'])),
            ('PCR Bottlenecking Coefficient 1 (PBC1)', QCGreaterThanEqualCheck('PCR Bottlenecking Coefficient 1 (PBC1)', 0.9)
                                                                              (encode_metrics['PBC1'])),
            ('PCR Bottlenecking Coefficient 2 (PBC2)', QCGreaterThanEqualCheck('PCR Bottlenecking Coefficient 2 (PBC2)', 3)
                                                                              (encode_metrics['PBC2'])),
            ('Replicated Peak Count', QCGreaterThanEqualCheck('Replicated Peak Count', 150000)
                                                             (rep_peak_count)),
            ('IDR Peak Count', QCGreaterThanEqualCheck('IDR Peak Count', 70000)
                                                      (idr_peak_count)),
            ('TSS Enrichment Value', QCGreaterThanEqualCheck('TSS Enrichment Value', tss_cutoff)
                                                            (tss_enrichment_val)),
            ("Fraction of reads in called peak regions", QCNoCheck("Fraction of reads in called peak regions")
                                                                  (self.metrics['integrative']['annotation_enrichments']['Peaks']))
            ])

        table_header_peaks = ['Metric', 'Peak Set', 'Selected Peak', 'QC']
        flatten_peak_table = False
        fragment_len_peak_table = OrderedDict([
            ('Presence of NFR', QCHasElementInRange(*nuc_range_metrics[0])
                                                    (self.metrics['final_bam']['fragment_len'][2])),
            ('Presence of Mono-Nucleosomal Peak', QCHasElementInRange(*nuc_range_metrics[1])
                                                                    (self.metrics['final_bam']['fragment_len'][2]))
            ])


        # Setup dictionary for parsing into HTML
        self.qc['stats_table'] = {'qc': stats_table,
                                  'type': 'table',
                                  'header': 'Thresholding',
                                  'description': None,
                                  'table_header': table_header_stats,
                                  'flatten': flatten_stats_table}

        self.qc['peak_table'] = {'qc': fragment_len_peak_table,
                                 'type': 'table',
                                 'header': 'Peak Presence',
                                 'description': None,
                                 'table_header': table_header_peaks,
                                 'flatten': flatten_peak_table}
        


QCResult = namedtuple('QCResult', ['metric', 'value', 'qc_pass', 'message'])
INF = float("inf")


class QCCheck(object):
    def __init__(self, metric):
        self.metric = metric

    def check(self, value):
        return True

    def message(self, value, qc_pass):
        return 'OK' if qc_pass else 'Failed'

    def __call__(self, value):
        qc_pass = self.check(value)
        return QCResult(self.metric, value, qc_pass,
                        self.message(value, qc_pass))


class QCIntervalCheck(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCIntervalCheck, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, value):
        return self.lower <= value <= self.upper

    def message(self, value, qc_pass):
        return ('OK - Value {} in range [{}, {}]'.format(value,
                                                        self.lower,
                                                        self.upper) 
        if qc_pass else
                'Out of range [{}, {}]'.format(self.lower,
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
        return [elem for elem in elems
                     if self.lower <= elem <= self.upper]

    def message(self, elems, qc_pass):
        return ('OK - Peak @ {} in peak set {} in range [{}, {}]'.format(qc_pass,
                                                                        elems, 
                                                                        self.lower, 
                                                                        self.upper) 
        if qc_pass else
                'Cannot find element in range [{}, {}]'.format( self.lower, 
                                                                self.upper))

# Convert to list then string to preserve commas that are removed from the array to string conversion,
# and to avoid flattening during the HTML parse (strings aren't 'flattened')
    def __call__(self, value):
        qc_pass = self.check(value)
        return QCResult(self.metric, str(list(value)), qc_pass,    
                        self.message(str(list(value)), qc_pass))

class QCNoCheck(QCCheck):
    def __init__(self, metric):
        super(QCNoCheck, self).__init__(metric)

    def __call__(self, value):
        qc_pass = self.check(value)
        return QCResult(self.metric, value, None, None)


def run_qc(metrics, data_files, outprefix, sample_name, mode='all_metrics'):

    if mode == 'all_metrics':
        sample_info = SampleInfo(metrics, data_files, outprefix, sample_name)
        summary_read_stats = SummaryReadStats(metrics, data_files, outprefix)
        alignment_stats = AlignmentStats(metrics, data_files, outprefix)
        fingerprints_plot = Fingerprints(metrics, data_files, outprefix)
        filtering_stats = FilteringStats(metrics, data_files, outprefix)
        fragment_len_stats = FragmentLengthStats(metrics, data_files, outprefix)
        peak_stats = PeakStats(metrics, data_files, outprefix)
        gc_bias_plot = GCBias(metrics, data_files, outprefix)
        annot_enrich_stats = AnnotationQualStats(metrics, data_files, outprefix)
        
        qc_groups = [sample_info,
                summary_read_stats, 
                alignment_stats,
                fingerprints_plot,
                filtering_stats,
                fragment_len_stats,
                peak_stats,
                gc_bias_plot,
                annot_enrich_stats]

    elif mode == 'encode_metrics':
        encode_stats = EncodeStats(metrics, data_files, outprefix, 5)

        qc_groups = [encode_stats]

    for qc_group in qc_groups:
        qc_group.get_metrics()
        qc_group.run_qc()

    return qc_groups