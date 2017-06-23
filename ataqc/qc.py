from collections import namedtuple
from collections import OrderedDict
import abc
import six


# QC groups


class QCGroup():
    """Groups QC according to categories"""

    def __init__(self, metrics, data_files):
        self.metrics = metrics
        self.data_files = data_files
        self.qc = OrderedDict()

    def name(self):
        """Returns name"""
        pass
        
    def description(self):
        """Brief description of QC group"""
        pass

    def run_qc(self):
        """Runs QC checks"""
        pass

    def get_qc(self):
        """Get qc object"""
        return self.qc
        

class SampleInfo(QCGroup):

    def get_name(self):
        return "Sample information"

    def get_description(self):
        return None
    
    def run_qc(self):
        sample_table = OrderedDict([
            ('Sample name', self.data_files['sample_name']),
            ('Genome', self.data_files['genome']),
            ('Paired or single end', self.metrics['raw_bam']['is_paired']),
            ('Read length', self.metrics['raw_bam']['read_length'])
        ])
        
        self.qc['sample_table'] = {'qc': sample_table,
                                   'type': 'table',
                                   'header': None,
                                   'description': None}

        
class SummaryReadStats(QCGroup):

    def get_name(self):
        return "Read filtering summary"

    def get_description(self):
        return None
    
    def run_qc(self):
        read_summary_table = OrderedDict([
            ("Read count from sequencer", self.metrics['raw_bam']['read_count']),
            ("Read count successfully aligned", self.metrics['raw_bam']['mapped_count']), 
            ("Read count after filtering for mapping quality", self.metrics['raw_bam']['mapq'][0]),
            ("Read count after removing duplicate reads",
             int(self.metrics['raw_bam']['mapq'][0] - self.metrics['raw_bam']['picard_dups'][0])),
            ("Read count after removing mitochondrial reads (final read count)", self.metrics['integrative']['final_reads'][0])
        ])

        summary_description = """
Note that all these read counts are determined using 'samtools view' - as such,
these are all reads found in the file, whether one end of a pair or a single
end read. In other words, if your file is paired end, then you should divide
these counts by two. Each step follows the previous step; for example, the
duplicate reads were removed after reads were removed for low mapping quality.
        """

        self.qc['read_summary_table'] = {'qc': read_summary_table,
                                         'type': 'metric_table',
                                         'header': None,
                                         'description': summary_description}

        # TODO make the bar graph
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

    read_count_data = [first_read_count, first_read_count*fract_mapq,
                                              first_read_count*fract_mapq*(1-float(percent_dup)),
                                              final_read_count]


    read_count_vals = [
        self.metrics['raw_bam']['read_count'],
        self.metrics['raw_bam']['read_count']*self.metrics['raw_bam']['mapq'][1],
        self.metrics['raw_bam']['read_count']*self.metrics['raw_bam']['mapq'][1]*(1 - float(self.metrics['raw_bam']['picard_dups'][1])),
        self.metrics['final_bam']['read_count']
    ]
    
    read_count_labels = ['Start', 'q>30', 'dups removed',
                         'chrM removed (final)']
    
    read_summary_plot = track_reads([read_count_vals, read_count_labels])

    plot_description = """
This bar chart also shows the filtering process and where the reads were lost
over the process. Note that each step is sequential - as such, there may
have been more mitochondrial reads which were already filtered because of
high duplication or low mapping quality. Note that all these read counts are
determined using 'samtools view' - as such, these are all reads found in
the file, whether one end of a pair or a single end read. In other words,
if your file is paired end, then you should divide these counts by two.
    """
    
    self.qc['read_summary_plot'] = {'qc': read_summary_plot,
                                    'type': 'plot',
                                    'header': None,
                                    'description': plot_description}
    
    
class AlignmentStats(QCGroup):

    def run_qc(self):
        pass


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
final read fraction above 0.70 is a good library."""
        return description
    
    def run_qc(self):
        filter_table = OrderedDict([
            ('Mapping quality > q30 (out of total)', self.metrics['raw_bam']['mapq']),
            ('Duplicates (after filtering)', self.metrics['raw_bam']['picard_dups']),
            ('Mitochondrial reads (out of total)', self.metrics['raw_bam']['chrM']),
            ('Duplicates that are mitochondrial (out of all dups)', self.metrics['raw_bam']['chrM_dups']),
            ('Final reads (after all filters)', self.metrics['integrative']['final_reads'])
        ])

        self.qc['filter_table'] = {'qc': filter_table,
                                   'type': 'metric_fraction_table',
                                   'header': None,
                                   'description': None}
    
    
class EncodeStats(QCGroup):

    def run_qc(self):
        pass


def run_qc(metrics, data_files, subset='all'):
    """groups metrics and runs QC on them"""

    if subset == 'encode':
        qc_group_classes = [SampleInfo, EncodeStats]
    else:
        #qc_groups = [SampleInfo, SummaryReadStats, AlignmentStats, FilteringStats]
        qc_group_classes = [SampleInfo, FilteringStats]

    qc_groups = [group(metrics, data_files) for group in qc_group_classes]
    qc = [group.run_qc() for group in qc_groups]

    return qc_groups



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
        return ('OK' if qc_pass else
                'Out of range [{}, {}]'.format(self.lower, self.upper))


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






# Set up QC groups
