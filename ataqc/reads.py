

import matplotlib
matplotlib.use('Agg')

import os
import math
import gzip

import subprocess
import multiprocessing
import logging

import numpy as np
import pandas as pd
import scipy.stats
from scipy.signal import find_peaks_cwt

import pysam
import pybedtools
import metaseq

from matplotlib import pyplot as plt
from matplotlib import mlab
from base64 import b64encode
from io import BytesIO


class Reads():

    
    def __init__(self, fastq_files):
        if len(fastq_files) == 1:
            self.paired_ended = 'Single-ended'
            self.fastq1 = fastq_files[0]
        elif len(fastq_files) == 2:
            self.paired_ended = 'Paired-ended'
            self.fastq1 = fastq_files[0]
            self.fastq2 = fastq_files[1]


    def count(self):
        """count reads in fastq file"""
        with gzip.open(self.fastq1, 'r') as fp:
            for i, l in enumerate(fp):
                pass
        if self.paired_ended == 'Single-ended':
            return i+1
        else:
            return 2*(i+1)

        
    def fastqc():
        """Run fastqc and keep output files"""

        return None

    
    def run_qc():
        """Run all QC on FASTQ files"""
        
        return None

    
        
class AlignedReads():

    
    def __init__(self, key, data_files, is_filtered, species_files, outprefix):
        bam_file = data_files[key]

        assert(bam_file.endswith('.bam'))
        self.bam_file = bam_file
        self.data_files = data_files
        self.is_filtered = is_filtered
        self.species_files = species_files
        self.outprefix = outprefix
        self.paired_ended = self.is_paired()

        
    def is_paired(self):
        """Determine if reads file contains paired ends"""
        
        logging.info("Determining if paired-ended...")
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
            for read in bam.fetch():
                if read.is_paired:
                    logging.info('paired-ended')
                    return "Paired-ended"
                else:
                    logging.info('single-ended')
                    return "Single-ended"
        
                
    def count(self):
        """Get read count"""
        
        logging.info("Getting read count in BAM file")
        return int(subprocess.check_output(['samtools', 'view', '-c', self.bam_file]).strip())


    def read_length(self, sample_num=1000):
        '''
        Calculate read length from sample of reads
        '''
        logging.info('Getting read length from BAM file...')
        
        sample_count = 0
        total = 0
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
            for read in bam.fetch():
                if sample_count < sample_num:
                    total += len(read.seq)
                    sample_count += 1
                else:
                    break
        read_len = math.ceil(float(total) / sample_num)

        return read_len


    def mapped_count(self):
        """Runs samtools flagstat to get mapped count"""
        logging.info('samtools flagstat...')
        results = pysam.flagstat(self.bam_file, split_lines=True)
        flagstat = ''
        for line in results:
            logging.info(line.strip())
            flagstat += line
            if "mapped" in line and "mate" not in line:
                mapped_reads = int(line.split('+')[0].strip())
        return mapped_reads


    def samtools_flagstat(self):
        """Get flagstat output as string output"""
        logging.info('samtools flagstat...')
        results = pysam.flagstat(self.bam_file, split_lines=True)
        flagstat = ''
        for line in results:
            logging.info(line.strip())
            flagstat += line
        return flagstat

                
                
    def chr_m_stats(self):
        """Get stats on chrM fraction and total count"""

        logging.info('Getting mitochondrial chromosome fraction...')
        chrom_list = pysam.idxstats(self.bam_file, split_lines=True)
        tot_reads = 0
        for chrom in chrom_list:
            chrom_stats = chrom.split('\t')
            if chrom_stats[0] == 'chrM':
                chr_m_reads = int(chrom_stats[2])
                tot_reads += int(chrom_stats[2])
        fract_chr_m = float(chr_m_reads) / tot_reads
            
        return (chr_m_reads, fract_chr_m)

    
    def mapq_stats(self, q=30):
        """Runs samtools view to get the fraction of reads of a certain map quality."""
        
        logging.info('samtools mapq 30...')
        # There is a bug in pysam.view('-c'), so just use subprocess
        num_qreads = int(subprocess.check_output(['samtools',
                                                  'view', '-c',
                                                  '-q', str(q), self.bam_file]).strip())
        tot_reads = int(subprocess.check_output(['samtools',
                                                 'view', '-c',
                                                 self.bam_file]).strip())
        fract_good_mapq = float(num_qreads)/tot_reads
        
        return (num_qreads, fract_good_mapq)

    
    def gc_content(self):
        """Get GC content stats, plots"""
        '''
        Uses picard tools (CollectGcBiasMetrics). Note that the reference
        MUST be the same fasta file that generated the bowtie indices.
        Assumes picard was already loaded into space (module add picard-tools)
        '''
        logging.info('Getting GC bias...')
        output_file = '{0}_gc.txt'.format(self.outprefix)
        plot_file = '{0}_gcPlot.pdf'.format(self.outprefix)
        summary_file = '{0}_gcSummary.txt'.format(self.outprefix)
        get_gc_metrics = ('java -Xmx4G -jar '
                          '{5}/picard.jar '
                          'CollectGcBiasMetrics R={0} I={1} O={2} '
                          'VERBOSITY=ERROR QUIET=TRUE '
                          'ASSUME_SORTED=FALSE '
                          'CHART={3} S={4}').format(self.species_files['ref_fa'],
                                                    self.bam_file,
                                                    output_file,
                                                    plot_file,
                                                    summary_file,
                                                    os.environ['PICARDROOT'])
        logging.info(get_gc_metrics)
        os.system(get_gc_metrics)
        return output_file, plot_file, summary_file

    
    def picard_complexity(self):
        '''
        Picard EsimateLibraryComplexity
        '''
        out_file = '{0}.picardcomplexity.qc'.format(self.outprefix)
        get_gc_metrics = ('java -Xmx4G -jar '
                          '{2}/picard.jar '
                          'EstimateLibraryComplexity INPUT={0} OUTPUT={1} '
                          'VERBOSITY=ERROR '
                          'QUIET=TRUE').format(self.bam_file,
                                               out_file,
                                               os.environ['PICARDROOT'])
        os.system(get_gc_metrics)
        
        # Extract the actual estimated library size
        header_seen = False
        est_library_size = 0
        with open(out_file, 'rb') as fp:
            for line in fp:
                if header_seen:
                    est_library_size = int(float(line.strip().split()[-1]))
                    break
                if 'ESTIMATED_LIBRARY_SIZE' in line:
                    header_seen = True
                    
        return est_library_size


    def preseq_complexity(self):
        '''
        Runs preseq.
        '''
        # First sort because this file no longer exists...
        sort_bam = 'samtools sort -o {1}.sorted.bam -T {1} -@ 2 {0}'.format(
            self.bam_file, self.outprefix)
        os.system(sort_bam)

        logging.info('Running preseq...')
        preseq_data = '{0}.preseq.dat'.format(self.outprefix)
        preseq_log = '{0}.preseq.log'.format(self.outprefix)
        preseq = ('preseq lc_extrap '
                  '-P -B -o {0} {1}.sorted.bam -v 2> {2}').format(preseq_data,
                                                                  self.outprefix,
                                                                  preseq_log)
        logging.info(preseq)
        os.system(preseq)
        os.system('rm {0}.sorted.bam'.format(self.outprefix))
        
        return preseq_data, preseq_log

    
    def encode_complexity(self):
        """Look for PBC log - if doesn't exist, calculate again from file"""
        
        logging.info('Getting ENCODE metrics...')
        if 'pbc_log' in self.data_files.keys():
            with open(self.data_files['pbc_log'], 'rb') as fp:
                for line in fp:
                    l_list = line.strip().split('\t')
                    NRF = float(l_list[4])
                    PBC1 = float(l_list[5])
                    PBC2 = float(l_list[6])
                    break
        else:
            # TODO put in ENCODE calculations here
            pass        
        
        return {'NRF':NRF, 'PBC1':PBC1, 'PBC2':PBC2}


    def picard_dup_stats(self):
        '''
        Parse Picard's MarkDuplicates metrics file
        '''
        logging.info('Running Picard MarkDuplicates...')
        if 'dup_log' in self.data_files.keys():
            
            mark = 0
            dup_stats = {}
            with open(self.data_files['dup_log']) as fp:
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
                        if self.is_paired() == 'Paired-ended':
                            return (float(line_elems[5]), float(line_elems[7]))
                        else:
                            return (float(line_elems[4]), float(line_elems[7]))

                    if mark > 0:
                        mark += 1
                        
        else:
            # TODO call picard markdups if needed
            pass
        
        return None

    
    def get_sambamba_dup_stats(sambamba_dup_file, paired_status):
        '''
        Parse sambamba markdup's metrics file
        '''
        logging.info('Running sambamba markdup...')

        # TODO fix this too
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

        
    def mito_dups(self, use_sambamba=False):
        '''
        Marks duplicates in the original aligned bam file and then determines
        how many reads are duplicates AND from chrM

        To use sambamba markdup instead of picard MarkDuplicates, set
        use_sambamba to True (default False).
        '''

        out_file = '{0}.dupmark.ataqc.bam'.format(self.outprefix)
        metrics_file = '{0}.dup.ataqc'.format(self.outprefix)

        # Filter bam on the flag 0x002
        tmp_filtered_bam = '{0}.filt.bam'.format(self.outprefix)
        tmp_filtered_bam_prefix = tmp_filtered_bam.replace('.bam', '')
        if self.is_paired() == 'Paired-ended':
            filter_bam = ('samtools view -F 1804 -f 2 -u {0} | '
                          'samtools sort - {1}'.format(self.bam_file, tmp_filtered_bam_prefix))
        else:
            filter_bam = ('samtools view -F 1804 -u {0} | '
                          'samtools sort - {1}'.format(self.bam_file, tmp_filtered_bam_prefix))
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

        return (mito_dups, float(mito_dups) / total_dups)

        
    def tss_enrichment(self, bins=400, bp_edge=2000, processes=8, greenleaf_norm=True):
        '''
        Take bootstraps, generate tss plots, and get a mean and
        standard deviation on the plot. Produces 2 plots. One is the
        aggregation plot alone, while the other also shows the signal
        at each TSS ordered by strength.
        '''
        logging.info('Generating tss plot...')
        
        # TODO calculate read length
        sample_num = 1000
        sample_count = 0
        total = 0
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
            for read in bam.fetch():
                if sample_count < sample_num:
                    total += len(read.seq)
                    sample_count += 1
                else:
                    break
        read_len = math.ceil(float(total) / sample_num)
        print read_len
                    
        tss_plot_file = '{0}_tss-enrich.png'.format(self.outprefix)
        tss_plot_large_file = '{0}_large_tss-enrich.png'.format(self.outprefix)

        # Load the TSS file
        tss = pybedtools.BedTool(self.species_files['tss_enrich'])
        tss_ext = tss.slop(b=bp_edge, g=self.species_files['chrsz'])

        # Load the bam file
        bam = metaseq.genomic_signal(self.bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
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

        fig.savefig(tss_plot_file)

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
        fig.savefig(tss_plot_large_file)

        return tss_plot_file, tss_plot_large_file, tss_point_val


    def fragment_length_distribution(self):
        """Use Picard to collect the fragment length distribution metrics """

        logging.info('insert size distribution...')
        insert_data = '{0}.inserts.hist_data.log'.format(self.outprefix)
        insert_plot = '{0}.inserts.hist_graph.pdf'.format(self.outprefix)
        graph_insert_dist = ('java -Xmx4G -jar '
                             '{3}/picard.jar '
                             'CollectInsertSizeMetrics '
                             'INPUT={0} OUTPUT={1} H={2} '
                             'VERBOSITY=ERROR QUIET=TRUE '
                             'W=1000 STOP_AFTER=5000000').format(self.bam_file,
                                                                 insert_data,
                                                                 insert_plot,
                                                                 os.environ['PICARDROOT'])
        logging.info(graph_insert_dist)
        os.system(graph_insert_dist)

        # nucleosomal QC
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

            percent_nfr = data[:NFR_UPPER_LIMIT].sum() / data.sum()
            percent_nfr_vs_mono_nuc = (
                data[:NFR_UPPER_LIMIT].sum() /
                data[MONO_NUC_LOWER_LIMIT:MONO_NUC_UPPER_LIMIT + 1].sum())
            peaks = find_peaks_cwt(data[:, 1], np.array([25]))

            # Calculate peaks
            print "peaks found for frag len distribution"
            print peaks

            return (percent_nfr, percent_nfr_vs_mono_nuc, peaks)

        percent_nfr, percent_nfr_vs_mono_nuc, peaks = fragment_length_qc(read_picard_histogram(insert_data))

        # insert plot
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

        plot_b64 = fragment_length_plot(insert_data)
        
        return percent_nfr, percent_nfr_vs_mono_nuc, peaks, plot_b64


    def annotation_enrichments():
        """Annotation enrichments given species file annotations"""
        
        return None
    

    def compare_samples():
        """Compare to other available samples for quick cell type spot check"""


        return None


    def run_metrics(self):
        """Go through qc functions and generate metrics, save out into dictionary"""

        metrics = {}

        # always run (whether filtered or not)
        metrics['is_paired'] = self.is_paired()
        metrics['read_count'] = self.count()
        metrics['read_length'] = self.read_length()
        metrics['flagstat'] = self.samtools_flagstat()

        if not self.is_filtered:
            # Duplicates, chrM, mapQ filtering, mito dups
            metrics['mapped_count'] = self.mapped_count()
            metrics['picard_dups'] = self.picard_dup_stats()
            metrics['chrM'] = self.chr_m_stats()
            metrics['chrM_dups'] = self.mito_dups()
            metrics['mapq'] = self.mapq_stats()

            # GC content
            metrics['gc'] = self.gc_content()

            # complexity metrics on BAM file
            metrics['picard_complexity'] = self.picard_complexity()
            metrics['preseq'] = self.preseq_complexity()
            metrics['encode_complexity'] = self.encode_complexity()

        else:
            # TSS enrichment
            metrics['tss_enrich'] = self.tss_enrichment()

            # fragment length dist
            metrics['fragment_len'] = self.fragment_length_distribution()
        
        return metrics


class TagAlign():

    def __init__(self, key, data_files, species_files, outprefix):
        bed_file = data_files[key]

        assert(bed_file.endswith('.tagAlign') or bed_file.endswith('tagAlign.gz') or bed_file.endswith('bed') or bed_file.endswith('bed.gz'))
        self.bed_file = bed_file
        self.data_files = data_files
        self.species_files = species_files
        self.outprefix = outprefix


    def tss_enrichment(self):
        """get TSS enrichment from BED of reads"""
        
        return None

    
    def get_fract_reads_in_regions(self, reads_bed, regions_bed):
        """ Function that takes in bed file of reads and bed file of regions and 
        gets fraction of reads sitting in said regions
        """
        reads_bedtool = pybedtools.BedTool(reads_bed)
        regions_bedtool = pybedtools.BedTool(regions_bed)

        reads = regions_bedtool.sort().merge().intersect(reads_bedtool, c=True, nonamecheck=True)
        
        read_count = 0
        for interval in reads:
            read_count += int(interval[-1])
            fract_reads = float(read_count)/reads_bedtool.count()
                                            
        return read_count, fract_reads
    

    def annotation_enrichments(self):
        """Given annotation files (species file), calculate enrichment of reads"""
        
        return None


    def run_metrics(self):

        return None
    
