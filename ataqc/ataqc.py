#!/usr/bin/env python


"""Main executable for ATAqC
"""


import os
import argparse
import json
import logging
import timeit

import reads
import peaks
import integrative
import qc
import viz

def parse_args():
    """Parse arguments to go into ataqc"""
    
    parser = argparse.ArgumentParser(description='ATAC-seq QC package')

    # where to put things
    parser.add_argument('--outprefix', help='Output prefix')

    # run mode (QC for general purpose or ENCODE)
    parser.add_argument('--mode', help='Run QC with all_metrics or encode_metrics', default='all_metrics')
    
    # run information
    parser.add_argument('--genome', help='Genome build used')
    parser.add_argument('--cutoff', help='TSS enrichment quality cutoff value')
    parser.add_argument('--sample_name', help='Sample name')

    # annotation files
    parser.add_argument('--chromsizes', help='chromsizes file')
    parser.add_argument('--ref_fa', help='Reference fasta file')
    parser.add_argument('--tss', help='TSS file')
    parser.add_argument('--dnase', help='Open chromatin region file')
    parser.add_argument('--blacklist', help='Blacklisted region file')
    parser.add_argument('--prom', help='Promoter region file')
    parser.add_argument('--enh', help='Enhancer region file')
    parser.add_argument('--reg2map', help='file with cell type signals')
    parser.add_argument('--meta', help='Roadmap metadata')

    # put in processed files
    parser.add_argument('--fastq1', help='First set of reads if paired end, or the single end reads')
    parser.add_argument('--fastq2', help='Second set of reads if paired end')
    parser.add_argument('--aligned_bam', help='BAM file from the aligner')
    parser.add_argument('--alignment_log', help='Alignment log')
    parser.add_argument('--coordsort_bam', help='BAM file sorted by coordinate')
    parser.add_argument('--dup_log', help='Picard duplicate metrics file')
    parser.add_argument('--pbc_log', help='ENCODE library complexity metrics file')
    parser.add_argument('--final_bam', help='Final filtered BAM file')
    parser.add_argument('--final_bed', help='Final filtered alignments in BED format')
    parser.add_argument('--bigwig', help='Final bigwig')
    parser.add_argument('--peaks', help='Raw Peak file')
    parser.add_argument('--naive_overlap_peaks', default=None, help='Naive overlap peak file')
    parser.add_argument('--idr_peaks', default=None, help='IDR peak file')
    parser.add_argument('--use_sambamba_markdup', action='store_true', help='Use sambamba markdup instead of Picard')

    # OR put in spreadsheet with all input files
    # NOTE: make sure to allow clean version if someone puts in Excel sheet

    # More descriptive help message? Say what files are necessary for the input file set?
    parser.add_argument('--files',
                        default=None, help='Input file set (json or tab delimited)')

    # load species file
    parser.add_argument('--species',
                        default=None, help='Species json with annotation files')

    # ENCODE metrics only
    parser.add_argument('--encode_only',
                        action='store_true', help='only run ENCODE metrics')
    
    args = parser.parse_args()

    return args


def unicode_to_byte(data, ignore_dicts=False):
    # if data is a unicode string, return its string representation
    if isinstance(data, unicode):
        return data.encode('utf-8')
        
    # if data is a list of unicode values, return list of byte-encoded values
    if isinstance(data, list):
        return [unicode_to_byte(element, ignore_dicts=True) for element in data]

    # if data is a dictionary, return dictionary of byte-encoded keys and values,
    # but only if key-value pairs have not already been converted
    if isinstance(data, dict) and not ignore_dicts:
        return {
unicode_to_byte(key, ignore_dicts=True): unicode_to_byte(value, ignore_dicts=True) for key, value in data.iteritems()
        }
    
    return data


def run_ataqc(args):
    """Run all the metrics and save out into a metrics dictionary"""

    # Set up the log file, timing, directories
    logging.basicConfig(filename='{}.ataqc.log'.format(args.outprefix), level=logging.DEBUG)
    os.system('mkdir -p {}'.format(os.path.dirname(args.outprefix)))
    start = timeit.default_timer()

    # insert unicode-to-byte function here. json.load() output is in unicode format, and is not usable

    # by samtools or picard

    # set up species file
    if (args.species is not None):
        with open(args.species, 'r') as fp:
            u_species_files = json.load(fp)

        # convert from unicode to byte here
        species_files = unicode_to_byte(u_species_files)

    else:
        annotations = [args.dnase, args.blacklist, args.prom, args.enh, args.peaks]
        annotation_names = ['DNAse', 'Blacklist', 'Promoters', 'Enhancers', 'Peaks']

        species_files = {'chromsizes': args.chromsizes,
                         'ref_fa': args.ref_fa,
                         'annotations': annotations,
                         'annotation_names': annotation_names,
                         'tss': args.tss,
                         'reg2map': args.reg2map,
                         'meta': args.meta
        }

    # set up inputs file
    if (args.files is not None):
        with open(args.files, 'r') as fp:
            u_data_files = json.load(fp)

        # convert from unicode to byte here
        data_files = unicode_to_byte(u_data_files)

    else:
        peak_files = [args.peaks, args.naive_overlap_peaks, args.idr_peaks]
        peak_names = ['Raw Peaks', 'Naive Overlap Peaks', 'IDR Peaks']
        
        data_files = {'genome': args.genome,
                     'fastq1': args.fastq1,
                     'fastq2': args.fastq2,
                     'final_bed': args.final_bed,
                     'peak_files': peak_files,
                     'peak_names': peak_names,
                     'aligned_bam': args.aligned_bam,
                     'final_bam': args.final_bam,
                     'alignment_log': args.alignment_log,
                     'dup_log': args.dup_log,
                     'pbc_log': args.pbc_log
        }


    metrics = {} # easy to convert this into json or any other format
            
    # TODO: give ENCODE only version
    # run fastq QC (fastqc results)
    # TODO - fastq1, fastq2
    # TODO set up parallel mode to run all these metrics on files at the same time
    

    # run BAM file QC
    raw_aligned_bam = reads.AlignedReads('aligned_bam', data_files, False, species_files, args.outprefix)
    metrics['raw_bam'] = raw_aligned_bam.run_metrics(args.mode)

    final_aligned_bam = reads.AlignedReads('final_bam', data_files, True, species_files, args.outprefix)
    metrics['final_bam'] = final_aligned_bam.run_metrics(args.mode)

    print metrics


# run BED (of reads) file QC
    

    
    # run peaks QC
    metrics['peaks'] = {}
    for peak_file_idx in range(len(data_files['peak_files'])):
        peak_file = peaks.Peaks(data_files['peak_files'][peak_file_idx], data_files['peak_names'][peak_file_idx])
        metrics['peaks'][data_files['peak_names'][peak_file_idx]] = peak_file.run_metrics(args.mode)

    print metrics
        
    # run compare signal tracks (also clean this up - get most variable peaks and only use those - small subset)
    
    

    # run integrative QC
    # this module operates exclusively on the metrics dict
    metrics['integrative'] = integrative.run_metrics(metrics, data_files, species_files, args.mode)


    #print metrics

    # and run QC checks on the data
    # ie, organize into QC groups
    qc_groups = qc.run_qc(metrics, data_files, args.outprefix, args.sample_name, args.mode)

    print qc_groups
    
    # output json and text files
    #viz.write_summary(atac_files, metrics, type='json') # this will take care of align log
    #viz.write_summary(atac_files, metrics, type='text')
    
    return qc_groups


def viz_ataqc(qc_groups, outprefix, sample_name):
    """View results - html or multiQC"""

    # give option of native html

    viz.write_html(qc_groups, outprefix, sample_name)

    # also give option of multiQC
    

    return None


def main():

    args = parse_args()
    qc_groups = run_ataqc(args)
    viz_ataqc(qc_groups, args.outprefix, args.sample_name)



if __name__ == "__main__":
    main()
