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
    
    # put in processed files
    parser.add_argument('--fastq1',
                        help='First set of reads if paired end, \
                              or the single end reads')
    parser.add_argument('--fastq2',
                        help='Second set of reads if paired end')
    parser.add_argument('--alignedbam', help='BAM file from the aligner')
    parser.add_argument('--alignmentlog', help='Alignment log')
    parser.add_argument('--coordsortbam', help='BAM file sorted by coordinate')
    parser.add_argument('--duplog', help='Picard duplicate metrics file')
    parser.add_argument('--pbc', help='ENCODE library complexity metrics file')
    parser.add_argument('--finalbam', help='Final filtered BAM file')
    parser.add_argument('--finalbed',
                        help='Final filtered alignments in BED format')
    parser.add_argument('--bigwig',
                        help='Final bigwig')
    parser.add_argument('--peaks',
                        help='Peak file')
    parser.add_argument('--naive_overlap_peaks',
                        default=None, help='Naive overlap peak file')
    parser.add_argument('--idr_peaks',
                        default=None, help='IDR peak file')
    parser.add_argument('--use_sambamba_markdup', action='store_true',
                        help='Use sambamba markdup instead of Picard')

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


def run_ataqc(args):
    """Run all the metrics and save out into a metrics dictionary"""

    # Set up the log file, timing, directories
    logging.basicConfig(filename='{}.ataqc.log'.format(args.outprefix), level=logging.DEBUG)
    os.system('mkdir -p {}'.format(os.path.dirname(args.outprefix)))
    start = timeit.default_timer()

    # set up species file
    with open(args.species, 'r') as fp:
        species_files = json.load(fp)

    # set up inputs file
    if args.files:
        with open(args.files, 'r') as fp:
            data_files = json.load(fp)

    metrics = {} # easy to convert this into json or any other format
            
    # TODO: give ENCODE only version
    # run fastq QC (fastqc results)
    # TODO - fastq1, fastq2
    # TODO set up parallel mode to run all these metrics on files at the same time
    

    # run BAM file QC
    raw_aligned_bam = reads.AlignedReads('raw_bam', data_files, False, species_files, args.outprefix)
    metrics['raw_bam'] = raw_aligned_bam.run_metrics()

    final_aligned_bam = reads.AlignedReads('final_bam', data_files, True, species_files, args.outprefix)
    metrics['final_bam'] = final_aligned_bam.run_metrics()

    print metrics


# run BED (of reads) file QC
    

    
    # run peaks QC
    metrics['peaks'] = []
    for peak_file_idx in range(len(data_files['peaks'])):
        peak_file = peaks.Peaks(data_files['peaks'][peak_file_idx], data_files['peak_names'][peak_file_idx])
        metrics['peaks'].append(peak_file.run_metrics())

    print metrics
        
    # run compare signal tracks (also clean this up - get most variable peaks and only use those - small subset)
    
    

    # run integrative QC
    # this module operates exclusively on the metrics dict
    metrics['integrative'] = integrative.run_metrics(metrics, data_files, species_files)


    #print metrics

    # and run QC checks on the data
    # ie, organize into QC groups
    qc_groups = qc.run_qc(metrics, data_files, args.outprefix)

    print qc_groups
    
    # output json and text files
    #viz.write_summary(atac_files, metrics, type='json') # this will take care of align log
    #viz.write_summary(atac_files, metrics, type='text')
    
    return qc_groups


def viz_ataqc(qc_groups, outprefix):
    """View results - html or multiQC"""

    # give option of native html

    viz.write_html(qc_groups, outprefix)

    # also give option of multiQC
    

    return None


def main():

    args = parse_args()
    qc_groups = run_ataqc(args)
    viz_ataqc(qc_groups, args.outprefix)



if __name__ == "__main__":
    main()
