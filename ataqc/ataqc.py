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
    
    # run information
    parser.add_argument('--genome', help='Genome build used')
    parser.add_argument('--sample_name', help='Sample name')

    # annotation files
    parser.add_argument('--chromsizes', help='chromsizes file')
    parser.add_argument('--ref', help='Reference fasta file')
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
    parser.add_argument('--alignedbam', help='BAM file from the aligner')
    parser.add_argument('--alignmentlog', help='Alignment log')
    parser.add_argument('--coordsortbam', help='BAM file sorted by coordinate')
    parser.add_argument('--duplog', help='Picard duplicate metrics file')
    parser.add_argument('--pbc', help='ENCODE library complexity metrics file')
    parser.add_argument('--finalbam', help='Final filtered BAM file')
    parser.add_argument('--finalbed', help='Final filtered alignments in BED format')
    parser.add_argument('--bigwig', help='Final bigwig')
    parser.add_argument('--peaks', help='Peak file')
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
            species_files = json.load(fp)
        # convert from unicode to byte here

    else:
        annotations = [args.dnase, args.blacklist, args.prom, args.enh, args.peaks]
        annotation_names = ['DNAse', 'Blacklist', 'Promoters', 'Enhancers', 'Peaks']

        species_files = {'chrsz': args.chromsizes,
                         'ref_fa': args.ref,
                         'annotations': annotations,
                         'annotation_names': annotation_names,
                         'tss_enrich': args.tss,
                         'reg2map': args.reg2map,
                         'roadmap_meta': args.meta
        }

    # set up inputs file
    if (args.files is not None):
        with open(args.files, 'r') as fp:
            data_files = json.load(fp)

        # convert from unicode to byte here
    else:
        peaks = [args.peaks, args.naive_overlap_peaks, args.idr_peaks]
        peak_names = ['Raw Peaks', 'Naive Overlap Peaks', 'IDR Peaks']

        data_files = {'genome': args.genome,
                     'sample_name': args.sample_name,
                     'final_reads_bed': args.finalbed,
                     'peaks': peaks,
                     'peak_names': peak_names,
                     'raw_bam': args.alignedbam,
                     'final_bam': args.finalbam,
                     'dup_log': args.duplog,
                     'pbc_log': args.pbc
        }


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


def viz_ataqc(qc_groups, sample_name, outprefix):
    """View results - html or multiQC"""

    # give option of native html

    viz.write_html(qc_groups, sample_name, outprefix)

    # also give option of multiQC
    

    return None


def main():

    args = parse_args()
    qc_groups = run_ataqc(args)
    viz_ataqc(qc_groups, args.sample_name, args.outprefix)



if __name__ == "__main__":
    main()
