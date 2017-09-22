#!/usr/bin/env python

import subprocess
import pysam
import gzip

from collections import OrderedDict


def get_fract_reads_in_regions(reads_bed, regions_bed):
    """Function that takes in bed file of reads and bed file of regions and
    gets fraction of reads sitting in said regions
    """
    # sort the bed, then merge, then write to correctly format
    sort_process = subprocess.Popen(
        ('bedtools', 'sort', '-i', regions_bed),
        stdout=subprocess.PIPE)
    
    # make Popen object for merging bed
    merge_process = subprocess.Popen(
        ('bedtools', 'merge', '-i', 'stdin'),
        stdin=sort_process.stdout,
        stdout=subprocess.PIPE)

    # write intersection of beds to file to correctly format
    # (default output has many newline and tab characters)
    intersect_process = subprocess.Popen(
        ('bedtools', 'intersect', '-c', '-nonamecheck', '-a', 'stdin', '-b', reads_bed),
        stdin=merge_process.stdout,
        stdout=subprocess.PIPE)

    out, err = intersect_process.communicate()
    out = out.splitlines()
    
    read_count = 0
    for line in out:
        read_count += int(line.strip().split('\t')[3])
        intersection_read_count = read_count
    
    # get read count for final_bed
    zcat_process = subprocess.Popen(
        ('zcat', reads_bed), stdout=subprocess.PIPE)
    reads_bed_count = int(subprocess.check_output(
        ('wc', '-l'), stdin=zcat_process.stdout).split('\n')[0])

    fract_reads = float(intersection_read_count) / reads_bed_count

    return intersection_read_count, fract_reads
    

def get_annotation_enrichments(data_files, species_files, outprefix):
    final_reads_bed = data_files['final_bed']# Read 'final_reads_bed' field of data_files into variable
    annotation_enrichments = OrderedDict()
    for annotation_key, annotation_bed in species_files["annotations"].items():
        if annotation_key == "tss":
            continue
        read_count, read_fract = get_fract_reads_in_regions(final_reads_bed, annotation_bed)
        annotation_enrichments[annotation_key] = (read_count, read_fract)
    # Also add called peaks
    read_count, read_fract = get_fract_reads_in_regions(final_reads_bed, data_files['peaks'][data_files["peaks"].keys()[0]])
    annotation_enrichments["called_peaks"] = (read_count, read_fract)
        
    return annotation_enrichments


def get_final_read_count(metrics):
    assert 'raw_bam' in metrics.keys() # Check if 'raw_bam' key is in metrics dictionary
    assert 'final_bam' in metrics.keys() # Check if 'final_bam' key is in metrics dictionary

    num_raw_reads = metrics['raw_bam']['read_count'] # Read the read counts from the initial bam into variable
    num_final_reads = metrics['final_bam']['read_count'] # Read the read counts from the final bam into variable

    final_read_fraction = float(num_final_reads) / num_raw_reads # Find the number of reads that weren't filtered out

    return (num_final_reads, final_read_fraction)


def run_metrics(all_metrics, data_files, species_files, outprefix, encode_only=False):
    metrics = {}
    metrics['final_reads'] = get_final_read_count(all_metrics)
    metrics['annotation_enrichments'] = get_annotation_enrichments(data_files, species_files, outprefix)
    
    return metrics                

