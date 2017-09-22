#!/usr/bin/env python

import subprocess
import pysam
import gzip


def get_fract_reads_in_regions(reads_bed, regions_bed, outprefix):

    # make Popen object for sorting bed
    sort_process = subprocess.Popen(('bedtools', 'sort', '-i', regions_bed), 
                                    stdout=subprocess.PIPE)

    # make Popen object for merging bed
    merge_process = subprocess.Popen(('bedtools', 'merge', '-i', 'stdin'), 
                                     stdin=sort_process.stdout, 
                                     stdout=subprocess.PIPE)

    # write intersection of beds to file to correctly format
    # (default output has many newline and tab characters)
    text_file = open("{0}_intersect.txt".format(outprefix), "w")
    intersect_process = subprocess.call(('bedtools', 'intersect', '-c', '-nonamecheck', '-a', 'stdin', '-b', reads_bed), 
                                        stdin=merge_process.stdout, 
                                        stdout=text_file)
    text_file.close()

    # get read count for intersection
    intersection_read_count = int(subprocess.check_output(('wc', '-l', '{0}_intersect.txt'.format(outprefix))).split()[0])

    # get read count for final_bed
    zcat_process = subprocess.Popen(('zcat', reads_bed), stdout=subprocess.PIPE)
    reads_bed_count = int(subprocess.check_output(('wc', '-l'), stdin=zcat_process.stdout).split('\n')[0])

    fract_reads = float(intersection_read_count)/reads_bed_count

    return intersection_read_count, fract_reads


    
    

def get_fract_reads_in_regions_old(reads_bed, regions_bed):
    reads_bedtool = pybedtools.BedTool(reads_bed)
    regions_bedtool = pybedtools.BedTool(regions_bed)

    reads = regions_bedtool.sort().merge().intersect(reads_bedtool, c=True, nonamecheck=True)

    read_count = 0
    for interval in reads:
        read_count += int(interval[-1])
    fract_reads = float(read_count)/reads_bedtool.count()

    return read_count, fract_reads


def get_annotation_enrichments(data_files, species_files, outprefix):
    final_reads_bed = data_files['final_bed']# Read 'final_reads_bed' field of data_files into variable
    annotation_enrichments = {}
    for annotation_key in species_files["annotations"].keys():
        read_count, read_fract = get_fract_reads_in_regions(final_reads_bed, species_files['annotations'][annotation_key], outprefix)
        annotation_enrichments[annotation_key] = (read_count, read_fract)

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


def get_signal_to_noise(final_bed, dnase_regions, blacklist_regions,
                        prom_regions, enh_regions, peaks):
    # Dnase regions
    reads_dnase, fract_dnase = get_fract_reads_in_regions(final_bed, dnase_regions)
    
    # Blacklist regions
    reads_blacklist, fract_blacklist = get_fract_reads_in_regions(final_bed, blacklist_regions)
    
    # Prom regions
    reads_prom, fract_prom = get_fract_reads_in_regions(final_bed, prom_regions)
    
    # Enh regions
    reads_enh, fract_enh = get_fract_reads_in_regions(final_bed, enh_regions)
    
    # Peak regions
    reads_peaks, fract_peaks = get_fract_reads_in_regions(final_bed, peaks)
    
    return reads_dnase, fract_dnase, reads_blacklist, fract_blacklist, reads_prom, fract_prom, reads_enh, fract_enh, reads_peaks, fract_peaks

def test_integrative():
    print 'Integrative'
