# metrics that are run on multiple QC types

import pybedtools

def get_fract_reads_in_regions(reads_bed, regions_bed):
    """Function that takes in bed file of reads and bed file of regions and
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


def get_annotation_enrichments(data_files, species_files):
    """Run through annotation files (species file) and get enrichment of reads. Measure of signal to noise"""
    final_reads_bed = data_files['final_reads_bed']
    annotation_enrichments = {}
    for i in range(len(species_files['annotations'])):
        read_count, read_fract = get_fract_reads_in_regions(final_reads_bed, species_file['annotations'][i])
        annotation_enrichments[species_file['annotation_names'][i]] = (read_count, read_fract)
        
    return annotation_enrichments
    

def get_final_read_count(metrics):
    """Get final mapped reads compared to initial reads"""

    assert 'raw_bam' in metrics.keys()
    assert 'final_bam' in metrics.keys()

    num_raw_reads = metrics['raw_bam']['read_count']
    num_final_reads = metrics['final_bam']['read_count']

    final_read_fraction = float(num_final_reads) / num_raw_reads
    
    return (num_final_reads, final_read_fraction)


def run_metrics(data_files, all_metrics):
    """Run all QC"""
    
    metrics = {}
    metrics['final_reads'] = get_final_read_count(all_metrics)
    #metrics['annotation_enrichments'] = get_annotation_enrichments(data_files, all_metrics)

    return metrics


