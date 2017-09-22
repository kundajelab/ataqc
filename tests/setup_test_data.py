# description: script to set up test data

import os

def setup_mouse_data():
    """Set up test data
    """
    # first download data
    tar_file = "test_mouse_reads-250k.tar.gz"
    os.system("wget http://mitra.stanford.edu/kundaje/dskim89/ataqc/{}".format(tar_file))
    os.system("tar -xzf {}".format(tar_file))
    os.system("rm {}".format(tar_file))

    # TODO download annotations as part of ataqc
    
    # then set up data sheet
    data_prefix = tar_file.split(".tar")[0]
    with open("{0}/{0}.datasheet.txt".format(data_prefix), 'w') as fp:
        headers = [
            "name",
            "genome",
            "fastq1",
            "fastq2",
            "raw_bam",
            "alignment_log",
            "dup_log",
            "pbc_log",
            "final_bam",
            "final_bed",
            "bigwig",
            "peaks"]

        work_dir = "{}/{}".format(os.getcwd(), data_prefix)
        files = [
            data_prefix,
            "mm9",
            "{}/fastq/mouse_250Kreads_rep1.1.fastq.gz".format(work_dir),
            "{}/fastq/mouse_250Kreads_rep1.2.fastq.gz".format(work_dir),
            "{}/align/rep1/mouse_250Kreads_rep1.1.trim.PE2SE.bam".format(work_dir),
            "{}/qc/rep1/mouse_250Kreads_rep1.1.trim.PE2SE.align.log".format(work_dir),
            "{}/qc/rep1/mouse_250Kreads_rep1.1.trim.PE2SE.dup.qc".format(work_dir),
            "{}/qc/rep1/mouse_250Kreads_rep1.1.trim.PE2SE.nodup.pbc.qc".format(work_dir),
            "{}/align/rep1/mouse_250Kreads_rep1.1.trim.PE2SE.nodup.bam".format(work_dir),
            "{}/align/rep1/mouse_250Kreads_rep1.1.trim.PE2SE.nodup.tn5.tagAlign.gz".format(work_dir),
            "{}/signal/macs2/rep1/mouse_250Kreads_rep1.1.trim.PE2SE.nodup.tn5.pf.pval.signal.bigwig".format(work_dir),
            ("overlap={0}/peak/macs2/overlap/mouse_250Kreads_rep1.1.trim.PE2SE.nodup.tn5.pf.pval0.1.500K.naive_overlap.filt.narrowPeak.gz;"
             "raw={0}/peak/macs2/rep1/mouse_250Kreads_rep1.1.trim.PE2SE.nodup.tn5.pf.narrowPeak.gz").format(work_dir)
        ]

        fp.write("{}\n".format("\t".join(headers)))
        fp.write("{}\n".format("\t".join(files)))

    return None

setup_mouse_data()
