ATAqC
===================================================


This pipeline is designed for collecting advanced quality control metrics for ATAC-seq datasets.


===================================================

Annotation sets:
hg19 - Gencode v19
hg38 - Gencode v24
mm9 - vM1 (though I believe the ENCODE portal no longer supports mm9)
mm10 - vM7

The TSS bed files are generated directly from the Gencode full GTF files, with the following command:

```
zcat $GTF | 
grep -P '\tgene\t' | 
grep 'protein_coding' | 
grep -v 'level 3' | 
awk -F '[\t|\"]' '{ print $1"\t"$4"\t"$5"\t"$10"\t0\t"$7 }' | 
awk -F '\t' 'BEGIN{ OFS="\t" } { if ($6=="+") { $3=$2-1; $2=$2-2 } else { $2=$3; $3=$3+1 } print }' |
sort -k1,1 -k2,2n > $TSS
```

*Note that the TSS file is a point file, and is not the same as the promoter file (described below).

===================================================

The promoter/enhancer annotations are a little trickier (and likely should be updated, given that the annotations are based off the data that was in the ENCODE portal as of 03/27/2016):

These annotations should be viewed as preliminary and approximate, not as part of the ENCODE encyclopedia. Good for QC, but for deeper analysis please do consider carefully the process by which we got these annotations.

hg19 - we made use of the high stringency (-log_pval > 10) Reg2Map promoter and enhancer sets: https://personal.broadinstitute.org/meuleman/reg2map/HoneyBadger2_release/
hg38 - the promoter set is the union of ENCODE RAMPAGE peaks, the enhancer set is the remainder of the union of ENCODE open chromatin (ie DNase) peak sets after removing blacklist and promoter regions. Ie, any other site that was accessible in some ENCODE DNase experiment that was not labeled as a promoter by RAMPAGE data. There is no Reg2Map resource for hg38.

mm9 - after getting the union of all mouse mm9 DNase peaks available in the ENCODE portal, the promoter set is those peaks that overlap the TSS file, and the enhancer set is the rest.
mm10 - the promoters are the predicted promoters from the ENCODE portal (https://www.encodeproject.org/data/annotations/v3/). after getting the union of all mouse mm10 DNase peaks available in the ENCODE portal,  the enhancer set is the remainder after subtracting the promoter and blacklist, since one exists for mm10.

Whenever a blacklist is mentioned, it's the recorded file from the ENCODE portal.

===================================================

# Known issues
The pipeline is not currently compatible with samtools/1.3 - we are working on this incompatibility


# Contributors

* Daniel Kim - MD/PhD Student, Biomedical Informatics Program, Stanford University
* Chuan Sheng Foo - PhD Student, Computer Science Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
