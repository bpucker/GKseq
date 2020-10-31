# GKseq

These scripts were written for the analyis of Arabidopsis thaliana ONT sequencing data. All scripts are intended for the use with Python v2.7.


### Selection of ultra long ONT reads

This script was used to extract the longest reads of various sequencing datasets for a combined Col-0 genome assembly.

```
python select_long_reads.py
--in <FULL_PATH_TO_INPUT_FILE>
--out <FULL_PATH_TO_OUTPUT_FILE>

optional:
--cut <LENGTH_CUTOFF_IN_BP>
```          

`--in` specifies a gzip compressed FASTQ input file which contains the long reads.

`--out` specifies a gzip compressed FASTQ output file. The selected ultra long reads are written into this file.

`--cut` specifies a length cutoff for the selection of ultra long reads. Only reads passing this filter will be written into the output file. Default: 3001bp.


### Genome-wide distribution of sequence variants

This script is intended to generate a genome-wide distribution plot of SVIM results.

```
python genome_wide_distribution.py
--in <SVIM_VCF_FILE>
--ref <REFERENCE_SEQUENCE>
--fig <FIGURE_NAME_WITH_PDF_EXTENSION>

optional:
--score <INT, minimal SVIM score to consider variant>[10]
```

`--in` specifies a VCF file produced by SVIM.

`--ref` specifies a FASTA file which matches the supplied VCF file.

`--fig` specifies the figure output file. The file format is determined by the file name extension e.g. PDF, PNG, JPEG, SVG. Support for the file formats depends on the local system.

`--score` specifies a minimal score cutoff to filter out low quality variants.






### References

