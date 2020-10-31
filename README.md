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


### Filter SVIM results

This script was used to remove 'small' structural variants from the SVIM results.

```
python filter_SVIM_VCF.py
--in <INPUT_VCF>
--out <OUTPUT_VCF>
--minsize <MINIMAL_SV_SIZE>[1000]
```


`--in` specifies a VCF file produced by SVIM.

`--out` specifies a VCF file which will include the variants passing the specified length filter.

`--minsize` specifies the minimal size of structural variants to be included in the output file. Default: 1000 bp.



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

`--in` specifies a VCF file produced by SVIM. Analysed variants are 'DUP:TANDEM', 'INV', 'INS:NOVEL', 'DEL', and 'DUP:INT'.

`--ref` specifies a FASTA file which matches the supplied VCF file.

`--fig` specifies the figure output file. The file format is determined by the file name extension e.g. PDF, PNG, JPEG, SVG. Support for the file formats depends on the local system.

`--score` specifies a minimal score cutoff to filter out low quality variants.


### Genome-wide distribution of coverage (sequencing depth)

This script is based on previous studies ([10.3390/genes10090671](https://doi.org/10.3390/genes10090671), [10.1534/g3.119.400847](https://doi.org/10.1534/g3.119.400847)). The coverage of a read mapping is plotted along the pseudochromosomes. This script is customized for the Col-0 reference genome sequence TAIR9/TAIR10.

```
python cov_plot2.py
--in <FULL_PATH_TO_COVERAGE_FILE>
--out <FULL_PATH_TO_OUTPUT_FILE>
--ref <FULL_PATH_TO_REFERENCE_COVERAGE_FILE>
					
optional:
--name <NAME>
```

`--in` specifies a coverage file. This TAB-separated file contains the pseudochromosome name, position, and coverage of the respective position. This file can be generated based on a BAM file using https://github.com/bpucker/MGSE/blob/master/construct_cov_file.py.

`--ref` specifies a FASTA file which was used to generate the read mapping. It is important that pseudochromosome names match the names in the coverage file.

`--out` specifies an output folder, where the figure file and the data file will be placed. This folder will be created if it does not exist already.

`--name` specifies a name to be included in all file names.


### Analysis of long read per base quality

This script is intended to analyse the local quality along ONT reads.

```
python analyze_read_quality.py
--in <FASTQ_INPUT_FILE>
--info <INFO_FILE>
--out <OUTPUT_FOLDER>
```

`--in` specifies a FASTQ input file which contains long reads.

`--info` specifies a text file with details about the reads of interest. This file is TAB-separated with the following columns: SampleID, ReadID (need to match a read in the supplied FASTQ), Start, End, and Comment (optiona). Start and end can be used to analyse just a fraction of a read.

`--out` specifies the output folder.







### References

Pucker, B. (2019). Mapping-based genome size estimation. bioRxiv. doi:[10.1101/607390](https://doi.org/10.1101/607390).

Pucker B, Rückert C, Stracke R, Viehöver P, Kalinowski J, Weisshaar B. Twenty-Five Years of Propagation in Suspension Cell Culture Results in Substantial Alterations of the Arabidopsis Thaliana Genome. Genes. 2019. doi:[10.3390/genes10090671](https://doi.org/10.3390/genes10090671).

Busche M*, Pucker B*, Viehöver P, Weisshaar B, Stracke R. (2020). Genome Sequencing of Musa acuminata Dwarf Cavendish Reveals a Duplication of a Large Segment of Chromosome 2. G3: Genes|Genomes|Genetics. doi:[10.1534/g3.119.400847](https://doi.org/10.1534/g3.119.400847).

