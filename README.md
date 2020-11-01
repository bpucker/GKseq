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


### Reduce sequence names in assembly file to contig names

This script is intended to remove everything except the actual contig name from the header line in a multiple FASTA file.

```
python reduce_contig_names.py
--in <INPUT_FILE>
--out <OUTPUT_FILE>
```

`--in` specifies a FASTA input file which contains sequences with complex header lines.

`--out` specifies the output folder. Sequences will be written to the output file, but header lines will be reduced to the actual contig name.



### Clean assembly

This script is intended to remove contamination from a long read assembly file. The general concept is based on previous work ( [10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321), [10.3389/fmolb.2018.00062](https://doi.org/10.3389/fmolb.2018.00062), [10.3390/genes11030274](https://doi.org/10.3390/genes11030274)).

```
python clean_assembly.py
--out <OUTPUT_FOLDER>
--assembly <ASSEMBLY_FILE>
--cov <COVERAGE_FILE>
--organel <ORGANEL_SEQ_FILE>
--white <WHITE_LIST_SEQ_FILE>
--black <BLACK_LIST_SEQ_FILIE>
					
optional:
--minlen <INT, MINIMAL_CONTIG_LENGTH>
--mincov <INT, MINIMAL_READ_MAPPING_COVERAGE_FOR_FILTERING>
--maxcov <INT, MAXIMAL_READ_MAPPING_COVERAGE_FOR_FILTERING>
--wordsize <INT, WORD_SIZE_FOR_BLAST>
--cpus <INT, NUMBER_OF_THREADS_FOR_BLAST_SEARCH>
```

`--assembly` specifies a FASTA input file which contains sequences with complex header lines.

`--out` specifies the output folder. All output files will be stored in this folder. The folder will be created if it does not exist already.

`--cov` specifies a coverage file. The sequence names in this coverage file must match the sequence names in the assembly file. The average coverage per contig can be used as filter.

`--organel` specifies a FASTA file containing the organellar genome sequences of the species of interest or closely related species. This is used to identiy contigs belonging to these subgenomes.

`--white` specifies a FASTA file with white list sequences. Hits against these sequences will be used to preserve the sequences during the filter steps.

`--black` specifies a FASTA file with black list sequences. Hits against these sequences will be used to discard the sequences during the filter steps unless their are also hits against the white list sequences.

`--minlen` specifies a minimal contig length. All contigs shorter than this value will be removed in the filtering. This value should be set based on the read length distribution. Default: 100000 bp (100kb).

`--mincov` specifies a minimal coverage in the read mapping that is expected for all valid contigs. This filter can be used to identify artifacts with low coverage. This filter needs to be adapted to the average coverage. Default: 10.

`--maxcov` specifies a maximal coverage in the read mapping that is expected for all valid contigs. This filter can be used to identify and discard collapsed sequences or organellar sequences. This value should be set depending on the average coverage. Default: 200.

`--wordsize` specifies the '-word_size' parameter which is passed on to BLAST. Using a large value can substantially speed up the process, but decreases the sensitivity. Default: 30.

`--cpus` specifies the number of threads to use for BLAST. Default: 20.





### References

Pucker, B (2019). Mapping-based genome size estimation. bioRxiv. doi:[10.1101/607390](https://doi.org/10.1101/607390).

Pucker B, Rückert C, Stracke R, Viehöver P, Kalinowski J, Weisshaar B. Twenty-Five Years of Propagation in Suspension Cell Culture Results in Substantial Alterations of the Arabidopsis Thaliana Genome. Genes. 2019. doi:[10.3390/genes10090671](https://doi.org/10.3390/genes10090671).

Busche M*, Pucker B*, Viehöver P, Weisshaar B, Stracke R. (2020). Genome Sequencing of Musa acuminata Dwarf Cavendish Reveals a Duplication of a Large Segment of Chromosome 2. G3: Genes|Genomes|Genetics. doi:[10.1534/g3.119.400847](https://doi.org/10.1534/g3.119.400847).

Pucker, B, Holtgräwe, D, Rosleff Sörensen, T, Stracke, R, Viehöver, P, and Weisshaar, B (2016). A de novo Genome Sequence Assembly of the Arabidopsis thaliana Accession Niederzenz-1 Displays Presence/Absence Variation and Strong Synteny. PloS-ONE 11:e0164321. doi:[10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321).

Siadjeu C*, Pucker B*, Viehoever P, Albach D and Weisshaar B (2020). High contiguity de novo genome sequence assembly of Trifoliate yam (Dioscorea dumetorum) using long read sequencing. Genes. doi:[10.3390/genes11030274](https://doi.org/10.3390/genes11030274).

Haak, M, Vinke, S, Keller, W, Droste, J, Rückert, C, Kalinowski, J, & Pucker, B. (2018). High Quality de novo Transcriptome Assembly of Croton tiglium. Frontiers in Molecular Biosciences, 5. doi:[10.3389/fmolb.2018.00062](https://doi.org/10.3389/fmolb.2018.00062).
