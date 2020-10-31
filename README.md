# GKseq

These scripts were written for the analyis of Arabidopsis thaliana ONT sequencing data.


### Selection of ultra long ONT reads

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





### References

