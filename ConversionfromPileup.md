# Introduction #

In CONDEX, we use only per exon measurements. The pileup format contains the read depth and calls on a per base pair format. This method converts from a per base pair to a per exon format. The output format of the file is tab delimited and contains the following fields:

```
chr 
exon_start_position 
exon_end_position 
gene_description 
number_of_heterozygous_sites 
number_of_reads_in_exon 
state
```

Initially, the state is the normal state for each exon.
The final column is the ground truth state (only relevant for simulations)

# Options #
  * `-e <exon file>` : prefix name of file(s) with the list of exons in the chromosome. Each file should contain values from a single chromosome.
  * `-pileup <pileup file>` : name of pileup file
  * `-c <start> - <end>` : list of chromosomes being processed (not used)
  * `-t` : print timing metrics
  * `-o <output file name>` : name of file to output the exon-level measurements. Default is StdOut

# Classes used #
ConvertToExonFormat
Pileup
Exon


# Description #
ConvertToExonFormat contain the main function. It reads in the list of exons as described by the capture array file (which contains the list of exons captures by that array). It then parses the pileup file and adds each base pair to the appropriate exon and updates the measurements for the exon accordingly.

# Example #
  1. Download and unzip "captureExons.zip"
  1. If input is in pileup format:
```
java -Xmx100M -Xms100m -jar PileupToExon.jar 
-e "captureExons" 
-c "10-10" 
-o "sample.exon" 
-t 
-pileup "sample.chr10.pileup"
```
  1. If input is in BAM format:
```
samtools pileup -cf <REF> <BAM> 
| awk -f pileupFilter.awk 
| java -Xmx100M -Xms100m -jar PileupToExon.jar 
   -e "captureExons" 
   -c "1-22" 
   -o <OUT> 
   -t
```
where
    * OUT is the name of the output file
    * BAM is the name of the bam file to be converted
    * REF is the human reference is fasta format
pileupFilter.awk can be found in the downloads section and filters calls based on snp quality and caps the read depth at 1000.