# Introduction #
CONDEX reads in the sample to be called, the other samples used as a baseline, and finds the most probable state for each exon. The output format is tab delimited and contains the following fields:
```
chr 
exon_start_position 
exon_end_position 
gene_description 
number_of_heterozygous_sites 
number_of_reads_in_exon 
state
chi-square statistic
pValue
average_number_of_heterozygous_sites
average_number_of_reads_in_exon 
```
The final column is the ground truth state (only relevant for simulations)

# Options #
  * `-e <exon file>`: name of file with the list of exons in the chromosome. Each file should contain values from a single chromosome.
  * `-o <output file>` : name of output file. Default is StdOut.
  * `-c <start>-<end>` : start and ending chromosomes analyzed (not used)
  * `-b <baseline files>` : comma separated list of exon files to be used as baselines (no spaces separating the list)
  * `-p <parameter file>` : name of file with the values of parameters. [example](http://code.google.com/p/condr/downloads/detail?name=SimulationParameterFile.04102011.length200000.rate500000&can=2&q=). "//" can be used for comments within the file
  * `-threshold <value>` : threshold above which to call as CNV (allows more calling of normal state)
  * `-t` : prints timing metrics

# Classes Used #
CONDR
Exon
HiddenMarkovModel
State
Probability

# Description #
The rough outline of steps followed are:
  * read measurements from samples acting as baseline
  * calculate the average measurements
  * read measurements from sample to be called
  * using forward-backward algorithm calculate the most likely state
  * output the measurements with the most likely state

# States #
For now, the state encoding is:
```
0 ~ Normal state (copy number = 2)
1 ~ Homozygous deletion state (copy number = 0)
2 ~ Heterozygous deletion state (copy number = 1)
3 ~ Amplification state (copy number = 3)
4 ~ Amplification state (copy number = 4)
```

# Example #
  1. Download and unzip [reference samples](http://code.google.com/p/condr/downloads/detail?name=referenceSamples.zip&can=2&q=).
  1. Run CONDR
```
java -jar -Xmx500M -Xms500m CONDEX.jar 
-c "10-10" -e "sample.exon.chr10" 
-t 
-b "sample1.chr10.exon,sample10.chr10.exon,sample11.chr10.exon,
sample12.chr10.exon,sample13.chr10.exon,sample14.chr10.exon,
sample15.chr10.exon,sample16.chr10.exon,sample2.chr10.exon,
sample3.chr10.exon,sample4.chr10.exon,sample5.chr10.exon,
sample6.chr10.exon,sample7.chr10.exon,sample8.chr10.exon,sample9.chr10.exon" 
-p "SimulationParameterFile.04102011.length200000.rate500000" 
-o "sample.chr10.result" 
-threshold 0
```