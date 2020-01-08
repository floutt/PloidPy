# PloidPy

## Introduction
PloidPy is a program written in Python designed to infer ploidy from next-generation reads aligned to a haploid reference genome. The program makes use of a stastical model representing the distribution of specific nucleotide counts and selects the most probable ploidy on the basis of a minimum AIC.

## Installation
Installation of PloidPy is relatively simple and can be done easily using pip. The most recent version can be installed using the following command:
```
pip install git+git://github.com/floutt/PloidPy
```
### Dependencies
In order to run PloidPy, the following dependencies are required:
- `Python 3.6+`
- `NumPy`
- `SciPy`
- `Statsmodels`
- `matplotlib`
- `seaborn`

## Getting Started
Once installed, PloidPy can be run on an indexed BAM file. For this tutorial we will download some BAM files created using simulated data from the [ploidyNGS](https://github.com/diriano/ploidyNGS) repository.
```
wget https://github.com/diriano/ploidyNGS/raw/master/test_data/simulatedDiploidGenome/Ploidy2.bowtie2.sorted.bam
wget https://github.com/diriano/ploidyNGS/raw/master/test_data/simulatedTriploidGenome/Ploidy3.bowtie2.sorted.bam
wget https://github.com/diriano/ploidyNGS/raw/master/test_data/simulatedTetraploidGenome/Ploidy4.bowtie2.sorted.bam
```
After downloading these files, we have to index all of these files using samtools.
```
samtools index Ploidy2.bowtie2.sorted.bam
samtools index Ploidy3.bowtie2.sorted.bam
samtools index Ploidy4.bowtie2.sorted.bam
```
From here we can extract the minor allele count (MAC) and total read coverage from each of the BAM files using `PloidPy process_bam`. In the following example the Phred-like mapping threshold is set to 15 
```
PloidPy process_bam --bam Ploidy2.bowtie2.sorted.bam --out diploid.count --quality 15
PloidPy process_bam --bam Ploidy3.bowtie2.sorted.bam --out triploid.count --quality 15
PloidPy process_bam --bam Ploidy4.bowtie2.sorted.bam --out tetraploid.count --quality 15
```

In addition to providing raw data for 
