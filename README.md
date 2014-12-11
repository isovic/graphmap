## GraphMap - A very accurate and sensitive long-read, high error-rate sequence mapper

Version: v0.18b
Release date: 11 December 2014

Precompiled binary, built on Mint 16 (Ubuntu 13.10) x64.
Tested on Mint 17 (Ubuntu 14.04), Ubuntu Server 14.04, Fedora 20 and Gentoo.

### Description
GraphMap is a novel mapper targeted at aligning long, error-prone third-generation sequencing data.
It can handle Oxford Nanopore data with unmatched accuracy, and also presents a significant improvement to the state-of-the-art for PacBio read mappers (namely, compared to BLASR and BWA-MEM).

Our mapper was also aimed at ease-of-use: the default parameters should be able to address any data you present it.
This is an extremely important feature for technologies with highly varying error rates, especially when also the error-profile continues changing.
Although GraphMap was not designed for NGS data, it can also handle reads from devices such as Roche and Illumina sequencers with very high accuracy, albeit slower than the best existing methods.

Please keep in mind that this is an early development version, and some parameters might have to be tweaked in order to best suit the rapidly changing Oxford Nanopore data.

### Comparison to other mappers
![comparison1](isovic.github.com/graphmap/doc/comparison/comparison1.png)



### Usage

```
# Process all reads from a given FASTA/FASTQ file with maximum number of threads:
./graphmap -r escherichia_coli.fa -d reads.fastq -o alignments.sam

# Limit the number of threads to 8, and load reads in batches of 50MB:
./graphmap -t 8 -B 50 -r escherichia_coli.fa -d reads.fastq -o alignments.sam

# Process only the first 1000 reads:
./graphmap -B 0 -n 1000 -r escherichia_coli.fa -d reads.fastq -o alignments.sam

# Process all reads from a given folder.
./graphmap -r escherichia_coli.fa -D reads_folder -O alignments_folder

# Generate only the index.
./graphmap -I -r escherichia_coli.fa
```



### Contact information
For additional information, help and bug reports please send an email to one of the following:
ivan.sovic@irb.hr, mile.sikic@fer.hr, nagarajann@gis.a-star.edu.sg

