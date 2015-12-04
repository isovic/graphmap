## GraphMap - A highly sensitive and accurate mapper for long, error-prone reads  
**__Current Version: 0.22__**  
Release date: 12 November 2015  
  
Updates:  
- Many tiny bug fixes, mostly related to anchored alignment. It should be slightly more sensitive now.  
- Two overlap modes merged from the dev branch: ```-w owler``` (fast, uses a trimmed GraphMap pipeline, reports output in the MHAP format) and ```-w overlapper``` (full GraphMap pipeline including alignment, output in SAM format). For usage - check examples at the bottom.  
- GraphMap integration into marginAlign - we forked marginAlign and extended it to support GraphMap alongside to LAST and BWA-MEM ([https://github.com/isovic/marginAlign](https://github.com/isovic/marginAlign)). Use parameters ```--graphmap``` or ```--graphmapanchor``` with marginAlign to specify the mapper.  

For more information on overlapping, take a look at [overlap.md](overlap.md).  

GraphMap is also used as an overlapper in a new *de novo* genome assembly project called [Ra](https://github.com/mariokostelac/ra-integrate) [https://github.com/mariokostelac/ra-integrate](https://github.com/mariokostelac/ra-integrate).  
Ra attempts to create *de novo* assemblies from raw nanopore and PacBio reads without requiring error correction, for which a highly sensitive overlapper is required.  


### Quick start on Linux x64
```  
git clone https://github.com/isovic/graphmap.git  
cd graphmap  
make  
  
# To align:  
./bin/Linux-x64/graphmap -r reference.fa -d reads.fasta -o output.sam  

# To overlap:  
./bin/Linux-x64/graphmap -w owler -r reads.fasta -d reads.fasta -o output.mhap  
```  

### Description
GraphMap is a novel mapper targeted at aligning long, error-prone third-generation sequencing data.  
It is **designed to handle Oxford Nanopore MinION 1d and 2d reads** with very high sensitivity and accuracy, and also presents a significant improvement over the state-of-the-art for PacBio read mappers.

GraphMap was also designed for ease-of-use: the **default parameters** can handle a wide range of read lengths and error profiles, including: *Illumina*, *PacBio* and *Oxford Nanopore*.  
This is an especially important feature for technologies where the error rates and error profiles can vary widely across, or even within, sequencing runs.  

**The GraphMap algorithm** is structured to achieve high-sensitivity and speed using a five-stage
read-funneling approach. In stage I, GraphMap uses a novel adaptation of gapped spaced seeds to efficiently reduce the search space and get seed hits as a form of coarse alignment. These are then refined in stage II using graph-based vertex-centric processing of seeds to efficiently construct alignment anchors. GraphMap then chains anchors using a kmer
version of longest common subsequence (LCS) construction (stage III), refines
alignments with a form of L1 linear regression (stage IV) and finally evaluates the
remaining candidates to select the best location to reconstruct a final alignment (stage V).
GraphMap computes a BLAST-like E-value as well as a mapping quality for its alignments.

**Evaluation** on MinION sequencing datasets against short and long-read mappers indicates that GraphMap increases mapping sensitivity by at least 15-80%. GraphMap alignments are the first to demonstrate consensus calling with <1 error in 100,000 bases, variant calling on the human genome with 76% improvement in sensitivity over the next best mapper (BWA-MEM), precise detection of structural variants from 100bp to 4kbp in length and species and strain-specific identification of pathogens using MinION reads.

Further details about the algorithm, comparison with other mappers and usage applications can be found in the **preprint** of our paper:  
[Fast and sensitive mapping of error-prone nanopore sequencing reads with GraphMap](http://biorxiv.org/content/early/2015/06/10/020719)  

**Nanopore sequencing data** of E. Coli UTI89 generated in-house and used in the paper now available on ENA:  
[PRJEB9557](http://www.ebi.ac.uk/ena/data/view/PRJEB9557)  

### Features  
- Mapping position agnostic to alignment parameters.
- Consistently very high sensitivity and precision across different error profiles, rates and sequencing technologies even with default parameters.
- Circular genome handling to resolve coverage drops near ends of the genome.
- E-value.
- Meaningful mapping quality.
- Various alignment strategies (semiglobal bit-vector and Gotoh, anchored).

### Installation
To build GraphMap from source type:  
```
make
```  
Required libraries are prebuilt for Linux x64 systems.  
To rebuild them for other systems, type:  
```
make deps
```  

You will need a recent GCC/G++ version (>=4.7).

More installation instructions can be found in the INSTALL file.


### Usage examples
```
# Align all reads from a given FASTA/FASTQ file with default number of threads using semiglobal bit-vector alignment:  
./graphmap -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Align all reads from a given FASTA/FASTQ file using anchored alignment approach:  
./graphmap -a anchor -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Overlap all reads from a given FASTA/FASTQ file and report overlaps in MHAP format (fast):
./graphmap -w owler -r reads.fa -d reads.fa -o overlaps.mhap  

# Overlap all reads from a given FASTA/FASTQ in a full GraphMap mode with generating alignments (slow):
./graphmap -w overlapper -r reads.fa -d reads.fa -o overlaps.sam  

# Align reads using the Gotoh for semiglobal alignment:  
./graphmap -a gotoh -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Align reads using Gotoh alignment with anchored approach: 
./graphmap -a anchorgotoh -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Process reads from a circular genome:  
./graphmap -C -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Threshold the E-value of alignments to 1e-100. Alignments with E-value > 1e-100 will be called unmapped:  
./graphmap -z 1e-100 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Output all similarly good alignments (to within F*num_kmers_of_best_alnmnt) instead of only one best:  
./graphmap -Z -F 0.05 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Limit the number of threads to 8, and load reads in batches of 50MB:  
./graphmap -t 8 -B 50 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Align reads using more sensitive parameters for Illumina data (currently equivalent to "-a gotoh"):  
./graphmap -x illumina -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Load all reads in one batch and align only the first 1000 reads:  
./graphmap -B 0 -n 1000 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Process all reads from a given folder.  
./graphmap -r escherichia_coli.fa -D reads_folder -O alignments_folder  

# Generate only the index.  
./graphmap -I -r escherichia_coli.fa  

# Run a debug version of GraphMap (build with "make debug") and verbose the SAM output to see various info about alignment:  
./graphmap-debug -b 3 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

```  

### Contact information

For additional information, help and bug reports please send an email to one of the following:  
ivan.sovic@irb.hr, mile.sikic@fer.hr, nagarajann@gis.a-star.edu.sg

### Acknowledgement  
This work was supported by the IMaGIN platform (project No. 102 101 0025), through a grant from the Science and Engineering Research Council, funding to the Genome Institute of Singapore from the Agency for Science, Technology and Research (A*STAR), Singapore, and funding from the Croatian Science Foundation (Project no. UIP-11-2013-7353 - Algorithms for Genome Sequence Analysis).  
