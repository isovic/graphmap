## GraphMap - A highly sensitive and accurate mapper for long, error-prone reads  
**__Current Version: 0.5.0__**  
Release date: 28 February 2017  

**\*new\* - Minimizer index**  
The hash index has now been completely reimplemented.  
Index construction is now much faster than before.  

There are several new important features to mention:  
- The index supports minimizers - the minimizer window size can be specified via command line and by default is equal to ```5```. To switch off minimizers, use ```--minimizer-window 1```.  
- The number of hits for each seed lookup is now thresholded by a percentile value, e.g. ```--freq-percentile 0.99``` means that ```1%``` of the most repetitive seeds will be skipped. To turn off this feature, specify ```--freq-percentile 1.0```. This feature slightly reduces sensitivity, but plays a huge role when mapping to large references.  
- ```--fly-index``` enables building the index on the fly, instead of storing it to disk.  
- To switch off minimizers and percentile filtering, there is a composite option ```-x sensitive```. This will produce results similar to previous versions.
- In case the version of an index file is not compatible with GraphMap, previous versions would automatically overwrite the index file. This is now prevented by default, and GraphMap now simply halts with a command line message. To force automatic rebuild if necessary (but not rebuild if a valid index file exists), use: ```--auto-rebuild-index```. To rebuild the index in any case, specify: ```--rebuild-index```.  
- Renamed the ```--sensitive``` mode to ```--double-index```  

Owler overlapping mode was also enhanced with the new index, and now works much faster on larger datasets.  

Many bug fixes were also made, including the ones related to mapping to circular references, various segfaults, and extra new lines when outputting secondary alignments. Some of the segfaults were caused by the index, which is now addressed with the new version.

Also, this release fixes an issue with **transcriptome mapping**, where recall would drop (however, not precision).  

**Important** - When building an index, memory consumption is now larger than before. However, once minimizers have been collected and the index stored, the final size of the index is much smaller. Concrete estimate is: ```32 x reference_size``` for index construction, and ```~7 x reference_size``` for the final index. This means that for a human genome, constructing the index would peak at about ```~102 GB``` while the final index would have only ```~20 GB```.  

The large peak **can and will be addressed** in the next minor release. This should not require the updating of the index versions, and the index should be compatible to ```v0.5.0```.  


**Mapping to transcriptomes**  
GraphMap can now accept a GTF file to internally construct a transcriptome sequence from a given reference genome, and map RNA-seq data to it. The final alignments are converted back to genome space by placing ```N``` operations in the CIGAR strings.  

To use the new transcriptome mapping option simply specify a GTF file using the ```--gtf``` option:  
```graphmap align -r reference.fa --gtf reference.gtf -d reads.fastq -o out.sam```  


For a detailed change log from the previous release, take a look at [doc/changelog.md](doc/changelog.md).  

For more information on overlapping, take a look at [overlap.md](overlap.md).  
For detailed installation instructions, take a look at [INSTALL.md](INSTALL.md) file.  
Description of custom parameters in GraphMap's SAM output can be found at [doc/sam_output.md](doc/sam_output.md).  

### Features  
- Mapping position agnostic to alignment parameters.  
- Consistently very high sensitivity and precision across different error profiles, rates and sequencing technologies even with default parameters.  
- Circular genome handling to resolve coverage drops near ends of the genome.  
- E-value.  
- Meaningful mapping quality.
- Various alignment strategies (semiglobal bit-vector and Gotoh, anchored).  
- **Overlapping** of reads for *de novo* assembly.  
- **Transcriptome mapping** through internal construction of a transcriptome from a given genomic reference and a GTF file.  
- ...and much more.  

GraphMap is also used as an overlapper in a new *de novo* genome assembly project called [Ra](https://github.com/mariokostelac/ra-integrate) ([https://github.com/mariokostelac/ra-integrate](https://github.com/mariokostelac/ra-integrate)).  
Ra attempts to create *de novo* assemblies from raw nanopore and PacBio reads without requiring error correction, for which a highly sensitive overlapper is required.  


### Quick start on Linux x64
```  
git clone https://github.com/isovic/graphmap.git  
cd graphmap  
make modules  
make  

# To align:  
./bin/Linux-x64/graphmap align -r reference.fa -d reads.fasta -o output.sam  

# To overlap:  
./bin/Linux-x64/graphmap owler -r reads.fasta -d reads.fasta -o output.mhap  
```  

### Description
GraphMap is a novel mapper targeted at aligning long, error-prone third-generation sequencing data.  
It is **designed to handle Oxford Nanopore MinION 1d and 2d reads** with very high sensitivity and accuracy, and also presents a significant improvement over the state-of-the-art for PacBio read mappers.

GraphMap was also designed for ease-of-use: the **default parameters** can handle a wide range of read lengths and error profiles, including: *Illumina*, *PacBio* and *Oxford Nanopore*.  
This is an especially important feature for technologies where the error rates and error profiles can vary widely across, or even within, sequencing runs.  

**The GraphMap algorithm** is structured to achieve high-sensitivity and speed using a five-stage
read-funneling approach. In stage I, GraphMap uses a novel adaptation of gapped spaced seeds to efficiently reduce the search space and get seed hits as a form of coarse alignment. These are then refined in stage II using graph-based vertex-centric processing of seeds to efficiently construct alignment anchors. GraphMap then chains anchors using a kmer
version of longest common subsequence (LCSk) construction (stage III), refines
alignments by chaining anchors in the anchored mode or with a form of L1 linear regression in the semiglobal alignment mode (stage IV) and finally evaluates the
remaining candidates to select the best location to reconstruct a final alignment (stage V).
GraphMap computes a BLAST-like E-value as well as a mapping quality for its alignments.

**Evaluation** on MinION sequencing datasets against short and long-read mappers indicates that GraphMap increases mapping sensitivity by at least 15-80%. GraphMap alignments are the first to demonstrate consensus calling with <1 error in 100,000 bases, variant calling on the human genome with 76% improvement in sensitivity over the next best mapper (BWA-MEM), precise detection of structural variants from 100bp to 4kbp in length and species and strain-specific identification of pathogens using MinION reads.

Further details about the algorithm, comparison with other mappers and usage applications can be found in the **preprint** of our paper:  
[Fast and sensitive mapping of error-prone nanopore sequencing reads with GraphMap](http://biorxiv.org/content/early/2015/06/10/020719)  

**Nanopore sequencing data** of E. Coli UTI89 generated in-house and used in the paper now available on ENA:  
[PRJEB9557](http://www.ebi.ac.uk/ena/data/view/PRJEB9557)  

### Installation
To build GraphMap from source type:  
```
make modules	# This pulls the latest version of all required submodules
make
```  

You will need a recent GCC/G++ version (>=4.7).  

Run ```sudo make install``` to install the graphmap binary to ```/usr/bin```.  

More installation instructions can be found in the [INSTALL.md](INSTALL.md) file.


### Usage examples
```
# **Align** all reads from a given FASTA/FASTQ file using anchored alignment approach:  
./graphmap align -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# **Overlap** all reads from a given FASTA/FASTQ file and report overlaps in MHAP format (fast):  
./graphmap owler -r reads.fa -d reads.fa -o overlaps.mhap  

# **Align** all reads to a transcriptome sequence:  
./graphmap align -r scerevisiae.fa --gtf scerevisiae.gtf -d reads.fastq -o alignments.sam  


# Align all reads and report alignments using the extended CIGAR format.  
./graphmap align -r escherichia_coli.fa -d reads.fastq -o alignments.sam --extcigar  

# Align all reads from a given FASTA/FASTQ file with default number of threads using semiglobal bit-vector alignment:  
./graphmap align -a sg -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Overlap all reads from a given FASTA/FASTQ in a full GraphMap mode with generating alignments (slow):  
./graphmap align -x overlap -r reads.fa -d reads.fa -o overlaps.sam  

# Align reads using the Gotoh for semiglobal alignment:  
./graphmap align -a sggotoh -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Align reads using Gotoh alignment with anchored approach:  
./graphmap align -a anchorgotoh -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Process reads from a circular genome:  
./graphmap align -C -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Threshold the E-value of alignments to 1e-100. Alignments with E-value > 1e-100 will be called unmapped:  
./graphmap align --evalue 1e-100 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Output all secondary alignments instead of only one best:  
./graphmap align --secondary -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Control the similarity for secondary alignments. All alignments to within F*num_covered_bases from the best will be output.  
./graphmap align --secondary -F 0.05 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Limit the number of threads to 8, and load reads in batches of 50MB:  
./graphmap align -t 8 -B 50 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Align reads using more sensitive parameters for Illumina data:  
./graphmap align -x illumina -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Load all reads in one batch and align only the first 1000 reads:  
./graphmap align -B 0 -n 1000 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Rebuild the index if it already exists:  
./graphmap align --rebuild-index -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

# Generate only the index.  
./graphmap align -I -r escherichia_coli.fa  

# Run a debug version of GraphMap (build with "make debug") and verbose the SAM output to see various info about alignment:  
./graphmap-debug align -b 3 -r escherichia_coli.fa -d reads.fastq -o alignments.sam  

```  

### Contact information

For additional information, help and bug reports please send an email to one of the following:  
ivan.sovic@irb.hr, mile.sikic@fer.hr, nagarajann@gis.a-star.edu.sg

### Acknowledgement  
This work was supported by the IMaGIN platform (project No. 102 101 0025), through a grant from the Science and Engineering Research Council, funding to the Genome Institute of Singapore from the Agency for Science, Technology and Research (A*STAR), Singapore, and funding from the Croatian Science Foundation (Project no. UIP-11-2013-7353 - Algorithms for Genome Sequence Analysis).  

We would like to acknowledge the contribution of [Ivan Krpelnik](https://github.com/Krpa) for his help and involvement in development of the transcriptome mapping option.
