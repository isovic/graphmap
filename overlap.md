## GraphMap Owler - Overlap With Long Erroneous Reads
GraphMap implements two overlap modes:  
- ```./graphmap owler``` - fast, uses a trimmed GraphMap pipeline, reports output in MHAP or PAF formats, and  
- ```./graphmap align -x overlap``` - full GraphMap pipeline including alignment, output in SAM format.  
  
Owler mode (Overlap With Long Erroneous Reads) skips the graph-mapping and alignment steps. The full pipeline consists of the following steps:  
1. Construct a gapped spaced index of the reads for only one shape (6-mers, "1111110111111").  
2. For a read, collect all gapped spaced seed hits.  
3. LCSk++.  
4. Filtering seeds reported by LCSk++.  
5. Output overlaps in MHAP-like or PAF format. For details, see below.  

Currently, no seed hits are discarded, which can make overlapping slow on larger or more repetitive datasets, but very sensitive.  

Note that the overlappers are still experimental, and require thorough testing.  

### Output formats
**MHAP** format is described here: [http://mhap.readthedocs.org/en/latest/quickstart.html#output](http://mhap.readthedocs.org/en/latest/quickstart.html#output).  
GraphMap's output uses the same columns, but the meaning of columns 3 and 4 (```Jaccard score``` and ```# shared min-mers``` respectively) is different in our context.  
Instead of ```Jaccard score``` the fraction of bases covered by seeds is reported.  
Instead of ```# shared min-mers``` the number of seeds which survived filtering is reported.  

GraphMap can also output overlaps to **PAF** format. Specification of the format can be found here: [https://github.com/lh3/miniasm/blob/master/PAF.md](https://github.com/lh3/miniasm/blob/master/PAF.md).  

### Comparison to other methods  
We are working on scripts to benchmark various overlapping tools on simulated and real (later) data.  
An initial functioning version can be found here: [https://github.com/isovic/overlap-benchmark](https://github.com/isovic/overlap-benchmark).  

### Examples   
```  
# Overlap all reads from a given FASTA/FASTQ file and report overlaps in MHAP format (fast):  
./graphmap owler -r reads.fa -d reads.fa -o overlaps.mhap  

# Overlap all reads from a given FASTA/FASTQ file and report overlaps in PAF format:  
./graphmap owler -r reads.fa -d reads.fa -o overlaps.paf -L paf  

# Overlap all reads from a given FASTA/FASTQ in a full GraphMap mode with generating alignments (slow):  
./graphmap align -x overlap -r reads.fa -d reads.fa -o overlaps.sam  
```  
