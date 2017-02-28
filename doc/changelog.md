## GraphMap - ChangeLog

**__Version 0.4.1 -> 0.5.0__**  
Release date: 28 February 2017
- Re-implemented the index. Removed all other indexes that were previously implemented, and cleaned up the code to only use the new index (MinimizerIndex). MinimizerIndex is implemented in a separate repo added to the codebase. It also uses a hash table to store the seeds, however instead of the perfect hash as before, Google's DenseHash is used. Seeds are first compiled in a giant list (each sequence in its space, in parallel), and afterwards the list is sorted (also multithreaded). Basic statistics on seed key distribution are calculated (mean, median, standard deviation). The index also allows thresholding the amount of hits during lookup (keys with a count higher than a user-specified percentil are skipped) which is very significant for large, repetitive genomes. The index can also generate minimizers (also user specified). Index also allows for custom indexing shapes to be defined, and creates the lookup shapes automatically.  
- Changed the command line parameters to allow for new features, concretely:  
  1. Removed the parameter ```max-hits``` which is now obsolete.
  2. Added parameter ```minimizer-window``` to specify the length of the minimizer window to choose minimizers from. If equal to 1, minimizers won't be used.  
  3. Added parameter ```freq-percentil``` to specify the percentil of key occurances which will be kept. E.g. if 0.99, then 1% of most repetitive keys will be skipped. If 1.0, no filtering will be used.  
  4. Added parameter ```fly-index``` which will generate index on the fly and won't store it to disk. If the index already exists on disk, it will be loaded. To completely generate a new index on the fly, use ```--fly-index --rebuild-index```.  
  5. Renamed the parameter which was previously known as ```sensitive``` to ```double-index```.
  6. Added a composite parameter called ```-x sensitive``` which will turn off minimizers and key frequency filtering.  

- Fixed an issue with RNA-seq transcriptome mapping, where recall would be lower than expected. There was a bug when checking if alignment is sane - the check would occur *after* the alignment was converted from transcriptome space to genome space, instead still on the transcriptome. This could not have caused false positives, but definitely caused many reads to be unmapped.  
- The reimplemented index now fixes the issue of segmentation fault on the human genome.  



**__Version 0.4.0 -> 0.4.1__**  
Release date: 28 January 2017  
- Fixed the SAM headers for transcriptome mapping. In the last version, the headers corresponded to the transcriptome headers, although the alignments are in the genome space.

**__Version 0.3.2 -> 0.4.0__**  
Release date: 22 January 2017  
- GraphMap can now accept a GTF file for mapping to a transcriptome. Transcriptome is internally generated using the reference file and the GTF file, and index built from the transcriptome. Reads are then mapped to the transcriptome, and final alignments converted back to the genome coordinate space by introducing 'N' operations at splice sites.  
- Transcriptome mapping is only available in anchored alignment modes.  
- Updated Edlib to the newest version. Previous version had a bug in the traceback.  
- Recent changes in Edlib produced leading and trailing deletions in some cases. This is now handled by removing the deletions and shifting the alignment start position.  
- Fixed several (possible) memory leaks and invalid reads/writes. Generating the MD tag in SAM files had an invalid read which for some reason caused strange artifacts in CIGAR strings.  

**__Version 0.3.1 -> 0.3.2__**  
Release date: 19 December 2016  
- There were segfaults caused by recently-introduced bugs to Edlib. It has since been updated, and this version of GraphMap now includes the fixed version of Edlib.
- There was a memory leak when generating clusters.
- Minor fixes to some syntax.

**__Version 0.3.0 -> 0.3.1__**  
Release date: 12 October 2016  
- Important: Fixed MD field issues
- Minor bug fixes: composite command line parameter ```-x illumina``` depended on a parameter which wasn't defined properly, filtered empty SAM lines, etc.

**__Version 0.22 -> 0.3.0__**  
Release date: 15 April 2016  
If you are using versions 0.3.x please update to the most recent commit. There were several important memory access issues which are now resolved.  
GraphMap's command line has changed significantly between version 0.3.x and 0.2x - although many options remain similar, the usage is incompatible with the older releases due to explicit tool specification.  
The first parameter is now mandatory, and specifies whether the **mapping/alignment** (```./graphmap align```) or **overlapping** (```./graphmap owler```) should be used.  
**Important change in command line parameters.** The new version is not completely compatible to the previous one. For this reason, the minor version number has changed.  
- Changed the version numbering from: ```x.yz``` to ```x.y.z```
- Implemented a new argument parser.
- Fixed a bug with overhanging base (Issue #14), commit: 41ae30b0d8603469c62794cba1960dc42f739d4e
- Fixed the extensions of alignment to read ends when near an overhang (Issue #18).
- Fixed Issue #19 - inconsistent behaviour for parameter ```-F```.
- Cleaned up the code a bit.
- Restructured the code. Majority of the code was extracted from the repository to be used as the codebase for this and other projects. GraphMap's main code is left in this repo, while the rest is linked via git submodules.
- Added support for reading SAM and GFA files as the input sequences. Gzipped versions of all formats are supported as well. By default the format is chosen by the extension of the fila (--infmt auto), but can be specified manually.
- Added support for the M5 output format.
- Added the MD field to the SAM output.
- New and better anchor filtering (anchored modes only) using chaining of anchors that passed the LCSk.
- New and better clustering of anchor stretches. This will be used for implementing RNA-seq alignment.
- No need to precompile libraries for your system anymore. Libraries are now included in the source, or in the submodules. To initialize submodules, either clone recursively, or call ```make modules``` once GraphMap repo has been cloned.
- Anchored alignment is now the default one.  

Important command line changes:
- Long argument names are now provided.
- Extended CIGAR format can now be used via commandline through the --extcigar parameter (unlike before, where the code needed to be recompiled).
- By default, GraphMap now uses only one gapped spaced index (previously, two were used by default; one could have been used by specifying the parsimonious mode). The defaults now are the ex parsimonious mode. To use two indexes, specify the parameter: --sensitive
- The ```-w owler``` and ```-w overlapper``` have been moved. The alignment/owler mode is chosen as the first parameter in the commandline now (a "subprogram"; e.g. run ```graphmap owler```. To use the ex ```-w overlapper```, specify ```-x overlap``` instead. This mode has now been simply converted to a composite parameter. There is also a command line parameter ```--overlapper``` which only controls the counting of hits in order to skip self-hits.
- There is now a default E-value filter set at ```1e0```
- There is now a default MAPQ filter set at ```1```
- It is now possible to switch off extension of alignments to read ends (parameter: ``--no-end2end```).
- If the index needs to be rebuilt, it can now be done using a sinle command line with parameter: ```--rebuild-index``
