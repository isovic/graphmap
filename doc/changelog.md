## GraphMap - ChangeLog

**__Version 0.3.0 -> 0.3.1__**  
Release date: 12 October 2016  
- Important: Fixed MD field issues
- Minor bug fixes: composite command line parameter ```-x illumina``` depended on a parameter which wasn't defined properly, filtered empty SAM lines, etc.

**__Version 0.22 -> 0.3.0__**  
Release date: 15 April 2016  
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
