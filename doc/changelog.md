## GraphMap - Changelog

**__Version on dev branch (will probably be 0.30 due to major argument parsing change)__**
- Fixed a bug with overhanging base (Issue #14), commit: 41ae30b0d8603469c62794cba1960dc42f739d4e
- Implemented a new argument parser.
- Cleaned up the code a bit.
- Restructured the code. Majority of the code was extracted from the repository to be used as the codebase for this and other projects. GraphMap's main code is left in this repo, while the rest is linked via git submodules.
- Added support for reading SAM and GFA files as the input sequences. Gzipped versions of all formats are supported as well. By default the format is chosen by the extension of the file, but can be specified manually.
- Added support for the M5 output format.
- Added the MD field to the SAM output.
- Extended CIGAR format can now be used via commandline through the --extcigar parameter (unlike before, where the code needed to be recompiled).
