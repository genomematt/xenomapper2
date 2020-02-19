[![Build Status](https://travis-ci.org/genomematt/xenomapper2.svg?branch=master)
](https://travis-ci.org/genomematt/xenomapper2)
[![Coverage Status](https://coveralls.io/repos/genomematt/xenomapper/badge.svg)
](https://coveralls.io/r/genomematt/xenomapper2)
[![pypi](https://img.shields.io/pypi/v/xenomapper2.svg)](https://pypi.python.org/pypi/xenomapper2/)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![PyPy 3.8](https://img.shields.io/badge/PyPy-3.8-brightgreen)](https://www.pypy.org)

# Xenomapper2
![logo](/logo.png "logo" =200x200)

A utility for splitting mixed origin NGS reads with secondary or alt mappings

Xenomapper2
==========

Xenomapper2 is a utility for post processing mapped reads that have been aligned to a primary genome and a secondary 
genome and binning reads into species specific, multimapping in each species, unmapped and unassigned bins.  
It can be used on single end or paired end sequencing data.  
In paired end data evidence of sequence specificity for either read can be used to assign both reads.

Use cases include xenografts of human cancers and host pathogen interactions.

Xenomapper2 is a complete rewrite of [xenomapper](https://github.com/genomematt/xenomapper) with fundamental changes in 
how the reads are handled internally.
Support for SAM files has been removed, and BAM files are read and written using 
[pylazybam](https://github.com/genomematt/pylazybam), a pure python BAM file parser with basic write support.

Xenomapper2 can be used with most common aligners including 
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), 
[HISAT2](https://ccb.jhu.edu/software/hisat2/), and 
[BWA-MEM](https://github.com/lh3/bwa).

Running Xenomapper2 generally results in the same calls for data containing only primary alignments, with the 
simplification that paired end status is automatically detected.

Xenomapper2 improves handling of cases where the second read of a template is unmapped, and the mapped read is reported 
first. Xenomapper 1.0 would interpret these cases as the first entry in the SAM file as being the forward read, 
creating a mismatch in the comparison of forward and reverse between primary and secondary. By using flag status this 
potential source of error is now eliminated.

The major change in Xenomapper2 is treating all alignments from the same read as a group. This allows the new mode 
`--max` that selects the maximum AS score from all alignments, and the maximum XS score that occurs with this AS score.
This allows assignment based on high scoring secondary mappings.

These feature greatly enhance support for BWA-MEM based pipelines including mapping to GRCh38 with alt contigs. This
will be useful especially in the analysis of complex immune regions such as HLA in humanised mice and xenografts.

![Schematic of Xenomapper Use](/schematic.jpg "Schematic of Xenomapper Use")

Installation
============
Xenomapper2 requires python 3.6 or higher and is tested on linux and MacOS with CPython and pypy3.

<!--Installing from the Python Package Index with pip is the easiest option:

    pip3 install xenomapper2
    
Alternatively -->if you would like to install from the github repository this can be done by cloning the repository

    git clone https://github.com/genomematt/xenomapper2
    pip3 install --upgrade xenomapper2
    
or directly with pip

    pip3 install git+https://github.com/genomematt/xenomapper2.git
	
Although the repository tests by continuous integration with TravisCI its good practice to run the tests locally and 
check your install works correctly.  The tests are run with the following command:

    python3 -m xenomapper.tests.test_all

Usage
=====

    Usage:
      xenomapper2 --primary=<file>  --secondary=<file>
                  [ --primary-specific=<file> --primary-multi=<file>
    .              --secondary-specific=<file> --secondary-multi=<file>
    .              --unassigned=<file> --unresolved=<file>
                  | --basename=<str> ]
                  [ --min-score=<int> ]
                  [ --zs | --cigar]
                  [ --max ]
                  [ --conservative ]
      xenomapper2 --version
      xenomapper2 [ -h | --help ]
    
    Options:
      -h --help                  Show this screen.
      --version                  Show version.
    
      Input files
      --primary=<file>           A BAM format file of primary species alignments
      --secondary=<file>         A BAM format file of secondary species alignments
    
      Output options
      --primary-specific=<file>  filename for primary specific unique alignments
      --primary-multi=<file>     filename for primary specific multimap alignments
      --secondary-specific=<file> filename for secondary specific unique alignments
      --secondary-multi=<file>   filename for secondary specific multimap alignments
      --unassigned=<file>        filename for unassigned alignments
      --unresolved=<file>        filename for unresolved alignments
      --basename=<str>           prefix for creating all other output files
                                 only valid if no other output options provided
    
      Processing options
      --min-score=<int>          minimum AS score required. Lower scores unassigned.
                                 [ Default : None (implemented as -2^31) ]
      --zs                       use ZS scores for spliced aligner (HISAT2)
      --cigar                    use cigar scores to calculate AS score
      --max                      use the maximum score for any alignment
                                 [ Default : Use score of primary alignment ]
      --conservative             require both ends of paired reads to support the
                                 assignment

Note that unlike prior xenomapper versions there is no --pair option as forward
and reverse reads are automatically extracted based on their flag. Files of
mixed paired and single end reads are now fully supported.

Most of the time you will want to invoke:

    xenomapper2  --primary <primary.bam> --secondary <secondary.bam> --basename <prefix>
    
This will produce six output files (eg prefix_primary_specific.bam) and print a summary to standard error.

Process substitution will allow the use of input SAM files

    xenomapper2 --primary <(samtools view -bS input.sam)
     
A worked example of using xenomapper can be found in [example_usage.ipynb](example_usage.ipynb)

Contributing to Xenomapper2
=========================
Xenomapper2 is licensed under the BSD three clause license.  You are free to fork this repository under the terms of 
that license.  If you have suggested changes please start by raising an issue in the issue tracker.  Pull requests are 
welcome and will be included at the discretion of the author, but must have 100% test coverage.

Bug reports should be made to the issue tracker.  Difficulty in understanding how to use the software is a documentation
 bug, and should also be raised on the issue tracker and will be tagged `question` so your question and my response are 
easily found by others.

Xenomapper2 uses numpy style docstrings, python type annotations, Travis CI, coverage and coveralls. All code should be
compatible with python versions >= 3.6 and contain only pure python code.

License
=======
To ensure the least encumbrance to all users, including those in commercial environments Xenomapper2 is now licensed 
under the BSD three clause license (rather than the GPL that was used for xenomapper 1.0)

Acknowledgements
================
Justin Bedo, Alan Rubin and Tony Papenfuss provided helpful suggestions, early testing and code review.
This work was supported by the Stafford Fox Medical Research Foundation