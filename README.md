# xenomapper2
A utility for splitting mixed origin NGS reads with secondary or alt mappings

Xenomapper
==========

Xenomapper2 is a utility for post processing mapped reads that have been aligned to a primary genome and a secondary genome and binning reads into species specific, multimapping in each species, unmapped and unassigned bins.  It can be used on single end or paired end sequencing data.  In paired end data evidence of sequence specificity for either read will be used to assign both reads.

Use cases include xenografts of human cancers and host pathogen interactions.

Xenomapper2 is a complete rewrite with fundamental changes in how the reads are handled internally.
Support for SAM files has been removed, and bam files are read and written using 
[pylazybam](https://github.com/genomematt/pylazybam), a pure python BAM file parser with basic write support.

Running Xenomapper2 results in the same calls for data containing only primary alignments, with the simplification that 
paired end status is automatically detected.

The major change in Xenomapper2 is treating all alignments from the same read as a group. This allows the new mode 
`--max` that selects the maximum AS score from all alignments, and the maximum XS score that occurs with this AS score.

These feature greatly enhance support for BWA-MEM based pipelines including mapping to GRCh38 with alt contigs. This
will be useful especially in the analysis of complex immune regions such as HLA in humanised mice and xenografts.

SCHEMATIC STILL TO BE UPDATED TO REFLECT
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

License
=======
To ensure the least encumbrance to all users, including those in commercial environments Xenomapper2 is now licensed 
under the BSD three clause license (rather than the GPL that was used for xenomapper 1.0)