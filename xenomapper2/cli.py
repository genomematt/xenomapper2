#!/usr/bin/env python3
# encoding: utf-8
"""
██╗  ██╗███████╗███╗   ██╗ ██████╗
╚██╗██╔╝██╔════╝████╗  ██║██╔═══██╗
 ╚███╔╝ █████╗  ██╔██╗ ██║██║   ██║
 ██╔██╗ ██╔══╝  ██║╚██╗██║██║   ██║
██╔╝ ██╗███████╗██║ ╚████║╚██████╔╝
╚═╝  ╚═╝╚══════╝╚═╝  ╚═══╝ ╚═════╝
███╗   ███╗ █████╗ ██████╗ ██████╗ ███████╗██████╗   ██████╗
████╗ ████║██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔══██╗  ╚════██╗
████╗ ████║██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔══██╗  ╚════██╗
██╔████╔██║███████║██████╔╝██████╔╝█████╗  ██████╔╝   █████╔╝
██║╚██╔╝██║██╔══██║██╔═══╝ ██╔═══╝ ██╔══╝  ██╔══██╗  ██╔═══╝
██║ ╚═╝ ██║██║  ██║██║     ██║     ███████╗██║  ██║  ███████╗
╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝     ╚══════╝╚═╝  ╚═╝  ╚══════╝
xenomapper2

A script for parsing pairs of BAM format files of reads aligned to two different
genomes and returning BAM files containing only reads where no better mapping
exist in other genome (and optionally outputting additional categories).
Used for filtering reads where multiple species may contribute (eg human cancer
xenografted into mouse, humanised mice, and co-sequenced pathogen and host).
Designed to work with single and paired end reads with or without secondary/alt
alignments.


Created by Matthew Wakefield.
Copyright (c) 2011-2020  Matthew Wakefield
The Walter and Eliza Hall Institute and The University of Melbourne.
All rights reserved.


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

The dots at the start of the useage line have no meaning. They are just there to
keep the argument parser (docopt) from disliking lines that start with --


License BSD-3-Clause

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""

import sys, gzip, time
from docopt import docopt
from typing import Counter, Tuple

from pylazybam import bam
from xenomapper2.xenomapper2 import *

__author__ = "Matthew Wakefield"
__copyright__ = ("Copyright 2018-2020 Matthew Wakefield"
                 "The Walter and Eliza Hall Institute and "
                 "The University of Melbourne")
__credits__ = ["Matthew Wakefield",]
__license__ = "BSD-3-Clause"
__version__ = "2.0rc1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development/Beta"


def main(arguments: str = None,
         output = sys.stderr) -> Tuple[Counter, Counter]:
    """main function for orchestrating a xenomapper2 run

    Parameters
    ----------
    arguments : str
        a string of unix command line arguments - used for testing
    output : TextIO
        destination to print progress and results to - used for testing

    Returns
    -------
    Counter
        a counter of pair states
    Counter
        a counter of read states
    """

    start_time = time.time()

    if not arguments:#pragma: no cover
        args = docopt(__doc__,
                      version= f"xenomapper2 v{__version__}")
    else:
        args = docopt(__doc__,
                      argv= arguments,
                      version= f"xenomapper2 v{__version__}")

    #if run with no arguments print help
    if not sum([bool(x) for x in args.values()]): #pragma: no cover
        docopt(__doc__, argv='--help')
        sys.exit()


    if args["--cigar"]:
        AS_function = get_cigar_based_score
        XS_function = always_very_negative

    elif args["--zs"]:
        AS_function = bam.get_AS
        XS_function = bam.get_ZS
    else:
        AS_function = bam.get_AS
        XS_function = bam.get_XS

    if args["--max"]:
        score_function = get_max_AS_XS
    else:
        score_function = get_bamprimary_AS_XS

    if args["--min-score"]:
        min_score = int(args["--min-score"])
    else:
        min_score = MIN32INT

    primary_bam = AlignbatchFileReader(gzip.open(args["--primary"]))
    primary_header = primary_bam.raw_header
    primary_refs = primary_bam.raw_refs

    secondary_bam = AlignbatchFileReader(gzip.open(args["--secondary"]))
    secondary_header = secondary_bam.raw_header
    secondary_refs = secondary_bam.raw_refs

    cmdline = " ".join([ f"{x}={args[x]}" for x in args.keys() \
                         if x not in  [".","--help"]])

    print(f"\nxenomapper2 v{__version__} {cmdline}\n", file=output)

    xow = XenomapperOutputWriter(primary_header,
                                 primary_refs,
                                 secondary_header,
                                 secondary_refs,
                                 primary_specific=args["--primary-specific"],
                                 secondary_specific=args["--secondary-specific"],
                                 primary_multi=args["--primary-multi"],
                                 secondary_multi=args["--secondary-multi"],
                                 unassigned=args["--unassigned"],
                                 unresolved=args["--unresolved"],
                                 basename=args["--basename"],
                                 cmdline=cmdline
                                 )

    pair_counts, counts, writer = xenomap(primary_bam,
                                          secondary_bam,
                                          output_writer=xow,
                                          score_function=score_function,
                                          AS_function=AS_function,
                                          XS_function=XS_function,
                                          min_score=min_score,
                                          conservative=args["--conservative"],
                                          )

    writer.close()

    output_summary(category_counts=pair_counts,
                   title='Read Category Summary',
                   outfile=output)
    output_summary(category_counts=counts,
                   title='Read Pair Category Summary',
                   outfile=output)

    end_time = time.time()

    print(f'\n\nTotal templates assigned : {sum(pair_counts.values())} '
          f'in {end_time-start_time:.2f}s\n',
          file=output)

    if arguments:
        # return data if testing
        return pair_counts, counts
    else: #pragma: no cover
        pass


if __name__ == '__main__': #pragma: no cover
    main()