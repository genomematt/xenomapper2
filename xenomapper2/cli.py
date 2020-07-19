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
"""

import argparse, textwrap, time
from typing import Counter
from xenomapper2.xenomapper2 import *

__author__ = "Matthew Wakefield"
__copyright__ = ("Copyright 2018-2020 Matthew Wakefield"
                 "The Walter and Eliza Hall Institute and "
                 "The University of Melbourne")
__credits__ = ["Matthew Wakefield",]
__license__ = "BSD-3-Clause"
__version__ = "2.0rc2"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development/Beta"

__author__ = "Matthew Wakefield"
__copyright__ = ("Copyright 2011-2020 Matthew Wakefield"
                 "The Walter and Eliza Hall Institute and "
                 "The University of Melbourne")
__credits__ = ["Matthew Wakefield",]
__license__ = "BSD-3-Clause"
__version__ = "2.0rc2"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development/Beta"

MIN32INT: int = -2147483648


def command_line_interface(arguments: str = None):  # pragma: no cover
    parser = argparse.ArgumentParser(prog="xenomapper2",
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent("""\
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
                    """),

                        epilog=textwrap.dedent("""\
                    
                    Note that unlike prior xenomapper versions there is no --pair option as forward
                    and reverse reads are automatically extracted based on their flag. Files of
                    mixed paired and single end reads are now fully supported.

                    License BSD-3-Clause
                    
                    This program is distributed in the hope that it will be useful,
                    but WITHOUT ANY WARRANTY; without even the implied warranty of
                    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n

                    """),
                                     )
    input_group = parser.add_argument_group(title='Input Files')
    input_group.add_argument('--primary',
                        type=str,
                        default=None,
                        help='A BAM format file of primary species alignments \
                             Entries must be name sorted')
    input_group.add_argument('--secondary',
                        type=str,
                        default=None,
                        help='A BAM format file of secondary species alignments \
                             Entries must be name sorted')

    output_group = parser.add_argument_group(title='Output Files')
    output_group.add_argument('--primary_specific',
                        type=str,
                        default=None,
                        help='name for BAM format output file for reads mapping to a specific location in the primary species')
    output_group.add_argument('--secondary_specific',
                        type=str,
                        default=None,
                        help='name for BAM format output file for reads mapping to a specific location in the secondary species')
    output_group.add_argument('--primary_multi',
                        type=str,
                        default=None,
                        help='name for BAM format output file for reads multi mapping in the primary species')
    output_group.add_argument('--secondary_multi',
                        type=str,
                        default=None,
                        help='name for BAM format output file for reads multi mapping in the secondary species')
    output_group.add_argument('--unassigned',
                        type=str,
                        default=None,
                        help='name for BAM format output file for unassigned (non-mapping) reads')
    output_group.add_argument('--unresolved',
                        type=str,
                        default=None,
                        help='name for BAM format output file for unresolved (maps equally well in both species) reads')
    output_group.add_argument('--basename',
                        type=str,
                        default=None,
                        help='prefix for creating all output files \
                             only valid if no other output options provided')
    processing_group = parser.add_argument_group(title='Processing options')
    processing_group.add_argument('--conservative',
                        action='store_true',
                        help='conservatively allocate paired end reads with discordant category allocations. \
                              Only pairs that are both specific, or specific and multi will be allocated as specific. \
                              Pairs that are discordant for species will be deemed unresolved.  Pairs where any read \
                              is unassigned will be deemed unassigned.')
    processing_group.add_argument('--max',
                        action='store_true',
                        help='use the maximum score for any alignment \
                              [Default: Use score of primary alignment]')
    processing_group.add_argument('--min_score',
                        type=int,
                        default=MIN32INT,
                        help='the minimum mapping score.  Reads with scores less than or equal to min_score will be considered unassigned. \
                              Values should be chosen based on the mapping program and read length')
    score_functions = processing_group.add_mutually_exclusive_group()
    score_functions.add_argument('--cigar',
                        action='store_true',
                        help='Use the cigar line and the NM tag to calculate a score. For aligners that do not support the AS tag. \
                              No determination of multimapping state will be done.  Reads that are unique in one species and multimap \
                              in the other species may be misassigned as no score can be calculated in the multimapping species. \
                              Score is -6 * mismatches + -5 * indel open + -3 * indel extend + -2 * softclip. \
                              Treatment of multimappers will vary with aligner.  If multimappers are assigned a cigar line they \
                              will be treated as species specific, otherwise as unassigned.')
    score_functions.add_argument('--zs',
                        action='store_true',
                        help='Use the value of the ZS tag in place of XS for determining the mapping score of the next best \
                              alignment.  Used with HISAT as the XS:A tag is conventionally used for strand in spliced mappers.')
    parser.add_argument('--version',
                        action='store_true',
                        help='print version information and exit')
    if arguments:
        args = parser.parse_args(arguments.split(' '))
    else:
        args = parser.parse_args()
    if args.version:
        print(f"xenomapper2 v{__version__}")
        sys.exit()
    elif not args.primary or not args.secondary:
        print('ERROR: You must provide both a primary and a secondary alignment'
              ' in valid bam format\n'
              f'Primary:   {args.primary}\n'
              f'Secondary: {args.secondary}\n'
              'see -h or --help for usage details\n',
              file=sys.stderr)
        sys.exit()
    return args

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
        args = command_line_interface()
    else:
        args = command_line_interface(arguments)

    if args.cigar:
        AS_function = get_cigar_based_score
        XS_function = always_very_negative

    elif args.zs:
        AS_function = bam.get_AS
        XS_function = bam.get_ZS
    else:
        AS_function = bam.get_AS
        XS_function = bam.get_XS

    if args.max:
        score_function = get_max_AS_XS
    else:
        score_function = get_bamprimary_AS_XS

    if args.min_score:
        min_score = int(args.min_score)
    else:
        min_score = MIN32INT

    primary_bam = AlignbatchFileReader(gzip.open(args.primary))
    primary_header = primary_bam.raw_header
    primary_refs = primary_bam.raw_refs

    secondary_bam = AlignbatchFileReader(gzip.open(args.secondary))
    secondary_header = secondary_bam.raw_header
    secondary_refs = secondary_bam.raw_refs

    cmdline = " ".join(
            [ f"{x}={vars(args)[x]}" for x in list(vars(args).keys())[:-1]])

    print(f"\nxenomapper2 v{__version__} {cmdline}\n", file=output)
    #print(f"\nxenomapper2 v{__version__} {cmdline}\n", file=sys.stderr)

    xow = XenomapperOutputWriter(primary_header,
                                 primary_refs,
                                 secondary_header,
                                 secondary_refs,
                                 primary_specific=args.primary_specific,
                                 secondary_specific=args.secondary_specific,
                                 primary_multi=args.primary_multi,
                                 secondary_multi=args.secondary_multi,
                                 unassigned=args.unassigned,
                                 unresolved=args.unresolved,
                                 basename=args.basename,
                                 cmdline=cmdline
                                 )

    pair_counts, counts, writer = xenomap(primary_bam,
                                          secondary_bam,
                                          output_writer=xow,
                                          score_function=score_function,
                                          AS_function=AS_function,
                                          XS_function=XS_function,
                                          min_score=min_score,
                                          conservative=args.conservative,
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