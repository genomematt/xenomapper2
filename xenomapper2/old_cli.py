#!/usr/bin/env python3
# encoding: utf-8
"""
xenomapper2/old_cli.py

A script for parsing pairs of sam files and returning sam files containing only
reads where no better mapping exist in other files.
Used for filtering reads where multiple species may contribute
(eg human tissue xenografted into mouse).

Created by Matthew Wakefield.
Copyright (c) 2011-2020  Matthew Wakefield
The Walter and Eliza Hall Institute and
The University of Melbourne.
All rights reserved.


   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


"""
import argparse, textwrap, gzip
from pylazybam import bam
from xenomapper2 import *

__author__ = "Matthew Wakefield"
__copyright__ = ("Copyright 2011-2020 Matthew Wakefield"
                 "The Walter and Eliza Hall Institute and "
                 "The University of Melbourne")
__credits__ = ["Matthew Wakefield",]
__license__ = "BSD-3-Clause"
__version__ = "2.0a1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development/Beta"

MIN32INT: int = -2147483648

def command_line_interface(*args, **kw):  # pragma: no cover
    parser = argparse.ArgumentParser(prog="xenomapper_classic",
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent("""\
                    DEPRICATED - Will be removed in a later version
                    This is the old Xenomapper 1.0 command line interface sitting on top
                    of the new Xenomapper 2.0. It exists mainly for consistency testing.
                    In most cases it will give the same output.
                    --pair has no effect as pairs are automatically detected in 2.0
                    Output files will contain supplementary and unmapped reads.
                    In Xenomapper 1.0 these input files would have caused an error 
                    
                    A script for parsing pairs of sam files and returning sam files
                    containing only reads where no better mapping exist in other files.
                    Used for filtering reads where multiple species may contribute 
                    (eg human tissue xenografted into mouse, pathogen growing on plant).

                    Files should contain an AS and XS score and better matches must have
                    a higher alignment score (but can be negative).
                    Reads must be in the same order in both species.

                    In practice this is best acchieved by using Bowtie2 in --local mode.
                    If the -p option is used you must also use --reorder.

                    Limited support is provided for aligners that do not produce AS and XS
                    score tags via the --cigar_score option.

                    All input files must be seekable
                    (ie not a FIFO, process substitution or pipe)'
                    """),
                                     epilog=textwrap.dedent("""\

                    This program is distributed in the hope that it will be useful,
                    but WITHOUT ANY WARRANTY; without even the implied warranty of
                    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n

                    """),
                                     )
    parser.add_argument('--primary_bam',
                        type=str,
                        default=None,
                        help='a BAM format Bowtie2 mapping output file corresponding to the primary species of interest')
    parser.add_argument('--secondary_bam',
                        type=str,
                        default=None,
                        help='a BAM format Bowtie2 mapping output file corresponding to the secondary or contaminating species')
    parser.add_argument('--primary_specific',
                        type=str,
                        default=sys.stdout,
                        help='name for BAM format output file for reads mapping to a specific location in the primary species')
    parser.add_argument('--secondary_specific',
                        type=str,
                        default=None,
                        help='name for BAM format output file for reads mapping to a specific location in the secondary species')
    parser.add_argument('--primary_multi',
                        type=str,
                        default=None,
                        help='name for BAM format output file for reads multi mapping in the primary species')
    parser.add_argument('--secondary_multi',
                        type=str,
                        default=None,
                        help='name for BAM format output file for reads multi mapping in the secondary species')
    parser.add_argument('--unassigned',
                        type=str,
                        default=None,
                        help='name for BAM format output file for unassigned (non-mapping) reads')
    parser.add_argument('--unresolved',
                        type=str,
                        default=None,
                        help='name for BAM format output file for unresolved (maps equally well in both species) reads')
    parser.add_argument('--paired',
                        action='store_true',
                        help='DEPRICATED Now auto detected')
    parser.add_argument('--conservative',
                        action='store_true',
                        help='conservatively allocate paired end reads with discordant category allocations. \
                              Only pairs that are both specific, or specific and multi will be allocated as specific. \
                              Pairs that are discordant for species will be deemed unresolved.  Pairs where any read \
                              is unassigned will be deemed unassigned.')
    parser.add_argument('--min_score',
                        type=int,
                        default=MIN32INT,
                        help='the minimum mapping score.  Reads with scores less than or equal to min_score will be considered unassigned. \
                              Values should be chosen based on the mapping program and read length')
    parser.add_argument('--cigar_scores',
                        action='store_true',
                        help='Use the cigar line and the NM tag to calculate a score. For aligners that do not support the AS tag. \
                              No determination of multimapping state will be done.  Reads that are unique in one species and multimap \
                              in the other species may be misassigned as no score can be calculated in the multimapping species. \
                              Score is -6 * mismatches + -5 * indel open + -3 * indel extend + -2 * softclip. \
                              Treatment of multimappers will vary with aligner.  If multimappers are assigned a cigar line they \
                              will be treated as species specific, otherwise as unassigned.')
    parser.add_argument('--use_zs',
                        action='store_true',
                        help='Use the value of the ZS tag in place of XS for determining the mapping score of the next best \
                              alignment.  Used with HISAT as the XS:A tag is conventionally used for strand in spliced mappers.')
    parser.add_argument('--version',
                        action='store_true',
                        help='print version information and exit')
    args = parser.parse_args()
    if args.version:
        print(__version__)
        sys.exit()
    if (not args.primary_bam or not args.secondary_bam):
        print('ERROR: You must provide --primary_bam and --secondary_bam\n')
        parser.print_help()
        sys.exit(1)
    return args


def main():  # pragma: no cover
    args = command_line_interface()

    if args.cigar_scores:
        AS_function = get_cigarbased_score
        XS_function = always_very_negative

    elif args.use_zs:
        AS_function = bam.get_AS
        XS_function = bam.get_ZS
    else:
        AS_function = bam.get_AS
        XS_function = bam.get_XS

    primary_bam = AlignbatchFileReader(gzip.open(args.primary_bam))
    primary_header = primary_bam.raw_header
    primary_refs = primary_bam.raw_refs

    secondary_bam = AlignbatchFileReader(gzip.open(args.secondary_bam))
    secondary_header = secondary_bam.raw_header
    secondary_refs = secondary_bam.raw_refs

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
                            basename=None,
                            )

    pair_counts, counts, writer = xenomap(primary_bam,
                                        secondary_bam,
                                        output_writer = xow,
                                        score_function = get_bamprimary_AS_XS,
                                        AS_function = AS_function,
                                        XS_function = XS_function,
                                        min_score = MIN32INT,
                                        conservative = args.conservative,
                                        )

    writer.close()



    output_summary(category_counts=pair_counts, title= 'Read Category Summary')
    output_summary(category_counts=counts, title= 'Read Pair Category Summary')
    pass


if __name__ == '__main__':  # pragma: no cover
    main()
