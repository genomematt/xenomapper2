#!/usr/bin/env python3
# encoding: utf-8
"""
xenomapper.py

A script for parsing pairs of BAM format files of reads mapped to two different
genomes and returning BAM files containing only reads where no better mapping
exist in other files (and optionally additional categories.
Used for filtering reads where multiple species may contribute (eg human cancer
 xenografted into mouse, humanised mice, and co-sequenced pathogen and host).

Created by Matthew Wakefield.
Copyright (c) 2011-2020  Matthew Wakefield
The Walter and Eliza Hall Institute and The University of Melbourne.
All rights reserved.


   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


"""

from pylazybam.bam import *
import sys
import os
import argparse, textwrap
import subprocess
import re
from collections import Counter
from copy import copy

__author__ = "Matthew Wakefield"
__copyright__ = ("Copyright 2011-2020 Matthew Wakefield, "
                 "The Walter and Eliza Hall Institute and "
                 "The University of Melbourne")
__credits__ = ["Matthew Wakefield", ]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Beta"

# get batches of read from each bam file
# split into forward and reverse for each file
# get AS and XS score pairs for the reads for each file
#   - option one max()
#   - option two primary mapping score
#   - option three determine AS from cigar string
# Note: should act on full alignments so all data available to function


class AlignbatchFileReader(FileReader):
    """A BAM file reader that iterates batches of reads with the same name

    Parameters
    ----------
        ubam : BinaryIO
            An binary (bytes) file or stream containing a valid uncompressed
            bam file conforming to the specification.

    Yields
    ------
        alignbatch : List[bytes]
            A byte string of a bam alignment entry in raw binary format

    Attributes
    ----------
        raw_header : bytes
            The raw bytestring representing the bam header
        raw_refs : bytes
            The raw bytestring representing the reference sequences
        header : str
            The ASCII representation of the header
        refs : Dict[str:int]
            A dictionary of reference_name keys with reference_length
        ref_to_index : Dict[str:int]
            A dictionary mapping reference names to the bam numeric identifier
        index_to_ref : Dict[int:str]
            A dictionary mapping bam reference numeric identifiers to names

    Notes
    -----
        It is advisable not to call the private functions or operate directly
        on the underlying file object.

        The detailed specification for the BAM format can be found at
        https://samtools.github.io/hts-specs/SAMv1.pdf

    Example
    -------
        This class requires an uncompressed bam file as input.
        For this example we will use gzip to decompress the test file included
        as a resource in the package

        >>> import gzip
        >>> from pkg_resources import resource_stream
        >>> ubam = gzip.open('/tests/data/paired_end_testdata_human.bam'))

        A bam filereader object can then be created and headers inspected
        >>> mybam = bam.FileReader(ubam)
        >>> print(mybam.header)

        The filereader object is an iterator and yields batches of alignments
        in raw BAM binary format
        >>> align_batch = next(mybam)
        >>> for align in align_batch:
        >>>     print(mybam.index_to_ref[get_ref_index(align)],
        >>>           get_pos(align),
        >>>           get_AS(align))

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.alignment_batches = self._get_alignment_batches()

    def _get_alignment_batches(self) -> Generator[bytes, None, None]:
        previous_name = None
        alignbatch = []
        try:
            for align in self.alignments:
                name = get_raw_read_name(align, get_len_read_name(align))
                if previous_name and name != previous_name:
                    yield alignbatch
                    alignbatch = []
                alignbatch.append(align)
                previous_name = name
        except StopIteration:
            if alignbatch:
                return alignbatch
            else:
                return

    def __next__(self):
        return next(self.alignment_batches)



def get_mapping_state(AS1,XS1,AS2,XS2, min_score=float('-inf')):
    """Determine the mapping state based on scores in each species.
    Scores can be negative but better matches must have higher scores
    Arguments:
        AS1  - float or interger score of best match in primary species
        XS1  - float or interger score of other match in primary species
        AS2  - float or interger score of best match in secondary species
        XS2  - float or interger score of other match in secondary species
        min_score - the score that matches must exceed in order to be
                    considered valid matches. Note scores equalling this
                    value will also be considered not to match.
                    Default = -inf
    Returns
        state - a string of 'primary_specific', 'secondary_specific',
                'primary_multi', 'secondary_multi',
                'unresolved', or 'unassigned' indicating match state.
    """
    if AS1 <= min_score and  AS2 <= min_score:  #low quality mapping in both
        return 'unassigned'
    elif AS1 > min_score and (AS2 <= min_score or AS1 > AS2): #maps in primary better than secondary
        if not XS1 or AS1 > XS1:       #maps uniquely in primary better than secondary
            return 'primary_specific'
        else:            #multimaps in primary better than secondary
            return 'primary_multi'
    elif AS1 == AS2:                   #maps equally well in both
        return 'unresolved'
    elif AS2 > min_score and (AS1 <= min_score or AS2 > AS1): #maps in secondary better than primary
        if (not XS2) or AS2 > XS2:
            return 'secondary_specific'
        else:
            return 'secondary_multi' #multimaps in secondary better than primary
    else: raise RuntimeError('Error in processing logic with values {0} '.format((AS1,XS1,AS2,XS2))) # pragma: no cover

def output_summary(category_counts, outfile=sys.stderr):
    print('-'*80, file=outfile)
    print('Read Count Category Summary\n', file=outfile)
    print('|       {0:45s}|     {1:10s}  |'.format('Category','Count'), file=outfile)
    print('|:','-'*50,':|:','-'*15,':|',sep='', file=outfile)
    for category in sorted(category_counts):
        print('|  {0:50s}|{1:15d}  |'.format(str(category),category_counts[category]), file=outfile)
    print(file=outfile)
    pass
