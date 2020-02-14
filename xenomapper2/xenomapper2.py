#!/usr/bin/env python3
# encoding: utf-8
"""
xenomapper2.py

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

import sys, gzip
from typing import List, Tuple
from collections import Counter
from itertools import zip_longest

from pylazybam.bam import *

__author__ = "Matthew Wakefield"
__copyright__ = ("Copyright 2018-2020 Matthew Wakefield"
                 "The Walter and Eliza Hall Institute and "
                 "The University of Melbourne")
__credits__ = ["Matthew Wakefield",]
__license__ = "BSD-3-Clause"
__version__ = "2.0a1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development/Beta"

MIN32INT: int = -2147483648

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
        for align in self.alignments:
            name = get_raw_read_name(align, get_len_read_name(align))
            if previous_name and name != previous_name:
                yield alignbatch
                alignbatch = []
            alignbatch.append(align)
            previous_name = name
        if alignbatch:
            yield alignbatch

    def __next__(self):
        return next(self.alignment_batches)


def calc_cigar_based_score(cigar_string: bytes,
                         NM: int = None,
                         mismatch: int = -6,
                         gap_open: int = -5,
                         gap_extend: int = -3,
                         softclip: int = -2) -> int:
    """Return calculated values for a score equivalent to a
    rescaled AS tag using the cigar line and NM tags from the
    SAM entry.
    Score is rescaled such that a perfect match for the full
    length of the read acchieves a score of 0.0.
    Penalties are: mismatch -6, gap open -5, gap extend -3,
                   softclipped -2 (equiv to +2 match rescaled)
    Parameters
    ----------
    align : bytes

    mismatch : int

    gap_open : int

    gap_extend : int

    softclip : int

    Returns
    -------
    int
        The AS score calculated from the cigar string
    """
    if cigar_string == '*' or NM == None:
        return MIN32INT #either a multimapper or unmapped
    mismatches = NM
    cigar = re.findall(r'([0-9]+)([MIDNSHPX=])', cigar_string)
    deletions = [int(x[0]) for x in cigar if x[1] == 'D']
    insertions = [int(x[0]) for x in cigar if x[1] == 'I']
    softclips = [int(x[0]) for x in cigar if x[1] == 'S']
    score = (mismatch * mismatches) + (
            gap_open * (len(insertions) + len(deletions))) \
            + (gap_extend * (sum(insertions) + sum(deletions))) \
            + (softclip * sum(softclips))
    return score

def get_cigar_based_score(align: bytes,
                          no_tag: int = MIN32INT) -> int:
    cigar_string = decode_cigar(get_raw_cigar(align,
                                            get_len_read_name(align),
                                            get_number_cigar_operations(align),
                                            )
                                )
    NM = get_int_tag(align, b'NM', None)
    if cigar_string == None or cigar_string == '*' or NM == None:
        return no_tag #either a multimapper or unmapped
    else:
        return calc_cigar_based_score(cigar_string, NM)


def always_zero(*args,**kwargs):
    """A function that returns zero for any combination of parameters"""
    return 0

def always_very_negative(*args,**kwargs):
    """A function that returns zero for any combination of parameters"""
    return MIN32INT

def split_forward_reverse(alignments: List[bytes]
                          ) -> Tuple[List[bytes], List[bytes]]:
    """
    Split a list of BAM alignments into forward and reverse

    Parameters
    ----------
    alignments : List[bytes]
        a list of BAM alignments in raw binary format containing forward
        and optionally reverse alignments

    Returns
    -------
    (List[bytes], List[bytes])
        a tuple containing a list of forward and reverse BAM reads

    Raises
    ------
    ValueError
        if there are any reads that do not have the forward or reverse flag
    """
    forward = []
    reverse = []
    for align in alignments:
        if is_flag(align, FLAGS['forward']):
            forward.append(align)
        elif is_flag(align, FLAGS['reverse']):
            reverse.append(align)
        else: #pragma: no cover
            raise ValueError(f"The alignment {align} has neither the forward"
                             "or the reverse flag set")
    return (forward, reverse)

def get_bamprimary_AS_XS(alignments: List[bytes],
                      AS_function= bam.get_AS,
                      XS_function= bam.get_XS,
                      ) -> Tuple[int,int] :
    #named bamprimary to distinguish from primary/secondary species alignment
    """Get the AS and XS of the primary alignment

    Parameters
    ----------
    alignments : List[bytes]
        a list of binary format BAM alignments from the same read.
        Should contain only forward or reverse reads and only one primary flag

    AS_function
        a function that accepts a BAM alignment bytestring and returns the AS
        score (alignment score) as an integer
        eg pylazybam.bam.get_AS or

    XS_function
        a function that accepts a BAM alignment bytestring and returns the XS
        score as an integer

    Returns
    -------
    Tuple[int,int]
        a tuple of the alignment score AS and suboptimal alignment score XS

    Raises
    ------
    ValueError
        raises value error if there is no primary or more than one primary flag

    """
    bamprimary = [a for a in alignments if not is_flag(a,FLAGS['secondary'])]
    # We keep unmapped reads as this is an expected state in secondary species
    if len(bamprimary) == 1:
        AS = AS_function(bamprimary[0])
        XS = XS_function(bamprimary[0])
        return (AS,XS)
    elif len(bamprimary) > 1:
        raise ValueError("Multiple primary alignments in alignment batch")
    elif len(bamprimary) == 0:
        raise ValueError("No primary alignments in alignment batch")


def get_max_AS_XS(alignments: List[bytes],
                  AS_function=bam.get_AS,
                  XS_function=bam.get_XS,
                  sort_key=None,
                  default=MIN32INT,
                  ) -> Tuple[int, int]:
    """Get the maximum AS and XS from a group of alignments

    Parameters
    ----------
    alignments : List[bytes]
        a list of binary format BAM alignments.

    AS_function : Callable[[bytes], int]
        a function that accepts a BAM alignment bytestring and returns the AS
        score (alignment score) as an integer
        eg pylazybam.bam.get_AS or

    XS_function : Callable[[bytes], int]
        a function that accepts a BAM alignment bytestring and returns the XS
        score as an integer

    sort_key : Callable[[Tuple[int,int]], Any]
        a key function that converts tuples to items to be used by sorted()

    Returns
    -------
    Tuple[int,int]
        a tuple of the alignment score AS and suboptimal alignment score XS
    """
    if not alignments:
        return (default,default)
    scores = sorted([(AS_function(a),XS_function(a)) for a in alignments],
                    key=sort_key)
    return scores[-1]


def get_mapping_state(AS1, XS1, AS2, XS2, min_score=MIN32INT):
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
    if AS1 <= min_score and AS2 <= min_score:  # low quality mapping in both
        return 'unassigned'
    elif AS1 > min_score and (AS2 <= min_score or AS1 > AS2):
        # maps in primary better than secondary
        if not XS1 or AS1 > XS1:
            # maps uniquely in primary better than secondary
            return 'primary_specific'
        else:
            # multimaps in primary better than secondary
            return 'primary_multi'
    elif AS1 == AS2:  # maps equally well in both
        return 'unresolved'
    elif AS2 > min_score and (AS1 <= min_score or AS2 > AS1):
        # maps in secondary better than primary
        if (not XS2) or AS2 > XS2:
            return 'secondary_specific'
        else:
            # multimaps in secondary better than primary
            return 'secondary_multi'
    else: #pragma: no cover
        raise RuntimeError(
            f"Error in processing logic with values {(AS1, XS1, AS2, XS2)}"
        )


class DummyFile(): #pragma: no cover
    """A dummy io class looks like pylazybam.bam.FileWriter that does nothing"""
    def __init__(self, *args, **kwargs):
        self.raw_header = None
        self.raw_refs = None
        self.header_written = False
        pass

    def __repr__(self):
        return "XenomapperDummyFile"

    def write(self, data):
        pass

    def close(self):
        pass

    def writable(self):
        return True

    def seekable(self):
        return False

    def get_updated_header(self, *args, **kwargs):
        pass

    def update_header(self, *args, **kwargs):
        pass

    def update_header_length(self, *args, **kwargs):
        pass

    def write_header(self, *args, **kwargs):
        pass

    def __enter__(self):
        pass

    def __exit__(self):
        pass


class XenomapperOutputWriter():
    def __init__(self,
                 primary_raw_header,
                 primary_raw_refs,
                 secondary_raw_header,
                 secondary_raw_refs,
                 primary_specific=None,
                 primary_multi=None,
                 secondary_specific=None,
                 secondary_multi=None,
                 unresolved=None,
                 unassigned=None,
                 basename=None,
                 cmdline='',
                 compresslevel=6
                 ):
        file_arguments = dict(list(locals().items())[5:11])

        if basename == None:
            self._fileobjects = {'unassigned': DummyFile(),
                                'primary_specific': DummyFile(),
                                'primary_multi': DummyFile(),
                                'unresolved': DummyFile(),
                                'secondary_specific': DummyFile(),
                                'secondary_multi': DummyFile(),
                                }
            for key in file_arguments:
                if file_arguments[key]:
                    self._fileobjects[key] = bam.FileWriter(file_arguments[key],
                                                    compresslevel=compresslevel)
        else:
            self._fileobjects = {}
            for key in file_arguments:
                self._fileobjects[key] = bam.FileWriter(f"{basename}_{key}",
                                                    compresslevel=compresslevel)

        self.write_headers(primary_raw_header,
                             primary_raw_refs,
                             secondary_raw_header,
                             secondary_raw_refs,
                             cmdline=cmdline)


    def write_headers(self, primary_raw_header,
                             primary_raw_refs,
                             secondary_raw_header,
                             secondary_raw_refs,
                             cmdline=''):
        # write out the correct species header to all of the files
        # use primary for unresolved and unassigned
        # add @PO and @CO fields to the header
        for key in self._fileobjects:
            if key in ['secondary_specific','secondary_multi']:
                self._fileobjects[key].raw_header = secondary_raw_header
                self._fileobjects[key].raw_refs = secondary_raw_refs
            else:
                self._fileobjects[key].raw_header = primary_raw_header
                self._fileobjects[key].raw_refs = primary_raw_refs
            self._fileobjects[key].update_header(id = 'xenomapper',
                                                    program = 'xenomapper',
                                                    version = __version__,
                                                    description= key,
                                                    command=cmdline,)
            comment = f"@CO\tThis file contains alignments assigned as {key}\n"
            self._fileobjects[key].raw_header += comment.encode('utf-8')
            self._fileobjects[key].update_header_length()
            self._fileobjects[key].write_header()


    def __getitem__(self, key):
        return self._fileobjects[key]

    def keys(self):
        return self._fileobjects.keys()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        for fileobj in self._fileobjects:
            self._fileobjects[fileobj].close()


def xenomap_states(primary_aligns,
                   secondary_aligns,
                   score_function = get_bamprimary_AS_XS,
                   AS_function = get_AS,
                   XS_function = get_XS,
                   min_score=MIN32INT,
                   ):
    primary_name = get_raw_read_name(primary_aligns[0],
                                     get_len_read_name(primary_aligns[0]))
    secondary_name = get_raw_read_name(secondary_aligns[0],
                                     get_len_read_name(secondary_aligns[0]))
    if primary_name != secondary_name:
        raise ValueError("Primary and secondary read names do not match: "
                         f"{primary_name} != {secondary_name}")
    prim_f_aligns, prim_r_aligns = split_forward_reverse(primary_aligns)
    sec_f_aligns, sec_r_aligns = split_forward_reverse(secondary_aligns)
    prim_f_AS, prim_f_XS = score_function(prim_f_aligns,
                                          AS_function=AS_function,
                                          XS_function=XS_function)
    sec_f_AS, sec_f_XS = score_function(sec_f_aligns,
                                        AS_function = AS_function,
                                        XS_function = XS_function)
    forward_state = get_mapping_state(prim_f_AS, prim_f_XS,
                                      sec_f_AS, sec_f_XS,
                                      min_score)

    if not prim_r_aligns and not sec_r_aligns:
        reverse_state = None
    else:
        prim_r_AS, prim_r_XS = score_function(prim_r_aligns,
                                                AS_function = AS_function,
                                                XS_function = XS_function)
        sec_r_AS, sec_r_XS = score_function(sec_r_aligns,
                                            AS_function = AS_function,
                                            XS_function = XS_function)
        reverse_state = get_mapping_state(prim_r_AS, prim_r_XS,
                                          sec_r_AS, sec_r_XS,
                                          min_score)
    return forward_state, reverse_state


def state_map(forward_state, reverse_state):
    # logic is directly copied from xenomapper 1.0
    if reverse_state == None:
        return forward_state
    elif forward_state == 'primary_specific' or reverse_state == 'primary_specific':
        return 'primary_specific'
    elif forward_state == 'secondary_specific' or reverse_state == 'secondary_specific':
        return 'secondary_specific'
    elif forward_state == 'primary_multi' or reverse_state == 'primary_multi':
        return 'primary_multi'
    elif forward_state == 'secondary_multi' or reverse_state == 'secondary_multi':
        return 'secondary_multi'
    elif forward_state == 'unresolved' or reverse_state == 'unresolved':
        return 'unresolved'
    elif forward_state == 'unassigned' or reverse_state == 'unassigned':
        return 'unassigned'
    else:
        raise ValueError(f'Unexpected states forward:{forward_state}'
                         f'reverse:{reverse_state}')

def conservative_state_map(forward_state, reverse_state):
    # logic is directly copied from xenomapper 1.0
    if reverse_state == None:
        return forward_state
    elif forward_state == 'unassigned' or reverse_state == 'unassigned':
        return 'unassigned'
    elif forward_state == 'unresolved' or reverse_state == 'unresolved' \
            or (forward_state in ['primary_specific', 'primary_multi'] and \
                reverse_state in ['secondary_specific', 'secondary_multi']) \
            or (forward_state in ['secondary_specific', 'secondary_multi'] and \
                reverse_state in ['primary_specific', 'primary_multi']):
        return 'unresolved'
    elif forward_state == 'primary_specific' or reverse_state == 'primary_specific':
        return 'primary_specific'
    elif forward_state == 'secondary_specific' or reverse_state == 'secondary_specific':
        return 'secondary_specific'
    elif forward_state == 'primary_multi' or reverse_state == 'primary_multi':
        return 'primary_multi'
    elif forward_state == 'secondary_multi' or reverse_state == 'secondary_multi':
        return 'secondary_multi'
    else:
        raise ValueError(f'Unexpected states forward:{forward_state} '
                         f'reverse:{reverse_state}')  # pragma: no cover

def xenomap(primary_bam,
                secondary_bam,
                output_writer,
                score_function: get_bamprimary_AS_XS,
                AS_function=get_AS,
                XS_function=get_XS,
                min_score=MIN32INT,
                conservative = False,
                ):
    """

    Parameters
    ----------
    primary_bam
    secondary_bam
    output_writer
    score_function
    min_score
    conservative

    Returns
    -------
    category_pair_counts,
    category_counts,
    output_writer
    """
    category_pair_counts = Counter()
    category_counts = { 'primary_specific' : 0,
                        'secondary_specific' : 0,
                        'primary_multi' : 0,
                        'secondary_multi' : 0,
                        'unresolved' : 0,
                        'unassigned' : 0,
                        }
    if primary_bam.sort_order == 'coordinate':
        raise ValueError(f"{primary_bam._ubam.name} is coordinate sorted"
                         "BAM files must be ordered by readname for xenomapper")

    if secondary_bam.sort_order == 'coordinate':
        raise ValueError(f"{secondary_bam._ubam.name} is coordinate sorted"
                         "BAM files must be ordered by readname for xenomapper")


    for primary_aligns, secondary_aligns in zip_longest(primary_bam,
                                                        secondary_bam,
                                                        fillvalue = None):
        # TODO This code probably does not do what I think it should
        # Need to work out how to properly test for uneven file lengths
        #if primary_aligns == None or secondary_aligns == None:
        #    raise ValueError("BAM files are unequal lengths")
        forward_state, reverse_state = xenomap_states(primary_aligns,
                                                  secondary_aligns,
                                                  score_function,
                                                  AS_function=AS_function,
                                                  XS_function=XS_function,
                                                  min_score = min_score)
        category_pair_counts[(forward_state, reverse_state)] += 1
        if conservative:
            category = conservative_state_map(forward_state, reverse_state)
        else:
            category = state_map(forward_state, reverse_state)
        category_counts[category] += 1
        if category in ['secondary_specific', 'secondary_multi']:
            for align in secondary_aligns:
                output_writer[category].write(align)
        else:
            for align in primary_aligns:
                output_writer[category].write(align)

    return category_pair_counts, category_counts, output_writer

def output_summary(category_counts,
                   title = 'Read Count Category Summary\n',
                   outfile=sys.stderr):
    print('-' * 80, file=outfile)
    print(file=outfile)
    print(title, file=outfile)
    print(file=outfile)
    print('|       {0:45s}|     {1:10s}  |'.format('Category', 'Count'),
          file=outfile)
    print('|:', '-' * 50, ':|:', '-' * 15, ':|', sep='', file=outfile)
    for category in sorted(category_counts):
        if type(category) != str:
            category_name = ' & '.join(category)
        else:
            category_name = category
        print('|  {0:50s}|{1:15d}  |'.format(category_name,
                                             category_counts[category]),
              file=outfile)
    print(file=outfile)
    pass

