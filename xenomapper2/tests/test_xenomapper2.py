#!/usr/bin/env python3
# encoding: utf-8
"""
test_xenomapper2.py

Created by Matthew Wakefield on 2012-08-16.
Copyright (c) 2011-2020 Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""

import gzip
import io
import unittest
import warnings
from tempfile import TemporaryDirectory, NamedTemporaryFile
from hashlib import sha256
from itertools import combinations, permutations

from pkg_resources import resource_stream, resource_filename
import docopt


from xenomapper2.xenomapper2 import *
from xenomapper2 import cli


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

BAMPAIR1 = [
    b'\x12\x01\x00\x00\x0c\x00\x00\x00eB\xf0\x07(&\n2\x02\x00S\x00d\x00\x00\x00\x0c\x00\x00\x00\x1eB\xf0\x07S\xff\xff\xffHWI-ST960:96:COTO3ACXX:3:1101:1220:2089\x000\x06\x00\x00\x14\x00\x00\x00HB\x18"$H\x14\x82"\x14(\x12\x88\x14AD(AD!D\x14\x11\x82B\x88A\x12"\x84AD!\x11D\x88B\x14\x84\x14"AA\x82\x12\x12!(\x12\x1f#"##!\x1f####""!!\x1e##$$$##$"#"%%\'\'\'\'\')((&\')))))))(\'(()((\'#!))))))))((()))(&\'\'))))())()(&))(\'\'%%\'%####\x1c\x10\x02ASC\xc6XSC~XNC\x00XMC\x00XOC\x00XGC\x00NMC\x00MDZ99\x00YSC\xbdYTZCP\x00',
    b'\x15\x01\x00\x00\x0c\x00\x00\x00\x1eB\xf0\x07(&\n2\x02\x00\xa3\x00d\x00\x00\x00\x0c\x00\x00\x00eB\xf0\x07\xad\x00\x00\x00HWI-ST960:96:COTO3ACXX:3:1101:1220:2089\x004\x00\x00\x00\x10\x06\x00\x00\xff\xf8\x82\x11\x11"\x14"\x84\x12!\x12\x18D\x84\x11\x12"\x84\x82\x82\x81(\x11\x11\x18\x12\x11\x11\x88\x14(DHHD\x84HB\x18"\x84H\x14\x82"\x14(\x12\x88\x02\x02\x02\x13\x13 #%\'&\'%\')()))))))))))))$&\'())))))()))()))()\'(\'\'(((\'))))()(()))"\'$$%### !##$#\x1e\x0c\x0c\x19\x1d ""###"##"\x1f""#ASC\xbdXSC2XNC\x00XMC\x01XOC\x00XGC\x00NMC\x01MDZ79C17\x00YSC\xc6YTZCP\x00']
BAMPAIR42 = [
    b'\x1a\x01\x00\x00\t\x00\x00\x00\xd1B\xd2\x07(,\x921\x04\x00c\x00d\x00\x00\x00\t\x00\x00\x00%D\xd2\x07\xb9\x01\x00\x00HWI-ST960:96:COTO3ACXX:3:1101:1591:2148\x00\x14\x00\x00\x00\xe0\x02\x00\x00\x11\x00\x00\x00@\x03\x00\x00\xf1\x88"\x88D(\x11\x14A\x82\x14\x81AD\x82\x81\x81\x81\x81\x81\x81\x81\x81\x81\x88\x88\x88\x88\x88\x82"((\x88\x88\x11\x11\x82\x12A\x84\x81\x18\x12\x18\x18\x12\x14!\x11\x02\x10\x19#### %\'\'\'\'((&(!&&$& "\'\'\'&\n \'&\'((((((((\x18#&(%(((\x1f##&(((\x1d$#\x16\r\x1d\x1d"\x1e" \x1f""\x1e"" "\x1d\x1d \x1b\x10\x1e\x1e\x1e#"#"\x1d""#""\x19\x1d\x19 """ASC\xbcXSC<XNC\x00XMC\x00XOC\x01XGC\x01NMC\x01MDZ98\x00YSC\xc8YTZCP\x00',
    b'\x0f\x01\x00\x00\t\x00\x00\x00%D\xd2\x07(,\x921\x01\x00\x93\x00d\x00\x00\x00\t\x00\x00\x00\xd1B\xd2\x07G\xfe\xff\xffHWI-ST960:96:COTO3ACXX:3:1101:1591:2148\x00@\x06\x00\x00D\x82\x14AH\x88AA"\x14"\x84\x12!\x12\x18D\x84\x11\x12"$\x82\x88\x81(\x14\x11\x12\x18\x11\x11\x11\x88\x14(DB\x18D\x84B\x12\x12B(H\x11\x82"\x18\x1e\x1f\x19\x19\x13\x07\x19\x19\x19\x1d\x1a\x1d"\x19\x0b\x0b\x0b\x1e\x1b\x12\x1e\x1a\x1a\x1a \x1f\x1a\x1a"\x1d\x1d\x1f\x19\x1e\x16\x14\x06\x13!\x17\x0e\x1a\x16\x16\x1b!\x1c\x13\'&&&$%!\x17\'%\x1a%\'\'#(\'#\x18\x1f\x1b&&\'&$"\x18&&(&\'\x1f%$$\x1b#\x16#\x1e\'# \x1a##\x1e\x1f\x1eASC\xc8XSC\x8dXNC\x00XMC\x00XOC\x00XGC\x00NMC\x00MDZ100\x00YSC\xbcYTZCP\x00']


class test_main(unittest.TestCase):
    def setUp(self):
        pass

    def test_always_zero(self):
        self.assertEqual(always_zero(),0)

    def test_always_very_negative(self):
        self.assertEqual(always_very_negative(),MIN32INT)

    def test_split_forward_reverse(self):
        expected = ([
                        b'\x12\x01\x00\x00\x0c\x00\x00\x00eB\xf0\x07(&\n2\x02\x00S\x00d\x00\x00\x00\x0c\x00\x00\x00\x1eB\xf0\x07S\xff\xff\xffHWI-ST960:96:COTO3ACXX:3:1101:1220:2089\x000\x06\x00\x00\x14\x00\x00\x00HB\x18"$H\x14\x82"\x14(\x12\x88\x14AD(AD!D\x14\x11\x82B\x88A\x12"\x84AD!\x11D\x88B\x14\x84\x14"AA\x82\x12\x12!(\x12\x1f#"##!\x1f####""!!\x1e##$$$##$"#"%%\'\'\'\'\')((&\')))))))(\'(()((\'#!))))))))((()))(&\'\'))))())()(&))(\'\'%%\'%####\x1c\x10\x02ASC\xc6XSC~XNC\x00XMC\x00XOC\x00XGC\x00NMC\x00MDZ99\x00YSC\xbdYTZCP\x00'],
                    [
                        b'\x15\x01\x00\x00\x0c\x00\x00\x00\x1eB\xf0\x07(&\n2\x02\x00\xa3\x00d\x00\x00\x00\x0c\x00\x00\x00eB\xf0\x07\xad\x00\x00\x00HWI-ST960:96:COTO3ACXX:3:1101:1220:2089\x004\x00\x00\x00\x10\x06\x00\x00\xff\xf8\x82\x11\x11"\x14"\x84\x12!\x12\x18D\x84\x11\x12"\x84\x82\x82\x81(\x11\x11\x18\x12\x11\x11\x88\x14(DHHD\x84HB\x18"\x84H\x14\x82"\x14(\x12\x88\x02\x02\x02\x13\x13 #%\'&\'%\')()))))))))))))$&\'())))))()))()))()\'(\'\'(((\'))))()(()))"\'$$%### !##$#\x1e\x0c\x0c\x19\x1d ""###"##"\x1f""#ASC\xbdXSC2XNC\x00XMC\x01XOC\x00XGC\x00NMC\x01MDZ79C17\x00YSC\xc6YTZCP\x00'])
        self.assertEqual(split_forward_reverse(BAMPAIR1), expected)

    def test_get_bamprimary_AS_XS(self):
        self.assertRaises(ValueError, get_bamprimary_AS_XS, BAMPAIR1)
        self.assertRaises(ValueError, get_bamprimary_AS_XS, [])
        self.assertEqual(get_bamprimary_AS_XS([BAMPAIR1[0], ]), (198, 126))

    def test_get_max_AS_XS(self):
        self.assertEqual(get_max_AS_XS(BAMPAIR1), (198, 126))
        self.assertEqual(get_max_AS_XS([]), (-2147483648, -2147483648))

    def test_calc_cigar_based_score(self):
        input_and_output = [(('*', None), -2147483648),
                            (('50M',0),0),
                            (('1S49M',0),-2),
                            (('50M', 2), -12),
                            (('10M1I39M', 0), -8),
                            (('10M2D38M', 0), -11),
                            (('10M1I10M1D28M', 0), -16),
                            (('10M1234N40M', 0), 0),
                            ]
        for (inpt, outpt) in input_and_output:
            self.assertEqual(calc_cigar_based_score(*inpt), outpt)

    def test_get_cigar_based_score(self):
        pass

    def test_XenomapperOutputWriter(self):
        xow_keys = ['primary_specific', 'primary_multi',
                    'secondary_specific', 'secondary_multi',
                    'unresolved', 'unassigned']

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")


            xow = XenomapperOutputWriter(b"p", b"pr", b"s", b"sr")
            self.assertEqual([repr(x) for x in list(xow._fileobjects.values())],
                             [repr(DummyFile()), ] * 6)
            tempd = TemporaryDirectory()
            xow = XenomapperOutputWriter(b"p", b"s", b"s", b"sr",
                                         basename=f'{tempd.name}/foo')
            expected = [f'{tempd.name}/{f}.bam' for f in ['foo_primary_specific',
                                                     'foo_primary_multi',
                                                     'foo_secondary_specific',
                                                     'foo_secondary_multi',
                                                     'foo_unresolved',
                                                     'foo_unassigned']]
            self.assertEqual([x.name for x in xow._fileobjects.values()],
                             expected)
            self.assertEqual(list(xow.keys()),xow_keys)

            #test with keywords not in arg order
            xow = XenomapperOutputWriter(b"p", b"s", b"s", b"sr",
                             **{x:tempd.name+'/'+x for x in reversed(xow_keys)})
            for key in xow_keys:
                self.assertEqual(xow[key].name.split('/')[-1],
                                 key)
            xow.close()
            tempd.cleanup()

    def test_XenomapperOutputWriter_round_trip(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            test_bam = resource_stream(__name__, 'data/minitest.bam')
            the_bam = bam.FileReader(gzip.open(test_bam))
            out_file = NamedTemporaryFile(delete=False)
            out_file_name = out_file.name
            out_file.close()
            xow = XenomapperOutputWriter(the_bam.raw_header,
                                         the_bam.raw_refs,
                                         the_bam.raw_header,
                                         the_bam.raw_refs,
                                         primary_specific=out_file_name)
            for align in the_bam:
                xow['primary_specific'].write(align)
            xow.close()
            the_bam.close() # this should close test_bam and does below
            test_bam.close() # this should not be needed but get ResourceWarning
            # everything should be closed
            # open original file again
            retest_bam = resource_stream(__name__, 'data/minitest.bam')
            the_bam = bam.FileReader(gzip.open(retest_bam))
            # parse the new temporary file
            new_bam = bam.FileReader(gzip.open(out_file_name))
            self.assertEqual(the_bam.header,new_bam.header[:86])
            try:
                for align1,align2 in zip(new_bam,the_bam):
                    self.assertEqual(align1,align2)
            except StopIteration:
                pass
            new_bam.close()
            the_bam.close()

            with XenomapperOutputWriter(primary_raw_header = b'',
                                        primary_raw_refs = b'',
                                        secondary_raw_header = b'',
                                        secondary_raw_refs = b'',
                                        ) as xow:
                for key in xow.keys():
                    xow[key].write('foo')

    def test_get_mapping_state(self):
        very_negative = -2147483648
        inpt_and_outpt = [
            ((200, 199, 199, 198, very_negative), 'primary_specific'),
            ((200, 200, 199, 198, very_negative), 'primary_multi'),
            ((199, 198, 200, 198, very_negative), 'secondary_specific'),
            ((199, 198, 200, 200, very_negative), 'secondary_multi'),
            ((very_negative, very_negative, very_negative,
              very_negative, very_negative), 'unassigned'),
            ((200, 199, 200, 198, very_negative), 'unresolved'),
            ((200, 199, 199, 199, very_negative), 'primary_specific'),
            ((200, 200, 199, 199, very_negative), 'primary_multi'),
            ((199, 199, 200, 199, very_negative), 'secondary_specific'),
            ((199, 199, 200, 200, very_negative), 'secondary_multi'),
            ((9, 8, 8, 8, 10), 'unassigned'),
            ((200, 200, 200, 200, very_negative), 'unresolved'),
            ((-6, very_negative, very_negative, very_negative, very_negative),
             'primary_specific'),
            ((very_negative, very_negative, -6, very_negative, very_negative),
             'secondary_specific'),
            ((-6, very_negative, -2, very_negative, very_negative),
             'secondary_specific'),
            ((0, very_negative, -2, very_negative, very_negative),
             'primary_specific'),
            ((-2, very_negative, 0, very_negative, very_negative),
             'secondary_specific'),
            ]
        for inpt, outpt in inpt_and_outpt:
            self.assertEqual(get_mapping_state(*inpt), outpt)
        pass

    def test_state_map(self):
        self.assertEqual(state_map(*('primary_specific', 'primary_multi')), 'primary_specific')
        self.assertEqual(state_map(*('primary_specific', 'secondary_specific')), 'primary_specific')
        self.assertEqual(state_map(*('primary_specific', 'secondary_multi')), 'primary_specific')
        self.assertEqual(state_map(*('primary_specific', 'unassigned')), 'primary_specific')
        self.assertEqual(state_map(*('primary_specific', 'unresolved')), 'primary_specific')
        self.assertEqual(state_map(*('primary_specific', None)), 'primary_specific')
        self.assertEqual(state_map(*('primary_multi', 'primary_specific')), 'primary_specific')
        self.assertEqual(state_map(*('primary_multi', 'secondary_specific')), 'secondary_specific')
        self.assertEqual(state_map(*('primary_multi', 'secondary_multi')), 'primary_multi')
        self.assertEqual(state_map(*('primary_multi', 'unassigned')), 'primary_multi')
        self.assertEqual(state_map(*('primary_multi', 'unresolved')), 'primary_multi')
        self.assertEqual(state_map(*('primary_multi', None)), 'primary_multi')
        self.assertEqual(state_map(*('secondary_specific', 'primary_specific')), 'primary_specific')
        self.assertEqual(state_map(*('secondary_specific', 'primary_multi')), 'secondary_specific')
        self.assertEqual(state_map(*('secondary_specific', 'secondary_multi')), 'secondary_specific')
        self.assertEqual(state_map(*('secondary_specific', 'unassigned')), 'secondary_specific')
        self.assertEqual(state_map(*('secondary_specific', 'unresolved')), 'secondary_specific')
        self.assertEqual(state_map(*('secondary_specific', None)), 'secondary_specific')
        self.assertEqual(state_map(*('secondary_multi', 'primary_specific')), 'primary_specific')
        self.assertEqual(state_map(*('secondary_multi', 'primary_multi')), 'primary_multi')
        self.assertEqual(state_map(*('secondary_multi', 'secondary_specific')), 'secondary_specific')
        self.assertEqual(state_map(*('secondary_multi', 'unassigned')), 'secondary_multi')
        self.assertEqual(state_map(*('secondary_multi', 'unresolved')), 'secondary_multi')
        self.assertEqual(state_map(*('secondary_multi', None)), 'secondary_multi')
        self.assertEqual(state_map(*('unassigned', 'primary_specific')), 'primary_specific')
        self.assertEqual(state_map(*('unassigned', 'primary_multi')), 'primary_multi')
        self.assertEqual(state_map(*('unassigned', 'secondary_specific')), 'secondary_specific')
        self.assertEqual(state_map(*('unassigned', 'secondary_multi')), 'secondary_multi')
        self.assertEqual(state_map(*('unassigned', 'unresolved')), 'unresolved')
        self.assertEqual(state_map(*('unassigned', None)), 'unassigned')
        self.assertEqual(state_map(*('unresolved', 'primary_specific')), 'primary_specific')
        self.assertEqual(state_map(*('unresolved', 'primary_multi')), 'primary_multi')
        self.assertEqual(state_map(*('unresolved', 'secondary_specific')), 'secondary_specific')
        self.assertEqual(state_map(*('unresolved', 'secondary_multi')), 'secondary_multi')
        self.assertEqual(state_map(*('unresolved', 'unassigned')), 'unresolved')
        self.assertEqual(state_map(*('unresolved', None)), 'unresolved')
        self.assertEqual(state_map(*(None, 'primary_specific')), 'primary_specific')
        self.assertEqual(state_map(*(None, 'primary_multi')), 'primary_multi')
        self.assertEqual(state_map(*(None, 'secondary_specific')), 'secondary_specific')
        self.assertEqual(state_map(*(None, 'secondary_multi')), 'secondary_multi')
        self.assertEqual(state_map(*(None, 'unassigned')), 'unassigned')
        self.assertEqual(state_map(*(None, 'unresolved')), 'unresolved')
        self.assertRaises(ValueError, state_map, *('foo','bar'))

    def test_conservative_state_map(self):
        self.assertEqual(conservative_state_map(*('primary_specific', 'primary_multi')), 'primary_specific')
        self.assertEqual(conservative_state_map(*('primary_specific', 'secondary_specific')), 'unresolved')
        self.assertEqual(conservative_state_map(*('primary_specific', 'secondary_multi')), 'unresolved')
        self.assertEqual(conservative_state_map(*('primary_specific', 'unassigned')), 'unassigned')
        self.assertEqual(conservative_state_map(*('primary_specific', 'unresolved')), 'unresolved')
        self.assertEqual(conservative_state_map(*('primary_specific', None)), 'primary_specific')
        self.assertEqual(conservative_state_map(*('primary_multi', 'primary_specific')), 'primary_specific')
        self.assertEqual(conservative_state_map(*('primary_multi', 'secondary_specific')), 'unresolved')
        self.assertEqual(conservative_state_map(*('primary_multi', 'secondary_multi')), 'unresolved')
        self.assertEqual(conservative_state_map(*('primary_multi', 'unassigned')), 'unassigned')
        self.assertEqual(conservative_state_map(*('primary_multi', 'unresolved')), 'unresolved')
        self.assertEqual(conservative_state_map(*('primary_multi', None)), 'primary_multi')
        self.assertEqual(conservative_state_map(*('secondary_specific', 'primary_specific')), 'unresolved')
        self.assertEqual(conservative_state_map(*('secondary_specific', 'primary_multi')), 'unresolved')
        self.assertEqual(conservative_state_map(*('secondary_specific', 'secondary_multi')), 'secondary_specific')
        self.assertEqual(conservative_state_map(*('secondary_specific', 'unassigned')), 'unassigned')
        self.assertEqual(conservative_state_map(*('secondary_specific', 'unresolved')), 'unresolved')
        self.assertEqual(conservative_state_map(*('secondary_specific', None)), 'secondary_specific')
        self.assertEqual(conservative_state_map(*('secondary_multi', 'primary_specific')), 'unresolved')
        self.assertEqual(conservative_state_map(*('secondary_multi', 'primary_multi')), 'unresolved')
        self.assertEqual(conservative_state_map(*('secondary_multi', 'secondary_specific')), 'secondary_specific')
        self.assertEqual(conservative_state_map(*('secondary_multi', 'unassigned')), 'unassigned')
        self.assertEqual(conservative_state_map(*('secondary_multi', 'unresolved')), 'unresolved')
        self.assertEqual(conservative_state_map(*('secondary_multi', None)), 'secondary_multi')
        self.assertEqual(conservative_state_map(*('unassigned', 'primary_specific')), 'unassigned')
        self.assertEqual(conservative_state_map(*('unassigned', 'primary_multi')), 'unassigned')
        self.assertEqual(conservative_state_map(*('unassigned', 'secondary_specific')), 'unassigned')
        self.assertEqual(conservative_state_map(*('unassigned', 'secondary_multi')), 'unassigned')
        self.assertEqual(conservative_state_map(*('unassigned', 'unresolved')), 'unassigned')
        self.assertEqual(conservative_state_map(*('unassigned', None)), 'unassigned')
        self.assertEqual(conservative_state_map(*('unresolved', 'primary_specific')), 'unresolved')
        self.assertEqual(conservative_state_map(*('unresolved', 'primary_multi')), 'unresolved')
        self.assertEqual(conservative_state_map(*('unresolved', 'secondary_specific')), 'unresolved')
        self.assertEqual(conservative_state_map(*('unresolved', 'secondary_multi')), 'unresolved')
        self.assertEqual(conservative_state_map(*('unresolved', 'unassigned')), 'unassigned')
        self.assertEqual(conservative_state_map(*('unresolved', None)), 'unresolved')
        self.assertEqual(conservative_state_map(*(None, 'primary_specific')), 'primary_specific')
        self.assertEqual(conservative_state_map(*(None, 'primary_multi')), 'primary_multi')
        self.assertEqual(conservative_state_map(*(None, 'secondary_specific')), 'secondary_specific')
        self.assertEqual(conservative_state_map(*(None, 'secondary_multi')), 'secondary_multi')
        self.assertEqual(conservative_state_map(*(None, 'unassigned')), 'unassigned')
        self.assertEqual(conservative_state_map(*(None, 'unresolved')), 'unresolved')
        self.assertRaises(ValueError, conservative_state_map, *('foo','bar'))

    def test_output_summary(self):
        canned_output = '-'*80 + '\n' + \
                        '\nRead Count Category Summary\n\n\n' + \
                        '|       Category '+' '*36 + '|     Count       |\n' + \
                        '|:' + '-'*50 + ':|:---------------:|\n' + \
                        '|  bar           '+' '*36 + '|            101  |\n' + \
                        '|  foo           '+' '*36 + '|              1  |\n\n'
        canned_output2 = '-'*80 + '\n' + \
                        '\nRead Count Category Summary\n\n\n' + \
                        '|       Category '+' '*36 + '|     Count       |\n' + \
                        '|:' + '-'*50 + ':|:---------------:|\n' + \
                        '|  bar & foo     '+' '*36 + '|            101  |\n' + \
                        '|  foo & bar     '+' '*36 + '|              1  |\n\n'
        test_outfile = io.StringIO()
        output_summary({'foo': 1, 'bar': 101}, outfile=test_outfile)
        self.assertEqual(test_outfile.getvalue(), canned_output)
        test_outfile = io.StringIO()
        output_summary({('foo','bar'): 1, ('bar','foo'): 101},
                       outfile=test_outfile)
        self.assertEqual(test_outfile.getvalue(), canned_output2)
        pass

    def test_AlignbatchFileReader(self):
        test_bam = resource_stream(__name__, 'data/paired_end_testdata_human.bam')
        the_bam = AlignbatchFileReader(gzip.open(test_bam))
        self.assertEqual(next(the_bam), BAMPAIR1)
        for i in range(41):
            next(the_bam)
        self.assertEqual(next(the_bam), BAMPAIR42)
        for alignbatch in the_bam:
            # ensure we read to the end of the file as a test
            pass
        the_bam.close()
        test_bam.close()

    def test_xenomap_states(self):
        self.assertEqual(xenomap_states(BAMPAIR1,BAMPAIR1),
                         ('unresolved', 'unresolved'))
        self.assertRaises(ValueError, xenomap_states,*(BAMPAIR1,BAMPAIR42))
        self.assertEqual(xenomap_states([BAMPAIR1[0],],[BAMPAIR1[0],]),
                         ('unresolved', None))

    def test_cli_help(self):
        #TODO
        #out,err = cli.main("--help")
        #print(out)
        #print(sha256(out.encode()).hexdigest())
        #self.assertEqual('6525535e7f09dedd84a0d7dae93d4327942b2522613cf091daa8e7dc3fdb5f12',
        #                sha256(out.encode()).hexdigest())
        #self.assertEqual('',err)
        pass

    def test_cli(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            output =  io.StringIO()
            prime = resource_filename(__name__,
                                      'data/paired_end_testdata_human.bam')
            second = resource_filename(__name__,
                                       'data/paired_end_testdata_mouse.bam')
            arguments = f"--primary {prime} --secondary {second}"
            self.assertEqual(list(dict(cli.main(arguments,
                                                output
                                                )[1]).values()),
                             [134, 89, 7, 6, 1, 1])
            arguments = (f"--primary {prime} --secondary {second} "
                         "--max")
            self.assertEqual(list(dict(cli.main(arguments,
                                                output
                                                )[1]).values()),
                             [134, 89, 7, 6, 1, 1])
            arguments = (f"--primary {prime} --secondary {second} "
                         "--min_score 190")
            self.assertEqual(list(dict(cli.main(arguments,
                                                output
                                                )[1]).values()),
                             [124, 76, 4, 3, 0, 31])
            arguments = (f"--primary {prime} --secondary {second} "
                         "--conservative")
            self.assertEqual(list(dict(cli.main(arguments,
                                                output
                                                )[1]).values()),
                             [133, 89, 7, 6, 2, 1])
            arguments = (f"--primary {prime} --secondary {second} "
                         "--cigar")
            self.assertEqual(list(dict(cli.main(arguments,
                                                output
                                                )[1]).values()),
                             [141, 95, 0, 0, 1, 1])
            arguments = (f"--primary {prime} --secondary {second} "
                         "--zs")
            self.assertEqual(list(dict(cli.main(arguments,
                                                output
                                                )[1]).values()),
                             [141, 95, 0, 0, 1, 1])

    def test_xenomap(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            primary = AlignbatchFileReader(gzip.open(resource_stream(__name__,'data/paired_end_testdata_human.bam')))
            secondary = AlignbatchFileReader(gzip.open(resource_stream(__name__,'data/paired_end_testdata_mouse.bam')))
            xow = XenomapperOutputWriter(primary.raw_header,
                                         primary.raw_refs,
                                         secondary.raw_header,
                                         secondary.raw_refs,)
            output = xenomap(primary,secondary,xow,get_bamprimary_AS_XS)
            self.assertEqual(list(dict(output[1]).values()),
                             [134, 89, 7, 6, 1, 1])
            self.assertEqual(str(type(output[2])),
                             "<class 'xenomapper2.xenomapper2.XenomapperOutputWriter'>")
            xow.close()
            primary.close()
            secondary.close()

            primary = AlignbatchFileReader(gzip.open(resource_stream(__name__,'data/minitest.sorted.bam')))
            secondary = AlignbatchFileReader(gzip.open(resource_stream(__name__,'data/paired_end_testdata_mouse.bam')))
            self.assertRaises(ValueError,xenomap,*(primary,secondary,xow,get_bamprimary_AS_XS))

            primary = AlignbatchFileReader(gzip.open(resource_stream(__name__,'data/paired_end_testdata_human.bam')))
            secondary = AlignbatchFileReader(gzip.open(resource_stream(__name__,'data/minitest.sorted.bam')))
            self.assertRaises(ValueError,xenomap,*(primary,secondary,xow,get_bamprimary_AS_XS))

            primary = AlignbatchFileReader(gzip.open(resource_stream(__name__,'data/paired_end_testdata_human.bam')))
            secondary = AlignbatchFileReader(gzip.open(resource_stream(__name__,'data/minitest.bam')))
            self.assertRaises(ValueError,xenomap,*(primary,secondary,xow,get_bamprimary_AS_XS))



if __name__ == '__main__':
    unittest.main()
