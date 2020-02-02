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

from pkg_resources import resource_stream

from xenomapper2.xenomapper2 import *

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2020 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield", ]
__license__ = "GPLv3"
__version__ = "2.0.0"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Beta"

BAMPAIR1 = [
    b'\x12\x01\x00\x00\x0c\x00\x00\x00eB\xf0\x07(&\n2\x02\x00S\x00d\x00\x00\x00\x0c\x00\x00\x00\x1eB\xf0\x07S\xff\xff\xffHWI-ST960:96:COTO3ACXX:3:1101:1220:2089\x000\x06\x00\x00\x14\x00\x00\x00HB\x18"$H\x14\x82"\x14(\x12\x88\x14AD(AD!D\x14\x11\x82B\x88A\x12"\x84AD!\x11D\x88B\x14\x84\x14"AA\x82\x12\x12!(\x12\x1f#"##!\x1f####""!!\x1e##$$$##$"#"%%\'\'\'\'\')((&\')))))))(\'(()((\'#!))))))))((()))(&\'\'))))())()(&))(\'\'%%\'%####\x1c\x10\x02ASC\xc6XSC~XNC\x00XMC\x00XOC\x00XGC\x00NMC\x00MDZ99\x00YSC\xbdYTZCP\x00',
    b'\x15\x01\x00\x00\x0c\x00\x00\x00\x1eB\xf0\x07(&\n2\x02\x00\xa3\x00d\x00\x00\x00\x0c\x00\x00\x00eB\xf0\x07\xad\x00\x00\x00HWI-ST960:96:COTO3ACXX:3:1101:1220:2089\x004\x00\x00\x00\x10\x06\x00\x00\xff\xf8\x82\x11\x11"\x14"\x84\x12!\x12\x18D\x84\x11\x12"\x84\x82\x82\x81(\x11\x11\x18\x12\x11\x11\x88\x14(DHHD\x84HB\x18"\x84H\x14\x82"\x14(\x12\x88\x02\x02\x02\x13\x13 #%\'&\'%\')()))))))))))))$&\'())))))()))()))()\'(\'\'(((\'))))()(()))"\'$$%### !##$#\x1e\x0c\x0c\x19\x1d ""###"##"\x1f""#ASC\xbdXSC2XNC\x00XMC\x01XOC\x00XGC\x00NMC\x01MDZ79C17\x00YSC\xc6YTZCP\x00']
BAMPAIR42 = [
    b'\x1a\x01\x00\x00\t\x00\x00\x00\xd1B\xd2\x07(,\x921\x04\x00c\x00d\x00\x00\x00\t\x00\x00\x00%D\xd2\x07\xb9\x01\x00\x00HWI-ST960:96:COTO3ACXX:3:1101:1591:2148\x00\x14\x00\x00\x00\xe0\x02\x00\x00\x11\x00\x00\x00@\x03\x00\x00\xf1\x88"\x88D(\x11\x14A\x82\x14\x81AD\x82\x81\x81\x81\x81\x81\x81\x81\x81\x81\x88\x88\x88\x88\x88\x82"((\x88\x88\x11\x11\x82\x12A\x84\x81\x18\x12\x18\x18\x12\x14!\x11\x02\x10\x19#### %\'\'\'\'((&(!&&$& "\'\'\'&\n \'&\'((((((((\x18#&(%(((\x1f##&(((\x1d$#\x16\r\x1d\x1d"\x1e" \x1f""\x1e"" "\x1d\x1d \x1b\x10\x1e\x1e\x1e#"#"\x1d""#""\x19\x1d\x19 """ASC\xbcXSC<XNC\x00XMC\x00XOC\x01XGC\x01NMC\x01MDZ98\x00YSC\xc8YTZCP\x00',
    b'\x0f\x01\x00\x00\t\x00\x00\x00%D\xd2\x07(,\x921\x01\x00\x93\x00d\x00\x00\x00\t\x00\x00\x00\xd1B\xd2\x07G\xfe\xff\xffHWI-ST960:96:COTO3ACXX:3:1101:1591:2148\x00@\x06\x00\x00D\x82\x14AH\x88AA"\x14"\x84\x12!\x12\x18D\x84\x11\x12"$\x82\x88\x81(\x14\x11\x12\x18\x11\x11\x11\x88\x14(DB\x18D\x84B\x12\x12B(H\x11\x82"\x18\x1e\x1f\x19\x19\x13\x07\x19\x19\x19\x1d\x1a\x1d"\x19\x0b\x0b\x0b\x1e\x1b\x12\x1e\x1a\x1a\x1a \x1f\x1a\x1a"\x1d\x1d\x1f\x19\x1e\x16\x14\x06\x13!\x17\x0e\x1a\x16\x16\x1b!\x1c\x13\'&&&$%!\x17\'%\x1a%\'\'#(\'#\x18\x1f\x1b&&\'&$"\x18&&(&\'\x1f%$$\x1b#\x16#\x1e\'# \x1a##\x1e\x1f\x1eASC\xc8XSC\x8dXNC\x00XMC\x00XOC\x00XGC\x00NMC\x00MDZ100\x00YSC\xbcYTZCP\x00']


class test_main(unittest.TestCase):
    def setUp(self):
        pass

    # def test_process_headers(self):
    #     test_primary_specific_outfile = io.StringIO()
    #     test_secondary_specific_outfile = io.StringIO()
    #     test_primary_multi_outfile = io.StringIO()
    #     test_secondary_multi_outfile = io.StringIO()
    #     test_unassigned_outfile = io.StringIO()
    #     test_unresolved_outfile = io.StringIO()
    #     sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_human.sam'))
    #     sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_mouse.sam'))
    #     process_headers(sam1,sam2,
    #                  primary_specific=test_primary_specific_outfile,
    #                  secondary_specific=test_secondary_specific_outfile,
    #                  primary_multi=test_primary_multi_outfile,
    #                  secondary_multi=test_secondary_multi_outfile,
    #                  unresolved=test_unresolved_outfile,
    #                  unassigned=test_unassigned_outfile,
    #                  )
    #     self.assertEqual(len(test_primary_specific_outfile.getvalue()),695)
    #     self.assertEqual(len(test_secondary_specific_outfile.getvalue()),629)
    #     self.assertEqual(len(test_primary_multi_outfile.getvalue()),708)
    #     self.assertEqual(len(test_secondary_multi_outfile.getvalue()),642)
    #     self.assertEqual(len(test_unassigned_outfile.getvalue()),705)
    #     self.assertEqual(len(test_unresolved_outfile.getvalue()),705)
    #     sam1.close()
    #     sam2.close()
    #     pass
    #
    # def test_consistent_output_SE(self):
    #     test_primary_specific_outfile = io.StringIO()
    #     test_secondary_specific_outfile = io.StringIO()
    #     test_primary_multi_outfile = io.StringIO()
    #     test_secondary_multi_outfile = io.StringIO()
    #     test_unassigned_outfile = io.StringIO()
    #     test_unresolved_outfile = io.StringIO()
    #     sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/test_human_in.sam'))
    #     sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/test_mouse_in.sam'))
    #     process_headers(sam1,sam2,
    #                      primary_specific=test_primary_specific_outfile,
    #                      secondary_specific=test_secondary_specific_outfile,
    #                      primary_multi=test_primary_multi_outfile,
    #                      secondary_multi=test_secondary_multi_outfile,
    #                      unresolved=test_unresolved_outfile,
    #                      unassigned=test_unassigned_outfile,
    #                      )
    #     cat_counts = main_single_end(getReadPairs(sam1,sam2),
    #                      primary_specific=test_primary_specific_outfile,
    #                      secondary_specific=test_secondary_specific_outfile,
    #                      primary_multi=test_primary_multi_outfile,
    #                      secondary_multi=test_secondary_multi_outfile,
    #                      unresolved=test_unresolved_outfile,
    #                      unassigned=test_unassigned_outfile,
    #                      )
    #     self.assertEqual(cat_counts['primary_specific'],
    #                      len(test_primary_specific_outfile.getvalue().split('\n'))-4) #29 lines of header in this file
    #     self.assertEqual(cat_counts['primary_multi'],
    #                      len(test_primary_multi_outfile.getvalue().split('\n'))-4) #29 lines of header in this file
    #     self.assertEqual(cat_counts['secondary_specific'],
    #                      len(test_secondary_specific_outfile.getvalue().split('\n'))-4) #26 lines of header in this file
    #     self.assertEqual(cat_counts['secondary_multi'],
    #                      len(test_secondary_multi_outfile.getvalue().split('\n'))-4) #26 lines of header in this file
    #     self.assertEqual(cat_counts['unassigned'],
    #                      len(test_unassigned_outfile.getvalue().split('\n'))-4) #26 lines of header in this file
    #     self.assertEqual(cat_counts['unassigned'],
    #                      len(test_unassigned_outfile.getvalue().split('\n'))-4) #26 lines of header in this file
    #     self.assertEqual(hashlib.sha224(test_primary_specific_outfile.getvalue().encode('latin-1')).hexdigest(),'381325b12dd9a9cd3afdd72eeb16b23cc92ddd16f675bb21bb21e08e')
    #     sam1.close()
    #     sam2.close()
    #     pass
    #
    # def test_consistent_output_PE(self):
    #     test_primary_specific_outfile = io.StringIO()
    #     test_secondary_specific_outfile = io.StringIO()
    #     test_primary_multi_outfile = io.StringIO()
    #     test_secondary_multi_outfile = io.StringIO()
    #     test_unassigned_outfile = io.StringIO()
    #     test_unresolved_outfile = io.StringIO()
    #     sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_human.sam'))
    #     sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_mouse.sam'))
    #     process_headers(sam1,sam2,primary_specific=test_primary_specific_outfile, secondary_specific=test_secondary_specific_outfile)
    #     cat_counts = main_paired_end(getReadPairs(sam1,sam2),
    #                                  primary_specific=test_primary_specific_outfile,
    #                                  secondary_specific=test_secondary_specific_outfile,
    #                                  primary_multi=test_primary_multi_outfile,
    #                                  secondary_multi=test_secondary_multi_outfile,
    #                                  unresolved=test_unresolved_outfile,
    #                                  unassigned=test_unassigned_outfile,
    #                                  )
    #     self.assertEqual(sum([cat_counts[x] for x in cat_counts if 'primary_specific' in x])*2,
    #                      len(test_primary_specific_outfile.getvalue().split('\n'))-30) #29 lines of header in this file
    #     self.assertEqual(sum([cat_counts[x] for x in cat_counts if 'secondary_specific' in x and not 'primary_specific' in x])*2,
    #                      len(test_secondary_specific_outfile.getvalue().split('\n'))-27) #26 lines of header in this file
    #     self.assertEqual(sum([cat_counts[x] for x in cat_counts if 'primary_multi' in x and not 'primary_specific' in x and not 'secondary_specific' in x])*2,
    #                     len(test_primary_multi_outfile.getvalue().split('\n'))-1)
    #     self.assertEqual(sum([cat_counts[x] for x in cat_counts if 'secondary_multi' in x \
    #                     and not 'primary_multi' in x and not 'primary_specific' in x and not 'secondary_specific' in x])*2,
    #                     len(test_secondary_multi_outfile.getvalue().split('\n'))-1)
    #     self.assertEqual(hashlib.sha224(test_primary_specific_outfile.getvalue().encode('latin-1')).hexdigest(),'64c0e24bf141c5aa3bb0993c73b34cdfe630a504ac424843f746918d')
    #     sam1.close()
    #     sam2.close()
    #     pass
    #
    # def test_consistent_output_conservative_PE(self):
    #     test_primary_specific_outfile = io.StringIO()
    #     test_secondary_specific_outfile = io.StringIO()
    #     test_primary_multi_outfile = io.StringIO()
    #     test_secondary_multi_outfile = io.StringIO()
    #     test_unassigned_outfile = io.StringIO()
    #     test_unresolved_outfile = io.StringIO()
    #     sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_human.sam'))
    #     sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_mouse.sam'))
    #     process_headers(sam1,sam2,primary_specific=test_primary_specific_outfile, secondary_specific=test_secondary_specific_outfile)
    #     cat_counts = conservative_main_paired_end(getReadPairs(sam1,sam2),
    #                                              primary_specific=test_primary_specific_outfile,
    #                                              secondary_specific=test_secondary_specific_outfile,
    #                                              primary_multi=test_primary_multi_outfile,
    #                                              secondary_multi=test_secondary_multi_outfile,
    #                                              unresolved=test_unresolved_outfile,
    #                                              unassigned=test_unassigned_outfile,
    #                                              )
    #     self.assertEqual(cat_counts[('primary_specific', 'secondary_specific')]*4 +\
    #                      cat_counts[('unresolved', 'unresolved')]*4, len(test_unresolved_outfile.getvalue().split('\n'))-1) #this test is not exhastive. Only states in test data
    #     self.assertEqual(cat_counts[('primary_multi', 'primary_multi')]*2, len(test_primary_multi_outfile.getvalue().split('\n'))-1)
    #     self.assertEqual(cat_counts[('secondary_multi', 'secondary_multi')]*2, len(test_secondary_multi_outfile.getvalue().split('\n'))-1)
    #     self.assertEqual(cat_counts[('primary_specific', 'primary_specific')]*2 + cat_counts[('primary_specific', 'primary_multi')]*2 + cat_counts[('primary_multi', 'primary_specific')]*2,
    #                      len(test_primary_specific_outfile.getvalue().split('\n'))-30) #29 lines of header in this file
    #     self.assertEqual(cat_counts[('secondary_specific', 'secondary_specific')]*2 + cat_counts[('secondary_specific', 'secondary_multi')]*2 + cat_counts[('secondary_multi', 'secondary_specific')]*2,
    #                      len(test_secondary_specific_outfile.getvalue().split('\n'))-27)  #26 lines of header in this file
    #     self.assertEqual(cat_counts[('unassigned', 'unassigned')]*2, len(test_unassigned_outfile.getvalue().split('\n'))-1)
    #
    #     self.assertEqual(hashlib.sha224(test_primary_specific_outfile.getvalue().encode('latin-1')).hexdigest(),'c4de3de755092c8f9ff1eb2cd360a502d74ebd4c1e65ed282515ed3e')
    #     sam1.close()
    #     sam2.close()
    #     pass

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
        self.assertEqual(get_max_AS_XS([]), (float('-inf'), float('-inf')))

    def test_get_cigar_based_score(self):
        input_and_output = [(('*', None), float('-inf')),
                            (('50M',0),0),
                            (('1S49M',0),-2),
                            (('50M', 2), -12),
                            (('10M1I39M', 0), -8),
                            (('10M2D38M', 0), -11),
                            (('10M1I10M1D28M', 0), -16),
                            (('10M1234N40M', 0), 0),
                            ]
        for (inpt, outpt) in input_and_output:
            self.assertEqual(get_cigarbased_score(*inpt), outpt)

    def test_XenomapperOutputWriter(self):
        xow = XenomapperOutputWriter("p", "s")
        self.assertEqual([repr(x) for x in list(xow._fileobjects.values())],
                         [repr(DummyFile()), ] * 6)
        # need to set up a tempfile dir
        xow = XenomapperOutputWriter("p", "s", basename='foo')
        self.assertEqual([x._handle.name for x in xow._fileobjects.values()],
                         ['foo_primary_specific', 'foo_primary_multi',
                          'foo_secondary_specific', 'foo_secondary_multi',
                          'foo_unresolved', 'foo_unassigned'])
        xow.close()

    def test_get_mapping_state(self):
        inpt_and_outpt = [
            ((200, 199, 199, 198, float('-inf')), 'primary_specific'),
            ((200, 200, 199, 198, float('-inf')), 'primary_multi'),
            ((199, 198, 200, 198, float('-inf')), 'secondary_specific'),
            ((199, 198, 200, 200, float('-inf')), 'secondary_multi'),
            ((float('-inf'), float('-inf'), float('-inf'), float('-inf'), float('-inf')), 'unassigned'),
            ((200, 199, 200, 198, float('-inf')), 'unresolved'),
            ((200, 199, 199, 199, float('-inf')), 'primary_specific'),
            ((200, 200, 199, 199, float('-inf')), 'primary_multi'),
            ((199, 199, 200, 199, float('-inf')), 'secondary_specific'),
            ((199, 199, 200, 200, float('-inf')), 'secondary_multi'),
            ((9, 8, 8, 8, 10), 'unassigned'),
            ((200, 200, 200, 200, float('-inf')), 'unresolved'),
            ((-6, float('-inf'), float('-inf'), float('-inf'), float('-inf')), 'primary_specific'),
            ((float('-inf'), float('-inf'), -6, float('-inf'), float('-inf')), 'secondary_specific'),
            ((-6, float('-inf'), -2, float('-inf'), float('-inf')), 'secondary_specific'),
            ((0, float('-inf'), -2, float('-inf'), float('-inf')), 'primary_specific'),
            ((-2, float('-inf'), 0, float('-inf'), float('-inf')), 'secondary_specific'),

        ]

        for inpt, outpt in inpt_and_outpt:
            self.assertEqual(get_mapping_state(*inpt), outpt)
        pass

    # def test_get_tag(self):
    #     inpt_and_outpt = [
    #                     ((['HWI-ST960:63:D0CYJACXX:4:1101:21264:2228', '4', '*', '0', '0', '*', '*', '0', '0', 'TGGTAGTATTGGTTATGGTTCATTGTCCGGAGAGTATATTGTTGAAGAGG', 'BBCBDFDDHHHGFHHIIIIIJIJJJIGJJJGIAF:CFEGHGGHEEEG@HI', 'YT:Z:UU'],'AS'),
    #                     float('-inf')),
    #                     ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:101', 'XS:i:99'],'AS'),101),
    #                     ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:i:99'],'XS'),99),
    #                     ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:i:99'],'NM'),0),
    #                     ]
    #     for inpt, outpt in inpt_and_outpt:
    #         self.assertEqual(get_tag(*inpt),outpt)
    #     pass
    #
    # def test_get_tag_with_ZS_as_XS(self):
    #     inpt_and_outpt = [
    #                     ((['HWI-ST960:63:D0CYJACXX:4:1101:21264:2228', '4', '*', '0', '0', '*', '*', '0', '0', 'TGGTAGTATTGGTTATGGTTCATTGTCCGGAGAGTATATTGTTGAAGAGG', 'BBCBDFDDHHHGFHHIIIIIJIJJJIGJJJGIAF:CFEGHGGHEEEG@HI', 'YT:Z:UU'],'AS'),
    #                     float('-inf')),
    #                     ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:101', 'XS:A:+', 'ZS:i:99'],'AS'),101),
    #                     ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:A:+', 'ZS:i:99'],'XS'),99),
    #                     ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:A:+', 'ZS:i:99'],'NM'),0),
    #                     ]
    #     for inpt, outpt in inpt_and_outpt:
    #         self.assertEqual(get_tag_with_ZS_as_XS(*inpt),outpt)
    #     pass
    #


def test_output_summary(self):
    test_outfile = io.StringIO()
    canned_output = '--------------------------------------------------------------------------------\n' + \
                    'Read Count Category Summary\n\n' + \
                    '|       Category                                     |     Count       |\n' + \
                    '|:--------------------------------------------------:|:---------------:|\n' + \
                    '|  bar                                               |            101  |\n' + \
                    '|  foo                                               |              1  |\n\n'
    output_summary({'foo': 1, 'bar': 101}, outfile=test_outfile)
    self.assertEqual(test_outfile.getvalue(), canned_output)
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


if __name__ == '__main__':
    unittest.main()
