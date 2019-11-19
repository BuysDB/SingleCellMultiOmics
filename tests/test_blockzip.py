#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
from singlecellmultiomics.utils import BlockZip
import os
"""
These tests check if the BlockZip module is working correctly
"""

class Test_BlockZip(unittest.TestCase):

    def test_read_write(self):
        zip_path = './data/test.bgzf'
        with BlockZip(zip_path,'w') as f:
            f.write('chr1',100,False,'yes')
            f.write('chr1',101,False,'yes')
            f.write('chr2',2101,True,'yes')
            f.write('chrX',0,False,'X')

        with BlockZip(zip_path,'r') as f:
            self.assertEqual( f[ ('chr1',100,False) ], 'yes' )
            self.assertEqual( f[ ('chr1',0,False) ], None)
            self.assertEqual( f[ ('chrX',0,False) ], 'X' )
            self.assertEqual( f[ ('chr2',2101,True) ], 'yes' )
            self.assertEqual( f[ ('chr1',0,True) ], None)

        os.remove(zip_path)


if __name__ == '__main__':
    unittest.main()
