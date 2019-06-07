#!/usr/bin/env python3
# -*- coding: utf-8 -*-


if __name__ == "__main__" :
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="""sSNV/gSNV single cell recalling""")

    argparser.add_argument('bamfiles', metavar='bamfile', type=str, nargs='+')

    argparser.add_argument('-variants', help='sSNV, gSNV list file: chromsome, sSNVlocation, sSNValt sSNVref gSNVpos gSNValt gSNVref')

    args = argparser.parse_args()
    
