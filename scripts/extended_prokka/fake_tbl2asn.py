#!/usr/bin/env python
"""A really fake script to simulate the tbl2asn api but returns empty 
output files. 

tbl2asn -V b -a s -N 1 -y 'Annotated using prokka 1.7 from http://www.vicbioinformatics.com' -Z prokka_output/PROKKA_11222014.err -i prokka_output/PROKKA_11222014.fsa 2> /dev/null
"""
import argparse
import os

def main(args):
    if args.version:
        print 'tbl2asn 24.2 arguments:'
        return
    with open(args.Z, 'a'):
        os.utime(args.Z, None)
    with open(args.i, 'a'):
        os.utime(args.i, None)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-V', help='This is a really b arg')
    parser.add_argument('-a', help='This should be s')
    parser.add_argument('-N', help='Give me one')
    parser.add_argument('-y', help="Long text goes here")
    parser.add_argument('-Z', help="Error output")
    parser.add_argument('-i', help="Fsa file")
    parser.add_argument('-', dest='version', action='store_true', help='Print version')
    args = parser.parse_args()
    main(args)
