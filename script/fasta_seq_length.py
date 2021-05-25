#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2021 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import argparse
import os
import sys

import pyfaidx

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Create tab file with sequence(s) length from FASTA file.')
    parser.add_argument('-i', '--input', dest='path_input', action='store', required=True, help='Input.')
    parser.add_argument('-o', '--output', dest='path_output', action='store', help='Output.')
    args = parser.parse_args(argv[1:])

    gfasta = pyfaidx.Fasta(args.path_input, as_raw=True)
    with open(args.path_output, 'wt') as fout:
        for chrom in sorted(gfasta.keys()):
            fout.write(f'{chrom}\t{len(gfasta[chrom])}\n')

if __name__ == '__main__':
    sys.exit(main())
