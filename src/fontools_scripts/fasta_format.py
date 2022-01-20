#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import argparse
import os
import shutil
import subprocess
import sys

import pyfaidx

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Sort/Format FASTA file.')
    parser.add_argument('-i', '--input', dest='path_input', action='store', required=True, help='Input.')
    parser.add_argument('-o', '--output', dest='path_output', action='store', help='Output.')
    parser.add_argument('-s', '--sort', dest='sort', action='store_true', help='Sort FASTA sequences.')
    parser.add_argument('-l', '--seq_length', dest='seq_length', action='store', type=int, default=-1, help='Length of sequence line(s)')
    args = parser.parse_args(argv[1:])

    # Prepare unformatted FASTA
    path_output_unft = args.path_output + '_unformatted'
    if args.path_input.endswith('.gz'):
        shutil.copy(args.path_input, path_output_unft + '.gz')
        subprocess.run(['gzip', '-d', path_output_unft + '.gz'], check=True)
    else:
        shutil.copy(args.path_input, path_output_unft)

    # Format FASTA
    gfasta = pyfaidx.Fasta(path_output_unft, as_raw=True, read_ahead=200000)
    with open(args.path_output, 'wt', buffering=1024 * 1024 * 12) as fout:
        chroms = list(gfasta.records)
        if args.sort:
            chroms.sort()
        for chrom in chroms:
            fout.write(f'>{chrom}\n')
            if args.seq_length > 0:
                for i in range(0, len(gfasta[chrom]), args.seq_length):
                    fout.write(gfasta[chrom][i:i+args.seq_length] + '\n')
            else:
                for l in gfasta[chrom]:
                    fout.write(l + '\n')

    # Remove unformatted FASTA and index
    os.remove(gfasta.faidx.indexname)
    os.remove(path_output_unft)

if __name__ == '__main__':
    sys.exit(main())
