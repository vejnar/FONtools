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
import json
import sys

import pyfaidx

import fontools as ft
import pyfnutils as pfu
import pyfnutils.log

def group_coords_by_chrom(feats, key, extension=0, exterior_extension=0):
    by_chrom = {}
    esup = max(0, exterior_extension - extension)
    for feat in feats:
        chrom = feat['chrom']
        coords = feat[key]
        if chrom not in by_chrom:
            by_chrom[chrom] = []
        if extension > 0:
            for ics, cs in enumerate(coords):
                coords[ics] = [max(0, cs[0]-extension), cs[1]+extension]
        if esup > 0:
            coords[0][0] = max(0, coords[0][0]-esup)
            coords[-1][1] = coords[-1][1] + esup
        by_chrom[chrom].extend(coords)
    return by_chrom

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Mask part(s) of sequence (FASTA).')
    parser.add_argument('-i', '--input_fon', dest='input_fon', action='store', required=True, help='Path to input FON1 file.')
    parser.add_argument('-f', '--input_fasta', dest='input_fasta', action='store', help='Path to input FASTA file.')
    parser.add_argument('-o', '--output_fasta', dest='output_fasta', action='store', help='Path to output FASTA file.')
    parser.add_argument('-k', '--key', dest='key', action='store', default='exons', help='Key (default: exons).')
    parser.add_argument('-e', '--extension', dest='extension', action='store', type=int, default=0, help='Extension.')
    parser.add_argument('-x', '--exterior_extension', dest='exterior_extension', action='store', type=int, default=-1, help='Extension for exterior boundaries of features.')
    parser.add_argument('-l', '--seq_length', dest='seq_length', action='store', type=int, default=1000, help='Length of sequence line(s)')
    parser.add_argument('-r', '--inverse', dest='inverse', action='store_true', default=False, help='Inverse mask (parts in input FON are masked)')
    args = parser.parse_args(argv[1:])

    # Logging
    logger = pfu.log.define_root_logger('mask')
    logger.info('Start')

    # If exterior_extension is unspecified then it's equal to extension
    if args.exterior_extension == -1:
        args.exterior_extension = args.extension

    # Parsing FON
    logger.info('Open FON')
    fon = json.load(open(args.input_fon))

    # Group by chromosome
    coords_by_chrom = group_coords_by_chrom(fon['features'], args.key, args.extension, args.exterior_extension)

    # Mask FASTA
    gfasta = pyfaidx.Fasta(args.input_fasta, as_raw=True, read_ahead=200000)
    logger.info(f'Writing masked FASTA {args.output_fasta}')
    nunmsk = 0
    ntot = 0
    with open(args.output_fasta, 'wb', buffering=1024 * 1024 * 12) as fout:
        for name, seq in gfasta.records.items():
            if args.inverse:
                mseq = bytearray(seq[:], 'ascii')
                if name in coords_by_chrom:
                    for cs in coords_by_chrom[name]:
                        mseq[cs[0]:cs[1]] = b'N' * (cs[1] - cs[0])
            else:
                mseq = bytearray(b'N') * len(seq)
                if name in coords_by_chrom:
                    for cs in coords_by_chrom[name]:
                        mseq[cs[0]:cs[1]] = bytearray(seq[cs[0]:cs[1]], 'ascii')
            fout.write(bytes(f'>{name}\n', 'utf8'))
            for i in range(0, len(mseq), args.seq_length):
                fout.write(mseq[i:i+args.seq_length] + b'\n')
            nunmsk += len(mseq) - mseq.count(b'N')
            ntot += len(mseq)
    logger.info(f'Unmasked sequence: {nunmsk}nt ({nunmsk/ntot*100.:.1f}%)')

if __name__ == '__main__':
    sys.exit(main())
