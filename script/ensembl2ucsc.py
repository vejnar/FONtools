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
import gzip
import sys

import fontools as ft

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Convert names from Ensembl to UCSC (in FASTA, GFF3 and tab).')
    parser.add_argument('-i', '--input', dest='input', action='store', default='-', help='Input.')
    parser.add_argument('-o', '--output', dest='output', action='store', default='-', help='Output.')
    parser.add_argument('-f', '--format', dest='format', action='store', help='Input format (fasta, gff3 or tab).')
    parser.add_argument('-m', '--path_mapping', dest='path_mapping', action='store', help='Path to chromosome name mapping.')
    args = parser.parse_args(argv[1:])

    # Input
    if args.input == '-':
        fin = sys.stdin
        input_format = args.format
    else:
        if args.format is None:
            if args.input.endswith('.fa') or args.input.endswith('.fa.gz'):
                input_format = 'fasta'
            elif args.input.endswith('.gff3') or args.input.endswith('.gff3.gz'):
                input_format = 'gff3'
            elif args.input.endswith('.tab') or args.input.endswith('.tab.gz'):
                input_format = 'tab'
        if args.input.endswith('.gz'):
            fin = gzip.open(args.input, 'rt')
        else:
            fin = open(args.input, 'rt')

    # Name mapping
    name_mapping = None
    if args.path_mapping is not None:
        name_mapping = {}
        with open(args.path_mapping, 'rt') as f:
            for l in f:
                fields = l.strip('\n').split('\t')
                if len(fields) == 2:
                    if fields[1] == '':
                        name_mapping[fields[0]] = None
                    else:
                        name_mapping[fields[0]] = fields[1]

    # Output
    if args.output is None:
        if args.input.endswith('.gz'):
            fout = open(args.input[:-3] + '.ucsc', 'wt')
        else:
            fout = open(args.input + '.ucsc', 'wt')
    elif args.output == '-':
        fout = sys.stdout
    else:
        fout = open(args.output, 'wt')

    # Convert
    if input_format == 'fasta':
        for l in fin:
            if l.startswith('>'):
                fout.write('>' + ft.ensembl.chrom_ensembl2ucsc(l[1:min(l.find(' '), len(l))], name_mapping) + '\n')
            else:
                fout.write(l)
    elif input_format == 'gff3':
        for l in fin:
            if l.startswith('##sequence-region'):
                fields = l.split()
                fields[1] = ft.ensembl.chrom_ensembl2ucsc(fields[1], name_mapping)
                # Write only mappable features
                if fields[1] is not None:
                    fout.write('  '.join(fields)+'\n')
            elif l.startswith('#'):
                fout.write(l)
            else:
                fields = l.split('\t')
                fields[0] = ft.ensembl.chrom_ensembl2ucsc(fields[0], name_mapping)
                # Write only mappable features
                if fields[0] is not None:
                    fout.write('\t'.join(fields))
    elif input_format == 'tab':
        for l in fin:
            if l.startswith('#'):
                fout.write(l)
            else:
                fields = l.split('\t')
                fields[0] = ft.ensembl.chrom_ensembl2ucsc(fields[0], name_mapping)
                # Write only mappable features
                if fields[0] is not None:
                    fout.write('\t'.join(fields))

if __name__ == '__main__':
    sys.exit(main())
