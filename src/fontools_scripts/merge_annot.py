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
import json
import os
import string
import sys

import fontools as ft
import pyfnutils as pfu
import pyfnutils.log

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Merge FON/GFF/FASTA files.')
    parser.add_argument('-p', '--species_prefix', dest='species_prefix', action='store', help='Species prefix (comma separated).')
    parser.add_argument('-f', '--input_fasta', dest='input_fasta', action='store', help='Path to input FASTA file(s) (comma separated).')
    parser.add_argument('-o', '--output_fasta', dest='output_fasta', action='store', help='Path to output FASTA file.')
    parser.add_argument('-g', '--input_gff', dest='input_gff', action='store', help='Path to input GTF file(s) (comma separated).')
    parser.add_argument('-a', '--output_gff', dest='output_gff', action='store', help='Path to output GTF file.')
    parser.add_argument('-i', '--input_fon', dest='input_fon', action='store', help='Path to input FON file(s) (template with version, method and biotype) (comma separated).')
    parser.add_argument('-n', '--output_fon', dest='output_fon', action='store', help='Path to output FON file (template with version, method and biotype).')
    parser.add_argument('-y', '--fon_version', dest='fon_version', action='store', default='1', help='FON version (comma separated).')
    parser.add_argument('-m', '--fon_method', dest='fon_method', action='store', default='', help='FON Fusion method: union2gene, longest (comma separated).')
    parser.add_argument('-b', '--fon_biotype', dest='fon_biotype', action='store', default='', help='FON biotype: all, protein_coding (comma separated).')
    args = parser.parse_args(argv[1:])

    # Logging
    logger = pfu.log.define_root_logger('merge')
    logger.info('Start')

    # Parse arguments
    if args.species_prefix is not None:
        species_prefix = args.species_prefix.strip().split(',')
    else:
        species_prefix = []
    if args.input_fasta is None:
        input_fastas = []
    else:
        input_fastas = args.input_fasta.strip().split(',')
    if args.input_gff is None:
        input_gffs = []
    else:
        input_gffs = args.input_gff.strip().split(',')
    input_fons = []
    if args.input_fon is not None:
        for p in args.input_fon.strip().split(','):
            if p.find('$') != -1:
                input_fons.append(string.Template(p))
            else:
                input_fons.append(p)
    output_fon = None
    if args.output_fon is not None:
        if args.output_fon.find('$') != -1:
            output_fon = string.Template(args.output_fon)
        else:
            output_fon = args.output_fon
    if args.fon_version is not None:
        fon_versions = args.fon_version.strip().split(',')
    if args.fon_method is not None:
        fon_methods = args.fon_method.strip().split(',')
    if args.fon_biotype is not None:
        fon_biotypes = args.fon_biotype.strip().split(',')
    if len(species_prefix) == 0:
        species_prefix = [''] * max(len(input_fastas), len(input_gffs), len(input_fons))

    # Checks
    if len(input_fastas) > 0 and len(species_prefix) != len(input_fastas):
        logger.error('Input FASTA and species prefix number differ')
        return 1
    if len(input_gffs) > 0 and len(species_prefix) != len(input_gffs):
        logger.error('Input GFF and species prefix number differ')
        return 1
    if len(input_fons) > 0 and len(species_prefix) != len(input_fons):
        logger.error('Input FON and species prefix number differ')
        return 1
    if args.output_fasta is not None and os.path.exists(args.output_fasta):
        logger.error('Output FASTA file already exists')
        return 1
    if args.output_gff is not None and os.path.exists(args.output_gff):
        logger.error('Output GFF file already exists')
        return 1

    # FASTA
    if args.output_fasta is not None:
        logger.info('FASTA: ' + args.output_fasta)
        new_names = set()
        for i, p in enumerate(input_fastas):
            logger.info('Merging: ' + p)
            with open(p) as f, open(args.output_fasta, 'at') as fout:
                for l in f:
                    if l.startswith('>'):
                        new_name = species_prefix[i] + l[1:].lstrip()
                        if new_name in new_names:
                            logger.error('Name collapsing ' + new_name)
                        else:
                            new_names.add(new_name)
                        fout.write('>' + new_name)
                    else:
                        fout.write(l)

    # GFF
    if args.output_gff is not None:
        logger.info('GFF: ' + args.output_gff)
        for i, p in enumerate(input_gffs):
            logger.info('Merging: ' + p)
            with open(p) as f, open(args.output_gff, 'at') as fout:
                for l in f:
                    if l.startswith('#'):
                        fout.write(l)
                    else:
                        fout.write(species_prefix[i] + l)

    # FON
    if output_fon is not None:
        for fon_version in fon_versions:
            for fon_method in fon_methods:
                for fon_biotype in fon_biotypes:
                    # Input
                    fins = []
                    merged_features = None
                    for i, p in enumerate(input_fons):
                        if isinstance(p, string.Template):
                            fin = p.substitute(version=fon_version, method=fon_method, biotype=fon_biotype)
                        else:
                            fin = p
                        features = json.load(open(fin))['features']
                        if fon_version == '1':
                            for j in range(len(features)):
                                features[j]['chrom'] = species_prefix[i] + features[j]['chrom'] 
                            if merged_features is None:
                                merged_features = features
                            else:
                                merged_features.extend(features)
                        fins.append(fin)
                    # Output
                    if isinstance(output_fon, string.Template):
                        fout = output_fon.substitute(version=fon_version, method=fon_method, biotype=fon_biotype)
                    else:
                        fout = output_fon
                    logger.info('FON: ' + fout)
                    for f in fins:
                        logger.info('Merging: ' + f)
                    if fon_version == '1':
                        json.dump({'fon_version': 1, 'features': merged_features}, open(fout, 'wt'))

if __name__ == '__main__':
    sys.exit(main())
