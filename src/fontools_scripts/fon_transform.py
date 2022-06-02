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
import itertools
import json
import os
import string
import sys

import fontools as ft
import pyfnutils as pfu
import pyfnutils.log

def split_seq(exons, seq):
    return [seq[s:e] for s, e in exons]

def apply_strand(seqs, strand):
    return [s[::ft.utils.strand2direction(strand)] for s in seqs[::ft.utils.strand2direction(strand)]]

def transform_per_gene(features, union_fn):
    unions = []
    # Sort by genes
    features.sort(key=lambda e: e['gene_stable_id'])
    # Loop over genes
    current_id = ''
    current_union = []
    for a in features:
        if a['gene_stable_id'] != current_id:
            if len(current_union) > 0:
                unions.append(union_fn(current_union))
            current_id = a['gene_stable_id']
            current_union = [a]
        else:
            current_union.append(a)
    if len(current_union) > 0:
        unions.append(union_fn(current_union))
    return unions

def get_annot_union(aset):
    g = aset[0]
    # List of transcripts and proteins
    g['transcript_stable_ids'] = [c['transcript_stable_id'] for c in aset]
    g['protein_stable_ids'] = [c['protein_stable_id'] for c in aset if c['protein_stable_id'] is not None]
    # Union of features
    tlength = 0
    for f in ['exons', 'utr5_exons', 'cds_exons', 'utr3_exons']:
        ivs = ft.utils.flatten([c[f] for c in aset])
        if len(ivs) > 0:
            if f == 'exons':
                seqs = ft.utils.flatten([apply_strand(split_seq(c['exons_on_transcript'], c['seq']), g['strand']) for c in aset if len(c['seq']) > 0])
                g[f], g['seq'] = ft.utils.get_intervals_seq_union(ivs, seqs)
                if g['strand'] == '-':
                    g['seq'] = g['seq'][::-1]
                g[f + '_on_transcript'], _ = ft.utils.get_exons_on_transcript(g[f], g['strand'])
            else:
                g[f], _ = ft.utils.get_intervals_seq_union(ivs, [])
                g[f + '_on_transcript'], tlength = ft.utils.get_exons_on_transcript(g[f], g['strand'], tlength)
        else:
            g[f] = []
            g[f + '_on_transcript'] = []
    # Remove transcript and protein field(s)
    for k in list(g.keys()):
        if k.startswith('transcript_') and k != 'transcript_stable_ids':
            del g[k]
        if k.startswith('protein_') and k != 'protein_stable_ids':
            del g[k]
    return g

def get_flength(features):
    return sum([e-s for s, e in features])

def get_annot_longest(aset):
    mlen = -1
    ilongest = -1
    for i, a in enumerate(aset):
        l = get_flength(a['exons'])
        if l > mlen:
            ilongest = i
            mlen = l
    return aset[ilongest]

def add_intron_to_features(features, exons_key='exons', introns_key='introns'):
    for feat in features:
        exons = feat[exons_key]
        feat[introns_key] = [[exons[i][1], exons[i + 1][0]] for i in range(0, len(exons) - 1)]
    return features

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Transform FON file.')
    parser.add_argument('-f', '--fon', dest='fon', action='store', required=True, help='Input FON file.')
    parser.add_argument('-m', '--method', dest='method', action='store', default='union', help='Fusion method (i.e. union, longest, add_introns).')
    parser.add_argument('-o', '--output', dest='path_output', action='store', required=True, help='Path to output file(s) (comma separated).')
    args = parser.parse_args(argv[1:])

    # Logging
    logger = pfu.log.define_root_logger('fon_' + args.method)
    logger.info('Start')

    # Parsing FON
    logger.info('Open FON')
    fon = json.load(open(args.fon))

    logger.info(f'Transform FON (method:{args.method})')
    if args.method == 'union':
        # Merge using union method
        features = transform_per_gene(fon['features'], get_annot_union)
    elif args.method == 'longest':
        # Merge using longest method
        features = transform_per_gene(fon['features'], get_annot_longest)
        # Sort by transcript
        features.sort(key=lambda e: e['transcript_stable_id'])
    elif args.method == 'add_introns':
        # Add new key "introns" using key "exons"
        features = add_intron_to_features(fon['features'])

    # Path to output
    if args.path_output.find('$') != -1:
        ot = string.Template(args.path_output)
        path_output = [ot.substitute(version='1'), ot.substitute(version='2')]
    else:
        path_output = args.path_output.split(',')
    # Write
    logger.info('FON1 export to '+path_output[0])
    json.dump({'fon_version': 1, 'features': features}, open(path_output[0], 'wt'))

if __name__ == '__main__':
    sys.exit(main())
