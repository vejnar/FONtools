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
import collections
import gzip
import json
import os
import string
import sys
import tempfile

import pyfaidx

import fontools as ft
import pyfnutils as pfu
import pyfnutils.log

def add_seq(transcripts, path_fasta, cdna, exclude_no_seq):
    new_transcripts = []
    count_noseq = collections.Counter()
    with tempfile.TemporaryDirectory() as tmpd:
        seqs = []
        for ipf, pf in enumerate(path_fasta):
            if pf.endswith('.gz'):
                path_in = os.path.join(tmpd, f'{ipf}.fa')
                with gzip.open(pf, 'rt') as f, open(path_in, 'wt') as fout:
                    fout.write(f.read())
            else:
                path_in = pf
            seqs.append(pyfaidx.Fasta(path_in, as_raw=True))
        if cdna:
            for t in transcripts:
                name = t['transcript_stable_id']
                if 'transcript_version' in t and t['transcript_version'] is not None:
                    name += '.' + t['transcript_version']
                found = False
                for seq in seqs:
                    if name in seq:
                        t['seq'] = seq[name][:]
                        found = True
                if found:
                    new_transcripts.append(t)
                else:
                    count_noseq[t['transcript_biotype']] += 1
                    if not exclude_no_seq:
                        t['seq'] = None
                        new_transcripts.append(t)
        else:
            for t in transcripts:
                name = t['chrom']
                found = False
                for seq in seqs:
                    if name in seq:
                        eseq = [seq[name][s:e] for s, e in t['exons']]
                        if t['strand'] == '+':
                            t['seq'] = ''.join(eseq)
                        elif t['strand'] == '-':
                            t['seq'] = ''.join([pyfaidx.complement(s[::-1]) for s in eseq[::-1]])
                        found = True
                if found:
                    new_transcripts.append(t)
                else:
                    count_noseq[t['transcript_biotype']] += 1
                    if not exclude_no_seq:
                        t['seq'] = None
                        new_transcripts.append(t)
    return new_transcripts, count_noseq

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Import annotations to FON.')
    parser.add_argument('-a', '--annotation', dest='path_annot', action='store', required=True, help='Path to annotation file(s) (comma separated).')
    parser.add_argument('-o', '--output', dest='path_output', action='store', required=True, help='Path to output file (comma separated).')
    parser.add_argument('-t', '--output_format', dest='output_format', action='store', default='fon', help='Output format: fon or fon1.')
    parser.add_argument('-e', '--output_feature', dest='output_feature', action='store', default='transcript_stable_id', help='Output features.')
    parser.add_argument('-d', '--output_coord', dest='output_coord', action='store', default='exons', help='Output coordinates.')
    parser.add_argument('-s', '--data_source', dest='data_source', action='store', help='Data source: ensembl, or xenbase.')
    parser.add_argument('-b', '--biotype', dest='biotype', action='store', default='all', help='Biotype: all, protein_coding (comma separated).')
    parser.add_argument('-f', '--fasta', dest='path_fasta', action='append', help='Path to FASTA file(s).')
    parser.add_argument('-c', '--cdna', dest='cdna', action='store_true', help='FASTA file containing cDNA.')
    parser.add_argument('-u', '--ucsc_names', dest='convert_ucsc_names', action='store_true', help='Convert to UCSC chromosome/contig names.')
    parser.add_argument('-m', '--path_mapping', dest='path_mapping', action='store', help='Path to chromosome name mapping.')
    parser.add_argument('-l', '--filter', dest='source_filter', action='append', help='Filter source (ex: ensembl, havana).')
    parser.add_argument('-x', '--exclude_no_seq', dest='exclude_no_seq', action='store_true', help='Exclude transcript without cDNA sequence.')
    parser.add_argument('--bed_name_as_id', dest='bed_name_as_id', action='store_true', help='BED: Use name (BED 4th column) as ID.')
    parser.add_argument('--bed_feature_id', dest='bed_feature_id', action='store', default='transcript_stable_id', help='BED: Feature ID (default:transcript_stable_id).')
    parser.add_argument('--bed_interval_name', dest='bed_interval_name', action='store', default='exons', help='BED: Interval name (default:exons).')
    parser.add_argument('--table_transcript', dest='path_table_transcript', action='store', help='Path to table transcript.txt.gz'),
    parser.add_argument('--table_object_xref', dest='path_table_object_xref', action='store', help='Path to table object_xref.txt.gz'),
    parser.add_argument('--table_xref', dest='path_table_xref', action='store', help='Path to table xref.txt.gz'),
    parser.add_argument('--table_ontology', dest='path_table_ontology', action='store', help='Path to table ontology.txt.gz'),
    parser.add_argument('--table_term', dest='path_table_term', action='store', help='Path to table term.txt.gz')
    args = parser.parse_args(argv[1:])

    # Logging
    logger = pfu.log.define_root_logger('fon_import')
    logger.info('Start')

    logger.info('Get transcript structure')
    if args.path_annot.find('.gff3') != -1:
        transcripts = ft.gff.get_transcripts_gff3(args.path_annot.split(','), args.convert_ucsc_names, args.path_mapping, args.data_source, args.source_filter)
    elif args.path_annot.find('.bed') != -1:
        transcripts = ft.bed.get_transcripts_bed6(args.path_annot.split(','), args.convert_ucsc_names, args.path_mapping, args.bed_name_as_id, args.bed_feature_id, args.bed_interval_name)
    else:
        logger.error('No compatible input')
        return 1

    if args.path_fasta:
        logger.info('Get transcript sequence')
        transcripts, count_noseq = add_seq(transcripts, args.path_fasta, args.cdna, args.exclude_no_seq)
        if len(count_noseq) > 0:
            noseq_msg = ' '.join([f'{k}:{v}' for k, v in count_noseq.items()])
            logger.info('Missing transcript sequence: '+noseq_msg)

    if args.path_table_transcript and args.path_table_object_xref and args.path_table_xref and args.path_table_ontology and args.path_table_term:
        logger.info('Add gene ontology (GO)')
        db = ft.ensembl_db.EnsemblDBGO(args.path_table_transcript, args.path_table_object_xref, args.path_table_xref, args.path_table_ontology, args.path_table_term)
        missing_terms = db.add_go(transcripts)
        logger.info(f'Missing GO terms: {len(missing_terms)}')

    logger.info(f'Transcript(s) found: {len(transcripts)}')

    # Export
    for bt in args.biotype.split(','):
        # Path to output
        if args.output_format == 'fon':
            if args.path_output.find('$') != -1:
                ot = string.Template(args.path_output)
                path_output = [ot.substitute(biotype=bt, version='1')]
            else:
                path_output = args.path_output.split(',')
        else:
            path_output = [args.path_output]
        # Apply filter
        if bt == 'all':
            ts = transcripts
        else:
            ts = [t for t in transcripts if t['transcript_biotype'] == bt]

        # Save FON1
        if args.output_format == 'fon' or args.output_format == 'fon1':
            logger.info('FON1 export to '+path_output[0])
            json.dump({'fon_version': 1, 'features': ts}, open(path_output[0], 'wt'))

if __name__ == '__main__':
    sys.exit(main())
