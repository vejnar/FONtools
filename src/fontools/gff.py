# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2021 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import gzip
import urllib.parse

from . import ensembl
from . import utils

def gff_parse_attributes(attributes_raw):
    d = {}
    for a in attributes_raw.split(';'):
        if len(a) > 0:
            i = a.find('=')
            if i >= 0:
                b, c = a[:i], a[i+1:]
                d[b] = urllib.parse.unquote(c)
    return d

def gff_parse_seqname(seqname, convert_ucsc, name_mapping):
    if convert_ucsc:
        return ensembl.chrom_ensembl2ucsc(seqname, name_mapping)
    else:
        return seqname

def gff_parse_phase(phase):
    if phase == '1':
        return 2
    elif phase == '2':
        return 1
    else:
        return 0

def get_transcripts_gff3(path_annot, convert_ucsc=False, path_mapping=None, data_source=None, filter_source=[]):
    # Open GFF
    fgff = []
    for p in path_annot:
        if p.endswith('.gz'):
            fgff.append(gzip.open(p, 'rt'))
        else:
            fgff.append(open(p, 'rt'))

    # Open chromosome/scaffold name mapping
    name_mapping = None
    if path_mapping is not None:
        name_mapping = {}
        with open(path_mapping, 'rt') as f:
            for l in f:
                fields = l.strip('\n').split('\t')
                if len(fields) == 2:
                    if fields[1] == '':
                        name_mapping[fields[0]] = None
                    else:
                        name_mapping[fields[0]] = fields[1]

    # Init.
    genes = {}
    transcripts = {}
    # GFF3 specifications authorize exon and CDS to be directly the child of a gene (ie. without an mRNA)
    # wo_transcript_genes stores made-up transcript for these genes
    wo_transcript_genes = {}
    transcript_num_id = 1
    # Parse GFF
    for f in fgff:
        for l in f.readlines():
            if l.startswith('#'):
                continue
            seqname, source, feature, start, end, score, strand, phase, attributes_raw = l.strip().split('\t')
            if filter_source is None or len(filter_source) == 0 or source in filter_source:
                attributes = gff_parse_attributes(attributes_raw)
                # Gene (i.e. top of target hierarchy)
                if 'ID' in attributes:
                    # Parse ID type
                    if data_source == 'ensembl':
                        id_type = attributes['ID'].split(':')[0]
                    else:
                        id_type = feature
                    if id_type == 'gene':
                        if data_source == 'ensembl':
                            gene_id = attributes['gene_id']
                            gene_biotype_column = 'biotype'
                        elif data_source == 'xenbase':
                            gene_id = attributes['ID']
                            gene_biotype_column = 'gene_biotype'
                        else:
                            gene_id = attributes['ID']
                            gene_biotype_column = 'biotype'
                        # Add gene
                        genes[gene_id] = {'gene_name': attributes.get('Name'),
                                          'gene_version': attributes.get('version'),
                                          'gene_biotype': attributes.get(gene_biotype_column),
                                          'attributes': attributes}
                # Everything else under gene
                if 'Parent' in attributes:
                    if data_source == 'ensembl':
                        parent_id = attributes['Parent'].split(':')[1]
                    else:
                        parent_id = attributes['Parent']
                    # Transcript
                    if parent_id in genes and parent_id not in wo_transcript_genes:
                        if data_source == 'ensembl' and 'transcript_id' in attributes:
                            transcript_id = attributes['transcript_id']
                        elif feature == 'exon' or feature == 'CDS':
                            transcript_id = 't' + str(transcript_num_id)
                            transcript_num_id += 1
                            wo_transcript_genes[parent_id] = transcript_id
                        elif 'ID' in attributes:
                            transcript_id = attributes['ID']
                        else:
                            transcript_id = 't' + str(transcript_num_id)
                            transcript_num_id += 1
                        if 'biotype' in attributes:
                            transcript_biotype = attributes['biotype']
                        else:
                            transcript_biotype = genes[parent_id]['gene_biotype']
                        transcripts[transcript_id] = {'transcript_stable_id': transcript_id,
                                                      'gene_stable_id': parent_id,
                                                      'gene_name': genes[parent_id]['gene_name'],
                                                      'protein_stable_id': None,
                                                      'chrom': gff_parse_seqname(seqname, convert_ucsc, name_mapping),
                                                      'strand': strand,
                                                      'transcript_version': attributes.get('version'),
                                                      'gene_version': genes[parent_id]['gene_version'],
                                                      'transcript_biotype': transcript_biotype,
                                                      'gene_biotype': genes[parent_id]['gene_biotype'],
                                                      'exons': [],
                                                      'exons_on_transcript': [],
                                                      'cds_exons': [],
                                                      'cds_exons_on_transcript': [],
                                                      'cds_exons_frame': [],
                                                      'cds_exons_frame_on_transcript': [],
                                                      'utr5_exons': [],
                                                      'utr5_exons_on_transcript': [],
                                                      'utr3_exons': [],
                                                      'utr3_exons_on_transcript': []}
                    if parent_id in transcripts:
                        cid = parent_id
                    elif parent_id in wo_transcript_genes and (feature == 'exon' or feature == 'CDS'):
                        cid = wo_transcript_genes[parent_id]
                    else:
                        cid = None
                    if cid:
                        if feature == 'exon':
                            transcripts[cid]['exons'].append([int(start) - 1, int(end)])
                        elif feature == 'CDS':
                            transcripts[cid]['cds_exons'].append([int(start) - 1, int(end)])
                            transcripts[cid]['cds_exons_frame'].append(gff_parse_phase(phase))
                            # Add protein ID
                            if transcripts[cid]['protein_stable_id'] is None and 'protein_id' in attributes:
                                transcripts[cid]['protein_stable_id'] = attributes['protein_id']
                        elif feature == 'five_prime_UTR':
                            transcripts[cid]['utr5_exons'].append([int(start) - 1, int(end)])
                        elif feature == 'three_prime_UTR':
                            transcripts[cid]['utr3_exons'].append([int(start) - 1, int(end)])
        f.close()

    # Get transcript list
    transcripts = list(transcripts.values())
    for t in transcripts:
        # Sort exon(s)
        t['exons'].sort(key=lambda e: e[0])
        t['utr5_exons'].sort(key=lambda e: e[0])
        t['utr3_exons'].sort(key=lambda e: e[0])
        if len(t['cds_exons']) > 0:
            cds_exons_idx = sorted(range(len(t['cds_exons'])), key=lambda k: t['cds_exons'][k][0])
            t['cds_exons'] = [t['cds_exons'][i] for i in cds_exons_idx]
            t['cds_exons_frame'] = [t['cds_exons_frame'][i] for i in cds_exons_idx]

        # Add missing exons ('exon' type of feature for each 'CDS' are not mandatory in GFF3)
        if len(t['exons']) == 0 and (len(t['utr5_exons']) > 0 or len(t['cds_exons']) > 0 or len(t['utr3_exons']) > 0):
            t['exons'] = utils.get_intervals_union(t['utr5_exons'] + t['cds_exons'] + t['utr3_exons'])

        # Add missing UTR exons ('five_prime_UTR' and 'three_prime_UTR' are not mandatory in GFF3)
        if len(t['utr5_exons']) == 0 and len(t['cds_exons']) > 0:
            if t['strand'] == '+':
                t['utr5_exons'] = utils.get_intervals_left_disjoint(t['exons'], t['cds_exons'])
            elif t['strand'] == '-':
                t['utr5_exons'] = utils.get_intervals_right_disjoint(t['exons'], t['cds_exons'])
        if len(t['utr3_exons']) == 0 and len(t['cds_exons']) > 0:
            if t['strand'] == '+':
                t['utr3_exons'] = utils.get_intervals_right_disjoint(t['exons'], t['cds_exons'])
            elif t['strand'] == '-':
                t['utr3_exons'] = utils.get_intervals_left_disjoint(t['exons'], t['cds_exons'])

        # Add "_on_transcript" fields
        t['exons_on_transcript'], tlength = utils.get_exons_on_transcript(t['exons'], t['strand'])
        t['utr5_exons_on_transcript'], tlength = utils.get_exons_on_transcript(t['utr5_exons'], t['strand'])
        if len(t['cds_exons']) > 0:
            t['cds_exons_on_transcript'], tlength = utils.get_exons_on_transcript(t['cds_exons'], t['strand'], tlength)
            t['cds_exons_frame_on_transcript'] = t['cds_exons_frame'][::utils.strand2direction(t['strand'])]
        else:
            t['cds_exons_on_transcript'] = []
            t['cds_exons_frame_on_transcript'] = []
        t['utr3_exons_on_transcript'], tlength = utils.get_exons_on_transcript(t['utr3_exons'], t['strand'], tlength)

    # Sorting transcripts
    transcripts.sort(key=lambda t: t['transcript_stable_id'])
    return transcripts
