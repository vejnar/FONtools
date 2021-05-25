# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2021 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import gzip

from . import ensembl
from . import utils

def bed_parse_seqname(seqname, convert_ucsc, name_mapping):
    if convert_ucsc:
        return ensembl.chrom_ensembl2ucsc(seqname, name_mapping)
    else:
        return seqname

def get_transcripts_bed6(path_annot, convert_ucsc=False, path_mapping=None, bed_name_as_id=False, bed_feature_id='transcript_stable_id', bed_interval_name='exons'):
    # Open BED
    fbed = []
    for p in path_annot:
        if p.endswith('.gz'):
            fbed.append(gzip.open(p, 'rt'))
        else:
            fbed.append(open(p, 'rt'))

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
    transcripts = []
    transcript_num_id = 1
    # Parse BED
    for f in fbed:
        for l in f.readlines():
            if l.startswith('#'):
                continue
            chrom, start, end, name, score, strand = l.split('\t')[:6]
            # ID
            if bed_name_as_id:
                transcript_id = name
            else:
                transcript_id = 't' + str(transcript_num_id)
                transcript_num_id += 1
            # Transcript
            t = {bed_feature_id: transcript_id,
                 'name': name,
                 'chrom': bed_parse_seqname(chrom, convert_ucsc, name_mapping),
                 'strand': strand,
                 bed_interval_name: [[int(start), int(end)]],
                 bed_interval_name+'_on_transcript': [],
                 'score': int(score)}
            t[bed_interval_name+'_on_transcript'], _ = utils.get_exons_on_transcript(t[bed_interval_name], strand)
            transcripts.append(t)

    # Sorting transcripts
    transcripts.sort(key=lambda t: t[bed_feature_id])
    return transcripts
