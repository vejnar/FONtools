# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2021 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import collections
import csv
import gzip
import sys

from . import ensembl

def open_table(path_tsv, index_col=None, index_integer=True, names=['id'], name='Table', list_out=True, filters=[]):
    if index_col is None:
        index_col = names[0]
    ndt = collections.namedtuple(name, names)
    with gzip.open(path_tsv, mode='rt') as fcsv:
        table = {}
        csv.field_size_limit(sys.maxsize)
        for row in csv.reader(fcsv, delimiter='\t'):
            rec = ndt(*row[:len(names)])
            filter_ok = True
            if len(filters) > 0:
                for flt in filters:
                    if getattr(rec, flt[0]) != flt[1]:
                        filter_ok = False
                        break
            if filter_ok:
                idx = getattr(rec, index_col)
                if index_integer:
                    idx = int(idx)
                if idx in table:
                    if list_out:
                        table[idx].append(rec)
                    else:
                        raise ValueError('Duplicate entry')
                else:
                    if list_out:
                        table[idx] = [rec]
                    else:
                        table[idx] = rec
    return table

class EnsemblDBGO(object):
    def __init__(self, path_table_transcript, path_table_object_xref, path_table_xref, path_table_ontology, path_table_term):
        # Source
        source = ensembl.Ensembl()
        # Open species tables
        self.table_transcript = open_table(path_table_transcript,
                                           index_col = 'stable_id',
                                           index_integer = False,
                                           names = source.table_transcript_names,
                                           list_out = False)
        self.table_object_xref = open_table(path_table_object_xref,
                                            index_col = 'ensembl_id',
                                            names = source.table_object_xref_names,
                                            filters=[['ensembl_object_type', 'Transcript']])
        self.table_xref = open_table(path_table_xref,
                                     names = source.table_xref_names,
                                     list_out = False,
                                     filters = [['external_db_id', '1000']])
        # Open ontology tables
        self.table_ontology = open_table(path_table_ontology,
                                         names = source.table_ontology_names,
                                         filters = [['name', 'GO']],
                                         list_out = False)
        self.table_term = open_table(path_table_term,
                                     index_col = 'accession',
                                     index_integer = False,
                                     names = source.table_term_names,
                                     list_out = False)

    def add_go(self, transcripts):
        missing_terms = set()
        for transcript in transcripts:
            # New terms
            terms = {}
            # Get transcript internal ID
            transcript_id = int(self.table_transcript[transcript['transcript_stable_id']].transcript_id)
            # Has transcript external reference(s)?
            if transcript_id in self.table_object_xref:
                # Union to object_xref table
                for o in self.table_object_xref[transcript_id]:
                    xref_id = int(o.xref_id)
                    # Union to xref table
                    if xref_id in self.table_xref:
                        x = self.table_xref[xref_id]
                        # Union to term table
                        if x.dbprimary_acc in self.table_term:
                            ontology_id = int(self.table_term[x.dbprimary_acc].ontology_id)
                            # Get sources
                            if o.linkage_annotation == '\\N':
                                sources = []
                            else:
                                sources = [{'name':x.info_text, 'evidence':o.linkage_annotation}]
                            # New term record
                            rec = {'term': x.description,
                                'domain': self.table_ontology[ontology_id].namespace,
                                'sources': sources}
                            # Add to terms for transcript
                            if x.dbprimary_acc in terms:
                                if len(sources) > 0:
                                    terms[x.dbprimary_acc]['sources'].extend(sources)
                            else:
                                terms[x.dbprimary_acc] = rec
                        else:
                            missing_terms.add(x.dbprimary_acc)
            transcript['go'] = terms
        return missing_terms
