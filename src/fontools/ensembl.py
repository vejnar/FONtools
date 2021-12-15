# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2021 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import itertools
import json
import os
import re
import string
import subprocess

from . import remote

class DataSource(object):
    genome_levels = []
    taxons = []

class EnsemblSource(DataSource):
    genome_levels = ['primary_assembly', 'toplevel']
    path_genome_root = string.Template('fasta/${species}/dna/')
    path_genome_file = string.Template('${species_title}.(?P<genome_version>\\S+).dna.${genome_level}.fa.gz')
    path_gff = string.Template('gff3/${species}/${species_title}.${genome_version}.${release}.gff3.gz')
    path_cdna = string.Template('fasta/${species}/cdna/${species_title}.${genome_version}.cdna.all.fa.gz')
    path_ncrna = string.Template('fasta/${species}/ncrna/${species_title}.${genome_version}.ncrna.fa.gz')
    path_table_root = 'mysql/'
    rest_url = 'https://rest.ensembl.org/'
    genome_naming_exceptions = [{'re': r'^BDGP6\.\d+$', 'name':'BDGP6'}]

    # Species tables
    table_transcript_names = ['transcript_id',
                              'gene_id',
                              'analysis_id',
                              'seq_region_id',
                              'seq_region_start',
                              'seq_region_end',
                              'seq_region_strand',
                              'display_xref_id',
                              'source',
                              'biotype',
                              'description',
                              'is_current',
                              'canonical_translation_id',
                              'stable_id']
    table_object_xref_names = ['object_xref_id',
                               'ensembl_id',
                               'ensembl_object_type',
                               'xref_id',
                               'linkage_annotation']
    table_xref_names = ['xref_id',
                        'external_db_id',
                        'dbprimary_acc',
                        'display_label',
                        'version',
                        'description',
                        'info_type',
                        'info_text']
    # Ontology tables
    table_ontology_names = ['ontology_id',
                            'name',
                            'namespace']
    table_term_names = ['definition',
                        'term_id',
                        'ontology_id',
                        'subsets',
                        'accession']

    def query_rest(self, query):
        if query == 'species':
            rest_query = self.rest_url + 'info/species?content-type=application/json'
        p = subprocess.run(['wget', '--quiet', '--output-document=-', rest_query], stdout=subprocess.PIPE, check=True)
        return json.loads(p.stdout.decode('utf-8'))

    def get_genome_info(self, species, release):
        prod = [[release]]
        if len(self.taxons) > 0:
            prod.append(self.taxons)
        prod.append([species])
        prod.append([species.capitalize()])
        prod.append(self.genome_levels)
        for p in itertools.product(*prod):
            ip = 0
            if len(self.taxons) > 0:
                taxon = p[1]
                url_path = self.url_path.substitute(release=p[0], taxon=taxon)
                ip += 2
            else:
                taxon = None
                url_path = self.url_path.substitute(release=p[0])
                ip += 1
            path_genome_root = self.path_genome_root.substitute(species=p[ip])
            genome_level = p[ip + 2]
            path_genome_file = self.path_genome_file.substitute(species_title=p[ip + 1], genome_level=genome_level)
            ls, e = remote.rlist(self.url_protocol + url_path + path_genome_root)
            if e:
                for f in ls:
                    rm = re.search(path_genome_file, f)
                    if rm:
                        return self.url_protocol, url_path, taxon, path_genome_root, f, rm.groupdict()['genome_version'], genome_level

    def search_mapping_file(self, genome, path_mapping, suffix='_ensembl2UCSC.txt'):
        p = os.path.join(path_mapping, genome + suffix)
        if os.path.exists(p):
            return p
        for ex in self.genome_naming_exceptions:
            if re.match(ex['re'], genome) is not None:
                p = os.path.join(path_mapping, ex['name'] + suffix)
                if os.path.exists(p):
                    return p
        return None

    def get_species_database_path(self, url_path, species, release):
        path_db = url_path + self.path_table_root
        ls, e = remote.rlist(self.url_protocol + path_db)
        if e:
            for f in ls:
                db = f.strip('/')
                rm = re.search(f'{species}_core_{release}_', db)
                if rm:
                    return path_db + db
        return None

class Ensembl(EnsemblSource):
    taxons = []
    url_protocol = 'http://'
    url_path = string.Template('ftp.ensembl.org/pub/release-${release}/')
    url_path_ontology = string.Template('ftp.ensembl.org/pub/release-${release}/mysql/ensembl_ontology_${release}/')

    def get_ontology_database_path(self, release):
        return self.url_path_ontology.substitute(release=release)

class EnsemblGenomes(EnsemblSource):
    taxons = ['bacteria', 'fungi', 'metazoa', 'plants', 'protists']
    url_protocol = 'http://'
    url_path = string.Template('ftp.ensemblgenomes.org/pub/release-${release}/${taxon}/')
    url_path_ontology = string.Template('ftp.ensemblgenomes.org/pub/release-${release}/pan_ensembl/mysql/')

    def get_ontology_database_path(self, release):
        path_db = self.url_path_ontology.substitute(release=release)
        ls, e = remote.rlist(self.url_protocol + path_db)
        if e:
            for f in ls:
                db = f.strip('/')
                rm = re.search(f'ensembl_ontology_', db)
                if rm:
                    return path_db + db
        return None

def chrom_ensembl2ucsc(ensembl_chrom, name_mapping=None):
    if name_mapping is not None:
        return name_mapping[ensembl_chrom]
    elif ensembl_chrom == 'MT' or ensembl_chrom == 'Mito' or ensembl_chrom == 'MtDNA' or ensembl_chrom == 'mitochondrion_genome':
        return 'chrM'
    elif ensembl_chrom.isdigit() or all([l in ['I', 'V', 'X', 'Y'] for l in ensembl_chrom]):
        return 'chr' + ensembl_chrom
    elif ensembl_chrom == 'Chromosome':
        return 'chromosome'
    else:
        return ensembl_chrom
