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
import math
import os
import re
import string
import shutil
import subprocess
import sys

import fontools as ft
import pyfnutils as pfu
import pyfnutils.log

def run_cmd(cmd, logger, cwd=None):
    logger.info('Start ' + str(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)

def check_exe(names):
    for name in names:
        if shutil.which(name) == None:
            raise FileNotFoundError(name+' missing')

def main(argv=None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description='Import Ensembl data.')
    parser.add_argument('-o', '--path_main', '--fontools_path_main', dest='fontools_path_main', action='store', help='Path to FONtools directory.')
    parser.add_argument('-d', '--path_download', '--fontools_path_download', dest='fontools_path_download', action='store', help='Path to download directory.')
    parser.add_argument('-m', '--path_mapping', '--fontools_path_mapping', dest='fontools_path_mapping', action='store', help='Path to chromosome name mappings.')
    parser.add_argument('-l', '--path_log', dest='path_log', action='store', help='Log path.')
    parser.add_argument('-n', '--division', dest='division', action='store', default='ensembl', help='Ensembl division (ensembl or ensembl_genomes).')
    parser.add_argument('-s', '--species', dest='species', action='store', required=True, help='Ensembl species name: all, list or name (e.g. danio_rerio) (comma separated).')
    parser.add_argument('-r', '--release', dest='release', action='store', required=True, help='Ensembl release.')
    parser.add_argument('-a', '--path_rest_cache', dest='path_rest_cache', action='store', help='Path to REST cache.')
    parser.add_argument('-x', '--skip_rest_cache', dest='skip_rest_cache', action='store_true', help='Don\'t use REST cache.')
    parser.add_argument('-t', '--steps', dest='steps', action='store', default='all', help='Get all or gg (=genome,gene) or genome,gene,bowtie2,star (comma separated).')
    parser.add_argument('-u', '--ucsc_naming', dest='ucsc_naming', action='store_true', default=False, help='Convert to UCSC chromosome/scaffold naming.')
    parser.add_argument('-g', '--import_go', dest='import_go', action='store_true', default=False, help='Import Gene Ontology.')
    parser.add_argument('-z', '--compress', dest='compress', action='store_true', default=False, help='Compress output when possible.')
    parser.add_argument('-p', '--processor', dest='num_processor', action='store', type=int, default=1, help='Number of processor')
    parser.add_argument('--path_config', dest='path_config', action='store', help='Path to config')
    args = parser.parse_args(argv[1:])

    # Get config (JSON single file or all files in path_config)
    config = {}
    paths = []
    if args.path_config is None:
        if 'HTS_CONFIG_PATH' in os.environ:
            paths.append(os.environ['HTS_CONFIG_PATH'])
        elif 'XDG_CONFIG_HOME' in os.environ:
            paths.append(os.path.join(os.environ['XDG_CONFIG_HOME'], 'hts'))
    else:
        paths.append(args.path_config)
    for path in paths:
        if os.path.isdir(path):
            for f in sorted(os.listdir(path)):
                if f.endswith('.json'):
                    config = {**config, **json.load(open(os.path.join(path, f)))}
        elif os.path.isfile(path):
            config = {**config, **json.load(open(path))}

    # Input local config from args
    vargs = vars(args)
    for a, v in vargs.items():
        if v is not None and (a not in config or v != parser.get_default(a)):
            config[a] = v
        if a in ['steps', 'species']:
            if v is None:
                config[a] = []
            else:
                config[a] = [r.strip() for r in v.split(',')]

    # Steps
    if config['steps'][0] == 'all':
        steps = ['genome', 'gene', 'bowtie2', 'star']
    elif config['steps'][0] == 'gg':
        steps = ['genome', 'gene']
    else:
        steps = config['steps']

    # Check directories
    for p in ['fontools_path_main', 'fontools_path_download']:
        if p not in config:
            print(f'ERROR: --{p} required')
            return 1
        if not os.path.exists(config[p]):
            print(f'ERROR: {config[p]} not found')
            return 1
    # Check executables
    exes = ['wget']
    if len(config['species']) > 0 and config['species'][0] != 'list':
        if 'genome' in steps:
            exes.extend(['fasta_format', 'fasta_seq_length'])
        if 'gene' in steps:
            exes.extend(['fon_import', 'fon_transform'])
        if config['ucsc_naming'] or 'fontools_path_mapping' in config:
            exes.append('ensembl2ucsc')
        for idxname, idxc in ft.indexes.idx2classes.items():
            if idxname in steps:
                idx = idxc()
                exes.append(idx.get_exe())
        if config['compress']:
            exes.append('zstd')
    check_exe(exes)

    # Source
    if config['division'] == 'ensembl':
        source = ft.ensembl.Ensembl()
    elif config['division'] == 'ensembl_genomes':
        source = ft.ensembl.EnsemblGenomes()
    else:
        raise ValueError('Unknown division: '+config['division'])

    # Species
    if config['species'][0] == 'all' or config['species'][0] == 'list':
        if config['skip_rest_cache']:
            jp = source.query_rest('species')
        else:
            if 'path_rest_cache' in config:
                path_cache = config['path_rest_cache']
            else:
                path_cache = os.path.join(config['fontools_path_main'], 'cache')
            path_rest_species = os.path.join(path_cache, config['division']+config['release']+'_species.json')
            if os.path.exists(path_rest_species):
                jp = json.load(open(path_rest_species))
            else:
                jp = source.query_rest('species')
                if os.path.exists(path_cache):
                    json.dump(jp, open(path_rest_species, 'wt'))
        if config['species'][0] == 'all':
            lspecies = [s['name'] for s in jp['species']]
        elif config['species'][0] == 'list':
            for s in sorted(jp['species']):
                print(s['name'])
            return 0
    else:
        lspecies = config['species']

    # Log
    if 'path_log' in config:
        path_log = config['path_log']
    else:
        path_log = os.path.join(config['fontools_path_main'], 'log')
        if os.path.exists(path_log):
            path_log = os.path.join(path_log, config['division']+config['release']+'.log')
        else:
            path_log = None
    logger = pfu.log.define_root_logger('import_'+config['division'], filename=path_log)

    for species in lspecies:
        logger.info(f"Starting ({species},{config['release']})")
        species_abv = ft.naming.get_species_abv(species)

        # Get remote path and info
        url_protocol, url_path, taxon, path_genome_root, path_genome_file, genome_version, genome_level = source.get_genome_info(species, config['release'])

        # Genome path and assembly
        path_genome = path_genome_root + path_genome_file
        genome_version_std = ft.naming.get_genome_version_std(genome_version)
        logger.info(f'Found assembly {genome_version} ({genome_level})')

        # Local paths
        suffix = ''
        if config['ucsc_naming']:
            suffix += '_ucsc_names'
        path_genome_local = os.path.join(config['fontools_path_main'], 'seqs', f"{species_abv}_genome_all_{config['division']}_{genome_version_std}{suffix}.fa")
        path_gff_local = os.path.join(config['fontools_path_main'], 'annots', f"{species_abv}_cdna_all_{config['division']}{config['release']}{suffix}.gff3")
        path_genome_chrom_length_local = os.path.join(config['fontools_path_main'], 'annots', f"{species_abv}_genome_all_{config['division']}_{genome_version_std}{suffix}_chrom_length.tab")
        path_fon_local = os.path.join(config['fontools_path_main'], 'annots', f"{species_abv}_cdna_${{biotype}}_{config['division']}{config['release']}{suffix}.fon${{version}}.json")

        # Create main folder
        if not os.path.exists(config['fontools_path_main']):
            os.mkdir(config['fontools_path_main'])

        # Download: Genome
        if 'genome' in steps:
            if os.path.exists(path_genome_local):
                logger.info('Found ' + path_genome_local)
            else:
                p = os.path.join(config['fontools_path_download'], url_path, path_genome)
                if os.path.exists(p):
                    logger.info('Found ' + p)
                else:
                    logger.info('Downloading ' + path_genome)
                    ft.remote.rget(url_protocol + url_path + path_genome, cwd=config['fontools_path_download'])

        # Download: Genes
        path_gff = source.path_gff.substitute(species=species, species_title=species.capitalize(), genome_version=genome_version, release=config['release'])
        path_cdna = source.path_cdna.substitute(species=species, species_title=species.capitalize(), genome_version=genome_version)
        path_ncrna = source.path_ncrna.substitute(species=species, species_title=species.capitalize(), genome_version=genome_version)
        if 'gene' in steps:
            for p in [path_gff, path_cdna, path_ncrna]:
                path_file_full = os.path.join(config['fontools_path_download'], url_path, p)
                if os.path.exists(path_file_full):
                    logger.info('Found ' + path_file_full)
                else:
                    logger.info('Downloading ' + p)
                    ft.remote.rget(url_protocol + url_path + p, cwd=config['fontools_path_download'])

            # Download GO
            if config['import_go']:
                path_species_db = source.get_species_database_path(url_path, species, config['release'])
                path_go_db = source.get_ontology_database_path(config['release'])
                for path_db, table in [(path_species_db, 'transcript'),
                                       (path_species_db, 'object_xref'),
                                       (path_species_db, 'xref'),
                                       (path_go_db, 'ontology'),
                                       (path_go_db, 'term')]:
                    path_table_full = os.path.join(config['fontools_path_download'], path_db, f'{table}.txt.gz')
                    if os.path.exists(path_table_full):
                        logger.info('Found ' + path_table_full)
                    else:
                        logger.info('Downloading ' + p)
                        ft.remote.rget(url_protocol + f'{path_db}/{table}.txt.gz', cwd=config['fontools_path_download'])

        # Sort/Convert/Copy chromosome names
        for s, pr, p, pl in [('genome', 'seqs', path_genome, path_genome_local), ('gene', 'annots', path_gff, path_gff_local)]:
            if s in steps:
                if os.path.exists(pl):
                    logger.info('Found ' + pl)
                else:
                    # Output folder
                    if not os.path.exists(os.path.join(config['fontools_path_main'], pr)):
                        os.mkdir(os.path.join(config['fontools_path_main'], pr))
                    # Export
                    p = os.path.join(config['fontools_path_download'], url_path, p)
                    if config['ucsc_naming']:
                        logger.info('Exporting to ' + pl)
                        path_mapping = source.search_mapping_file(genome_version, config['fontools_path_mapping'])
                        if path_mapping is not None:
                            run_cmd(['ensembl2ucsc',
                                     '--input', p,
                                     '--output', pl,
                                     '--path_mapping', path_mapping], logger)
                        else:
                            logger.error('UCSC chromosome name(s) mapping file not found for '+genome_version)
                    else:
                        if s == 'genome':
                            logger.info('Sorting to ' + pl)
                            run_cmd(['fasta_format',
                                     '--sort',
                                     '--input', p,
                                     '--output', pl], logger)
                        elif s == 'gene':
                            run_cmd(['cp', p, pl + '.gz'], logger)
                            run_cmd(['gzip', '-d', pl + '.gz'], logger)

        # Chromosome index and lengths
        if 'genome' in steps:
            if os.path.exists(path_genome_chrom_length_local):
                logger.info('Found ' + path_genome_chrom_length_local)
            else:
                logger.info('Creating chromosome length file ' + path_genome_chrom_length_local)
                run_cmd(['fasta_seq_length',
                         '--input', path_genome_local,
                         '--output', path_genome_chrom_length_local], logger)
            if not config['ucsc_naming'] and 'fontools_path_mapping' in config:
                path_mapping = source.search_mapping_file(genome_version, config['fontools_path_mapping'])
                path_genome_chrom_length_local_ucsc = path_genome_chrom_length_local.replace('_chrom_length.tab', '_ucsc_names_chrom_length.tab')
                if os.path.exists(path_genome_chrom_length_local_ucsc):
                    logger.info('Found ' + path_genome_chrom_length_local_ucsc)
                else:
                    if path_mapping is not None:
                        logger.info('Creating chromosome length file for UCSC ' + path_genome_chrom_length_local_ucsc)
                        run_cmd(['ensembl2ucsc',
                                 '--input', path_genome_chrom_length_local,
                                 '--output', path_genome_chrom_length_local_ucsc,
                                 '--path_mapping', path_mapping], logger)
                    else:
                        logger.warning('UCSC chromosome name(s) mapping file not found for '+genome_version)

        # Import annotation
        if 'gene' in steps:
            path_fon = string.Template(path_fon_local).substitute(biotype='all', version='1')
            if os.path.exists(path_fon):
                logger.info('Found ' + path_fon)
            elif os.path.exists(path_fon+'.zst'):
                logger.info('Found ' + path_fon + '.zst')
            else:
                logger.info('Importing annotation')
                cmd = ['fon_import',
                       '--annotation', os.path.join(config['fontools_path_download'], url_path, path_gff),
                       '--data_source', 'ensembl',
                       '--fasta', os.path.join(config['fontools_path_download'], url_path, path_cdna),
                       '--fasta', os.path.join(config['fontools_path_download'], url_path, path_ncrna),
                       '--cdna',
                       '--exclude_no_seq',
                       '--biotype', 'all,protein_coding',
                       '--output', path_fon_local,
                       '--output_format', 'fon']
                if config['ucsc_naming']:
                    cmd.append('--ucsc_names')
                # Add GO annotation
                if config['import_go']:
                    cmd.extend(['--table_transcript', os.path.join(config['fontools_path_download'], path_species_db, 'transcript.txt.gz'),
                                '--table_object_xref', os.path.join(config['fontools_path_download'], path_species_db, 'object_xref.txt.gz'),
                                '--table_xref', os.path.join(config['fontools_path_download'], path_species_db, 'xref.txt.gz'),
                                '--table_ontology', os.path.join(config['fontools_path_download'], path_go_db, 'ontology.txt.gz'),
                                '--table_term', os.path.join(config['fontools_path_download'], path_go_db, 'term.txt.gz')])
                run_cmd(cmd, logger)

        # Select/merge transcript(s)
        if 'gene' in steps:
            for biotype, method, method_name in [('protein_coding', 'union', 'union2gene'),
                                                 ('protein_coding', 'longest', 'longest_transcript'),
                                                 ('all', 'union', 'union2gene'),
                                                 ('all', 'longest', 'longest_transcript')]:
                path_fon = string.Template(path_fon_local).substitute(biotype=method_name + '_' + biotype, version='1')
                if os.path.exists(path_fon):
                    logger.info('Found ' + path_fon)
                elif os.path.exists(path_fon+'.zst'):
                    logger.info('Found ' + path_fon + '.zst')
                else:
                    logger.info(f'Transform FON ({method},{biotype})')
                    run_cmd(['fon_transform',
                             '--fon', string.Template(path_fon_local).substitute(biotype=biotype, version='1'),
                             '--method', method,
                             '--output', string.Template(path_fon_local).safe_substitute(biotype=method_name + '_' + biotype)], logger)

        # Compress
        if 'gene' in steps:
            cmd = ['zstd', '--rm', '-T'+str(config['num_processor']), '-19']
            path_annot = os.path.dirname(path_fon_local)
            for f in os.listdir(path_annot):
                if f.endswith('.fon1.json'):
                    run_cmd(cmd + [os.path.join(path_annot, f)], logger)

        # Indexes
        for idxname, idxc in ft.indexes.idx2classes.items():
            if idxname in steps:
                idx = idxc()
                path_idx_local = idx.get_path(config['fontools_path_main'], species_abv, config['division'], config['release'], genome_version_std, suffix)
                if os.path.exists(path_idx_local):
                    logger.info('Found index ' + path_idx_local)
                else:
                    logger.info(f"Creating {idxname.title()} index in {path_idx_local}")
                    # Output folder
                    os.makedirs(path_idx_local)
                    # Create index
                    run_cmd(idx.get_create_cmd(path_genome_local, path_gff_local, path_genome_chrom_length_local, config['num_processor']), logger, cwd=path_idx_local)

if __name__ == '__main__':
    sys.exit(main())
