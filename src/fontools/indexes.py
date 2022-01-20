# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import math
import os

def get_genome_length(path_tab):
    genome_length = 0
    with open(path_tab, 'rt') as f:
        for l in f:
            genome_length += int(l.split('\t')[-1])
    return genome_length

def star_string_length(genome_length):
    return int(math.floor(min(14, math.log2(genome_length) / 2 - 1)))

class Index(object):
    pass

class Bowtie2(Index):
    def get_exe(self):
        return 'bowtie2-build'

    def get_path(self, fontools_path_main, species_abv, division, release, genome_version_std, suffix):
        return os.path.join(fontools_path_main, 'bowtie2_indices', f'{species_abv}_genome_all_{division}_{genome_version_std}{suffix}')

    def get_create_cmd(self, path_genome, path_gff, path_genome_chrom_length, num_processor):
        return ['bowtie2-build',
                '--threads', str(num_processor),
                path_genome,
                'seq']

class Star(Index):
    def get_exe(self):
        return 'STAR'

    def get_path(self, fontools_path_main, species_abv, division, release, genome_version_std, suffix):
        return os.path.join(fontools_path_main, 'star_indices', f"{species_abv}_genome_all_cdna_all_{division}{release}_{genome_version_std}{suffix}_sjdboverhang75")
    
    def get_create_cmd(self, path_genome, path_gff, path_genome_chrom_length, num_processor):
        return ['STAR',
                '--runMode', 'genomeGenerate',
                '--genomeDir', '.',
                '--genomeFastaFiles', path_genome,
                '--sjdbGTFfile', path_gff,
                '--sjdbGTFtagExonParentTranscript', 'Parent',
                '--sjdbOverhang', '75',
                '--genomeSAindexNbases', str(star_string_length(get_genome_length(path_genome_chrom_length))),
                '--runThreadN', str(num_processor)]

idx2classes = {'bowtie2': Bowtie2,
               'star': Star}
