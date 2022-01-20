# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import re

def get_species_abv(species, abl=3):
    return ''.join([s[:abl] for s in species.split('_')])

def get_genome_version_std(genome_version):
    return re.sub(r'[^a-zA-Z0-9]', '', genome_version.lower())
