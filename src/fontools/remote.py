# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2022 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import subprocess

def find_links(s, left='href="', right='"'):
    links = []
    i = 0
    while True:
        i = s.find(left, i)
        if i == -1:
            break
        j = s.find(right, i + len(left))
        if j >= 0:
            links.append(s[i+len(left):j])
        i += len(left)
    return links

def rlist(url):
    try:
        p = subprocess.run(['wget', '--quiet', '--output-document=-', url], stdout=subprocess.PIPE, text=True, check=True)
        return find_links(p.stdout), True
    except:
        return None, False

def rget(url, cwd=None):
    return subprocess.run(['wget', '-m', url], cwd=cwd, check=True)
