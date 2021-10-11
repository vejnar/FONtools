# -*- coding: utf-8 -*-

#
# Copyright (C) 2015-2021 Charles E. Vejnar
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
#

import itertools

def strand2direction(strand):
    if strand == '-':
        return -1
    else:
        return 1

def flatten(l1d):
    return list(itertools.chain.from_iterable(l1d))

def get_intervals_union(intervals):
    # Sort intervals
    intervals_idx = sorted(range(len(intervals)), key=lambda k: intervals[k][0])
    intervals = [intervals[i][:] for i in intervals_idx]
    # Union
    y = [intervals[0]]
    for x in intervals[1:]:
        if y[-1][1] < x[0]:
            y.append(x)
        elif y[-1][1] == x[0]:
            y[-1][1] = x[1]
        if x[1] > y[-1][1]:
            y[-1][1] = x[1]
    return y

def get_intervals_seq_union(intervals, seqs):
    # Sort intervals
    intervals_idx = sorted(range(len(intervals)), key=lambda k: intervals[k][0])
    intervals = [intervals[i][:] for i in intervals_idx]
    # Sort sequences
    if len(seqs) == len(intervals):
        seqs = [seqs[i] for i in intervals_idx]
        s = [seqs[0]]
        do_seq = True
    else:
        s = []
        do_seq = False
    # Union
    y = [intervals[0]]
    ix = 1
    for x in intervals[1:]:
        if y[-1][1] < x[0]:
            y.append(x)
            if do_seq:
                s.append(seqs[ix])
        elif y[-1][1] == x[0]:
            if do_seq:
                s[-1] += seqs[ix]
            y[-1][1] = x[1]
        if x[1] > y[-1][1]:
            if do_seq:
                s[-1] += seqs[ix][(x[1] - y[-1][1]) * -1:]
            y[-1][1] = x[1]
        ix += 1
    return y, ''.join(s)

def get_exons_on_transcript(exons, strand='+', tlength=0):
    exons_on_transcript = []
    for exon in exons[::strand2direction(strand)]:
        exon_len = exon[1] - exon[0]
        exons_on_transcript.append([tlength, tlength + exon_len])
        tlength += exon_len
    return exons_on_transcript, tlength

def get_intervals_left_disjoint(outer_intervals, inner_intervals):
    ds = []
    for x in outer_intervals:
        for y in inner_intervals:
            if x[0] == y[0]:
                return ds
            elif x[0] < y[0]:
                if y[0] < x[1]:
                    ds.append([x[0], y[0]])
                    return ds
                else:
                    ds.append([x[0], x[1]])
                    break
    return ds

def get_intervals_right_disjoint(outer_intervals, inner_intervals):
    ds = []
    for x in outer_intervals[::-1]:
        for y in inner_intervals[::-1]:
            if x[1] == y[1]:
                return ds[::-1]
            elif y[1] < x[1]:
                if x[0] < y[1]:
                    ds.append([y[1], x[1]])
                    return ds[::-1]
                else:
                    ds.append([x[0], x[1]])
                    break
    return ds[::-1]
