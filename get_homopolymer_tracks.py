#!/usr/bin/env python
# -*- coding: utf8 -*-

# Copyright (C) by Gaik Tamazian
# mail (at) gtamazian (dot) com

"""Report homopolymer tracks in the BED format.

The minimum length of a track to be reported can be specified as a command-line
option (the default value is 7 bp).
"""

import argparse
import sys
from pyfaidx import Fasta

assert sys.version_info >= (3, 6), 'Python 3.6 or higher required'


def parse_args(args=None):
    parser = argparse.ArgumentParser(description='Get homopolymer tracks')
    parser.add_argument('fasta', help='input FASTA file')
    parser.add_argument('-m', '--min', type=int, default=7,
                        help='minimum length of a track to be reported')
    parser.add_argument('--no-gaps', action='store_true',
                        help='exclude gaps from output')

    return parser.parse_args()


def homopolymer_iterator(seq, alphabet, min_size=0):
    """
    Iterate homopolymer tracks in the specified sequence.
    """
    prev_sym = ''
    start_pos = 0
    for pos, sym in enumerate(seq):
        if sym != prev_sym:
            if (pos - start_pos) >= min_size and prev_sym in alphabet:
                yield (start_pos, pos, prev_sym)
            prev_sym = sym
            start_pos = pos


def get_homopolymer_tracks(args):
    """
    Given command-line arguments, get homopolymer tracks.
    """
    alphabet = {'A', 'C', 'G', 'T', 'N'}
    if args.no_gaps:
        alphabet -= {'N'}

    assert args.min > 1, 'minimum homopolymer track size must be > 1'

    for seq in Fasta(args.fasta):
        for i in homopolymer_iterator(str(seq), alphabet, args.min):
            print('{}\t{}\t{}\t{}'.format(seq.name, i[0], i[1], i[2]))


if __name__ == '__main__':
    args = parse_args()
    get_homopolymer_tracks(args)
