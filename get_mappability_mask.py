#!/usr/bin/env python
# -*- coding: utf8 -*-

# Copyright (C) 2017 by Gaik Tamazian
# mail (at) gtamazian (dot) com

"""Implements the mappability mask as described in

Mallick, Swapan, et al. "The Simons Genome Diversity Project: 300 genomes from
142 diverse populations." Nature 538.7624 (2016): 201-206. Supplementary
Information section 4.

Note that unlike the paper, we use bowtie2 instead of BWA.
"""

import argparse
import collections
import logging
import os
import pysam
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pyfaidx import Fasta
from tqdm import tqdm

assert sys.version_info >= (3, 6), 'Python 3.6 or higher required'

logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
            '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.WARNING)


def parse_args(args=None):
    parser = argparse.ArgumentParser(description='Construct mappability mask')
    parser.add_argument('fasta', help='input FASTA file')
    parser.add_argument('seq', help='sequence name from input file')
    parser.add_argument('k', type=int, help='k-mer size')
    parser.add_argument('index', help='bowtie2 index filename prefix')
    parser.add_argument('-p', '--threads', type=int, default=1,
                        help='number of bowtie2 threads')
    parser.add_argument('-b', '--buffer', type=int, default=10000,
                        help='k-mer buffer size in sequences')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='show detailed log messages')

    return parser.parse_args(args)


def get_bowtie2_caller(index, threads=1):
    """
    Return a function that calls bowtie2 to align sequences from an input FASTA
    file to the specified index and write the alignments to an output SAM file.

    The optional argument specifies the number of bowtie2's threads.
    """
    call_prefix = ['bowtie2', '-f', '-k', '2', '-p', str(threads), '-x', index,
                   '--no-unal', '--omit-sec-seq']
    devnull = open(os.devnull, 'w')

    def align_reads(in_fname, out_fname):
        call_list = call_prefix + ['-U', in_fname, '-S', out_fname]
        logger.info('subprocess call: %s', ' '.join(call_list))
        return subprocess.call(call_list, stderr=devnull)

    return align_reads


def create_kmer_file(seq, k, out_fname):
    """
    Write k-mers from the sequence to the file in the FASTA format.

    K-mers contaning only Ns are skipped.
    """
    logger.info('writing %d-mers of %d bp sequence to %s', k, len(seq),
                out_fname)
    kmer_num = 0
    with open(out_fname, 'w') as output_file:
        for i in range(len(seq) - k + 1):
            if set(seq[i:i+k]) != {'N'}:
                output_file.write('>{}\n{}\n'.format(i, seq[i:i+k]))
                kmer_num += 1
    logger.info('%d %d-mers were written', kmer_num, k)
    return kmer_num


def get_temp_fname():
    """
    A secure workaround to get a name of a temporary file.
    """
    f = tempfile.NamedTemporaryFile(delete=False)
    fname = f.name
    f.close()
    return fname


def get_alignment_numbers(sam_fname, kmer_num=None):
    """
    Return a list of alignment numbers for k-mers from the SAM file.
    """
    logger.info('counting alignments from %s', sam_fname)
    num_alignments = [0] * kmer_num
    with pysam.AlignmentFile(sam_fname, 'r') as sam_file:
        for record in sam_file:
            num_alignments[int(record.qname)] += 1

    return num_alignments


@contextmanager
def kmer_alignment_uniqueness(seq, k, caller, buf_size=1000):
    """
    Iterate k-mers of the sequence and get for each whether it maps uniquely by
    the caller.
    """
    temp_fa = get_temp_fname()
    temp_sam = get_temp_fname()

    def kmer_alignment_uniqueness_iterator():
        buf_pos = buf_size
        buf = []
        for i in range(len(seq) - k + 1):
            if buf_pos == buf_size:
                kmer_num = create_kmer_file(seq[i:i + k + buf_size], k,
                                            temp_fa)
                if kmer_num > 0:
                    caller(temp_fa, temp_sam)
                    buf = get_alignment_numbers(temp_sam, len(seq) - k + 1)
                else:
                    logger.info('zero %d-mers for [%d, %d)', k, i,
                                i + k + buf_size)
                    buf = [0] * (len(seq) - k + 1)
                buf_pos = 0
            yield buf[buf_pos] == 1
            buf_pos += 1
    yield kmer_alignment_uniqueness_iterator

    if os.path.isfile(temp_fa):
        os.unlink(temp_fa)
    if os.path.isfile(temp_sam):
        os.unlink(temp_sam)


def mapscore_iterator(uniq_list, k):
    """
    Iterate mappability scores given k-mer mapping uniq values and the
    k-mer size.
    """
    if len(uniq_list) < k:
        raise StopIteration

    uniq_queue = collections.deque(maxlen=k)

    # process positions from 0 to k-1
    for i in range(k):
        uniq_queue.append(uniq_list[i])
        yield sum(uniq_queue)/len(uniq_queue)

    # process positions from k to n - k, where n is the length of uniq_list
    for i in range(len(uniq_list) - k):
        uniq_queue.append(uniq_list[i])
        yield sum(uniq_queue)/k

    # process last k positions
    for i in range(len(uniq_list) - k, len(uniq_list)):
        yield sum(uniq_queue)/len(uniq_queue)
        uniq_queue.popleft()


def print_mapbed(seq_name, scores):
    """
    Print mappability scores for sequence intervals in the BED format.
    """
    prev_value = next(scores)
    cur_pos = prev_pos = 0
    for i, v in enumerate(scores):
        if v != prev_value:
            print('{}\t{}\t{}\t{}'.format(seq_name, prev_pos, i,
                                          round(prev_value, 2)))
            prev_pos, prev_value = i, v
        cur_pos += 1

    # output the last line
    print('{}\t{}\t{}\t{}'.format(seq_name, prev_pos, cur_pos, prev_value))


def get_mappability_mask(args):
    """
    Given command-line arguments, get the mappability mask.
    """
    if args.verbose:
        logger.setLevel(logging.INFO)

    seq = str(Fasta(args.fasta)[args.seq])
    logger.info('%s %d bp length loaded', args.seq, len(seq))

    logger.info('%d-mer mapping started', args.k)
    kmer_uniqueness = []
    bw2 = get_bowtie2_caller(args.index, args.threads)
    with kmer_alignment_uniqueness(seq, args.k, bw2, args.buffer) as kau:
        for i in tqdm(kau(), total=len(seq) - args.k + 1,
                      desc='Mapping k-mers'):
            kmer_uniqueness.append(i)
    logger.info('%d-mer mapping completed', args.k)

    print_mapbed(args.seq, mapscore_iterator(kmer_uniqueness, args.k))


if __name__ == '__main__':
    args = parse_args()
    get_mappability_mask(args)
