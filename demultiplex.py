#!/usr/bin/env python
"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode

Output file includes the input file name and the barcode:
e.g., file.fastq demultiplexes into file_CGAT.fastq, etc.


Generalize to read in the input file and the barcode file

Log the number of reads in the input file and the number in each barcode
(and the number and percent that do not belong to one of those barcodes)

"""

import gzip
import sys # for output logging
import os
import argparse


def barcodes_prep(samples):

    # Extract barcode list from sample dictionary
    barcode_list = []
    for keys in samples:
        barcode_list.append(keys)

    barcode_length = len(barcode_list[0])

    print('\n=== Checking barcodes ===')

    for b in barcode_list:
        print(b)
        if len(b) != barcode_length:
            print('Error: barcodes are not the same length')
            exit() # How to really do error handling???
    print('n={0} barcodes of same length ({1} nt)'.format(len(barcode_list), barcode_length))

    return barcode_list, barcode_length



def demultiplex_fastq(fastq_file):
    if fastq_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    sample_list = { \
    'CGAT':'Input1', \
    'GCTA':'Input2', \
    'AGTC':'Output1', \
    'AAAA':'Output2'}

    # barcodes = list of barcodes (sequences only)
    # b_len = barcode length
    barcodes, b_len = barcodes_prep(sample_list)

    with opener(fastq_file, 'r') as f:

        print('\n=== Demultiplexing FASTQ input file by 5\' barcode ===')

        (file_root, file_ext) = (os.path.splitext(fastq_file))

        # fastq record
        # identifier =     record[0]
        # sequence =       record[1]
        # alt_identifier = record[2]
        # quality =        record[3]
        record = []
        for i,line in enumerate(f):
            record.append(line.rstrip('\n'))
            if i % 4 == 3:
                if record[1][0:b_len] in barcodes:
                    with open('{0}_{1}{2}'.format(file_root, record[1][0:b_len], file_ext), 'a') as fo:
                        fo.write('\n'.join(record) + '\n')
                record = []
            if (i+1) % 1E+5 == 0:
                if (i+1) % 1E+6 == 0:
                    print('\n=== Demultiplexing FASTQ input file by 5\' barcode ===')
                print('{0} records processed.'.format(i+1)) # index i starts at 0
        print('{0} total records processed.'.format(i+1))


def main():
    reads = 'E689_400k_lines.fastq'
    demultiplex_fastq(reads)

if __name__ == "__main__":
    main()
