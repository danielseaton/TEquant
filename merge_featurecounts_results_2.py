import pandas as pd
import numpy as np
import glob
import os
import argparse
import fileinput

'''Merge outputs from multiple featurecounts output files.

Example usage for merging outputs for all featurecount output files in a results directory:

"python merge_featurecounts_results.py -i /nfs/leia/research/stegle/dseaton/hipsci/TEquantification/data/featurecounts_all/*/*.txt -s /nfs/leia/research/stegle/dseaton/hipsci/TEquantification/data/featurecounts_all/*/*.summary -o outfile.txt"

'''

parser = argparse.ArgumentParser(description='A script to concatenate output from featurecounts.')
parser.add_argument("-i",
                    type=str,
                    nargs='+',
                    default='',
                    dest='raw_count_filenames',
                    action='store',
                    help='Path to input count files.'
                    )
parser.add_argument("-s",
                    type=str,
                    nargs='+',
                    default='',
                    dest='summary_filenames',
                    action='store',
                    help='Path to input summary files.'
                    )
parser.add_argument("-o",
                    type=str,
                    default='',
                    dest='output_prefix',
                    action='store',
                    help='Output file prefix. Output files will be {output_prefix}_counts.tsv and {output_prefix}_features.tsv'
                    )
parser.add_argument("--filterzeros",
                    action="store_true",
                    help="If specified, filter out rows with zero counts across all samples")

args = parser.parse_args()

assert(len(args.raw_count_filenames)==len(args.summary_filenames))


output_file = open(args.output_prefix+'_counts.tsv','w')

nF = len(args.raw_count_filenames)
file_iterators = [open(f,'r') for f in args.raw_count_filenames]
#First line is thrown away
first_lines = [x.next() for x in file_iterators]
#Second line gives info for the header
header_lines = [x.next() for x in file_iterators]
output_header = 'GeneID\t'+'\t'.join([x.strip().split('\t')[6] for x in header_lines]) + '\n'
output_file.write(output_header)
while True:
    try:
        next_lines = [x.next() for x in file_iterators]
        values = [x.strip().split('\t')[6] for x in next_lines]
        if args.filterzeros:
            if all([x=='0' for x in values]):
                continue
        gene_identifier = next_lines[0].strip().split('\t')[0]
        line = gene_identifier +'\t' + '\t'.join(values) + '\n'
        output_file.write(line)
    except StopIteration:
        #files are empty
        break
output_file.close()
