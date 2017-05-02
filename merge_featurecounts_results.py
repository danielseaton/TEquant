import pandas as pd
import numpy as np
import glob
import os
import argparse
import fileinput

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
                    dest='output_filename',
                    action='store',
                    help='Path of output file.'
                    )
parser.add_argument("--nofilter",
                    action="store_true",
                    help="If specified, don't filter out rows with zero counts across all samples")

args = parser.parse_args()

#Concatenate featurecounts files
list_of_dfs = []
for filename in args.raw_count_filenames:
    df = pd.read_csv(filename,sep='\t',index_col=0,header=1)
    #Select only the counts, which are in the last column (which has a variable column name)
    df = df.ix[:,-1]
    list_of_dfs.append(df)

eDF = pd.concat(list_of_dfs,axis=1)

#Parse featurecounts summary file to get total numbers of aligned reads used in each case
list_of_dfs = []
for filename in args.summary_filenames:
    df = pd.read_csv(filename,sep='\t',index_col=0)
    list_of_dfs.append(df)
summary_df = pd.concat(list_of_dfs,axis=1)

rows_to_sum = ['Assigned','Unassigned_Ambiguity','Unassigned_NoFeatures']
total_mapped_reads = summary_df.loc[rows_to_sum,:].apply(sum)

eDF.loc['_total_mapped_reads']=total_mapped_reads
assert(len(total_mapped_reads.index)==len(eDF.columns))
assert(len(set(total_mapped_reads.index)&set(eDF.columns))==len(total_mapped_reads.index))

if not args.nofilter:
    eDF = eDF.loc[~(eDF==0).all(axis=1)]

eDF.to_csv(args.output_filename,sep='\t')
