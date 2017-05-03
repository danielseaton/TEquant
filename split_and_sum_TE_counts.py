import argparse
import fileinput
import pandas as pd
import numpy as np
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i",
                    type=str,
                    default='',
                    dest='input',
                    action='store',
                    help='Path to input file (tab-separated).'
                    )
parser.add_argument("-o",
                    type=str,
                    default='',
                    dest='output_prefix',
                    action='store',
                    help='Output prefix'
                    )
parser.add_argument("-m",
                    type=str,
                    default='',
                    dest='TE_map_file',
                    action='store',
                    help='TE name mapping'
                    )
args = parser.parse_args()


#Set a pattern used to identify gene identifiers, so they can be separated from TEs
gene_id_regex_pattern = "ENS[A-Z]+[0-9]{11}"
gene_id_regex = re.compile(gene_id_regex_pattern)
footer_id_regex = re.compile("_[a-z]+")

#Input identification and classification information for TE identifiers
TE_df = pd.read_csv(args.TE_map_file,sep='\t')
TE_df.head()

#Load input expression data
eDF = pd.read_csv(args.input,sep='\t',index_col=0)

#Put gene data in one output file
gene_ids = []
footer_ids = []
for x in eDF.index:
    if gene_id_regex.search(x):
        gene_ids.append(x)
    if footer_id_regex.match(x):
        footer_ids.append(x)
assert(len(footer_ids)<15)

eDF.loc[gene_ids+footer_ids,:].to_csv(args.output_prefix+'_gene_id_counts.tsv',sep='\t')

#Cut down matrix to just TE elements
TE_ids = TE_df['transcript_id'].to_list()
eDF = eDF.loc[TE_ids+footer_ids,:]

#Output to file
output_filename = '{}_TE{}_counts.tsv'.format(args.output_prefix,'transcript_id')
eDF.to_csv(output_filename,sep='\t')

levels = ['gene_id','family_id','class_id'] #MUST be in increasing level in heirarchy
previous_level = 'transcript_id' #start with the most basic element
for level in levels:
    output_eDF = pd.DataFrame(0,columns=eDF.columns,index=IDs)
    for TE_id in eDF.index:
        TE_id_to_add_to = TE_df[TE_df[previous_level]==TE_id][level]
        output_eDF.loc[TE_id_to_add_to,:] = eDF.loc[TE_id,:]

    output_filename = '{}_TE{}_counts.tsv'.format(args.output_prefix,level)
    output_eDF.to_csv(output_filename,sep='\t')
    # Finish by resetting eDF to this output
    # This reduces the number of summations required
    eDF = output_eDF
    previous_level = level