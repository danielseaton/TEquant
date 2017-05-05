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
TE_mapping_dicts = dict()
levels=['transcript_id','gene_id','family_id','class_id']
for idx in range(len(levels)-1):
    current_level = levels[idx]
    next_level = levels[idx+1]
    temp_TE_df = TE_df[[current_level,next_level]].drop_duplicates()
    temp_TE_df.set_index(current_level,drop=True,inplace=True)
    TE_mapping_dicts[current_level] = temp_TE_df[next_level].to_dict()


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
TE_ids = list(set(TE_df['transcript_id'])&set(eDF.index))
eDF = eDF.loc[TE_ids+footer_ids,:]

#Output to file
output_filename = '{}_TE{}_counts.tsv'.format(args.output_prefix,'transcript_id')
eDF.to_csv(output_filename,sep='\t')

levels = ['gene_id','family_id','class_id'] #Must be in increasing level in heirarchy
previous_level = 'transcript_id' #start with the most basic element
for level in levels:
    # Make an index for the TE dataframe at this level of description
    new_TE_index = list(set(TE_df[level])) + footer_ids
    output_eDF = pd.DataFrame(0,columns=eDF.columns,index=new_TE_index)
    output_eDF.loc[footer_ids,:] = eDF.loc[footer_ids,:]

    # Each row from the lower level gets added to one row of the next level up
    for TE_id in eDF.index:
        if TE_id not in footer_ids:
            TE_id_to_add_to = TE_mapping_dicts[previous_level][TE_id]
            output_eDF.loc[TE_id_to_add_to,:] += eDF.loc[TE_id,:]

    output_filename = '{}_TE{}_counts.tsv'.format(args.output_prefix,level)
    output_eDF.to_csv(output_filename,sep='\t')
    # Finish by resetting eDF to this output
    # This reduces the number of summations required
    eDF = output_eDF
    previous_level = level

