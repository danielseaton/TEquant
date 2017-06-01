
# coding: utf-8

# In[35]:

import pandas as pd
import json


# In[33]:

gtf = '/nfs/leia/research/stegle/dseaton/genomes/hg19/annotation/hg19_rmsk_TE.gtf'

with open(gtf,'r') as f:
    #get part of line with gene identifiers
    lines = [x.strip().split('\t')[8].strip(';') for x in f.readlines()]
    lines = [x.replace('"','') for x in lines]
    #split lines
    split_lines = [x.split(';') for x in lines]
    name_tuples = [[tuple(y.strip().split(' ')) for y in x] for x in split_lines]
    #dicts
    gtf_dicts = [dict(x) for x in name_tuples]


# In[34]:

#Map from the class id through the heirarchies down to transcript_id (the individual elements)

levels = ['class_id','family_id','gene_id','transcript_id']
TE_dict = dict()

for d in gtf_dicts:
    TE_dict[d['class_id']] = dict()

for d in gtf_dicts:
    TE_dict[d['class_id']][d['family_id']] = dict()

for d in gtf_dicts:
    TE_dict[d['class_id']][d['family_id']][d['gene_id']] = []

for d in gtf_dicts:
    TE_dict[d['class_id']][d['family_id']][d['gene_id']].append(d['transcript_id'])


# In[36]:

#Output TE dict as a json file

with open('TE_class_to_id_dict.json', 'w') as f:
    json.dump(TE_dict, f)


# In[37]:

#Output TE info as a table

TEdf = pd.DataFrame(gtf_dicts)
TEdf = TEdf[['transcript_id','gene_id','family_id','class_id']]

TEdf.to_csv('TE_name_mapping.tsv',sep='\t',index=False)


# In[ ]:



