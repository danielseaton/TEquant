import pandas as pd
import numpy as np

def merge_featurecounts_by_dataframe(input_filenames,summary_filenames,output_prefix,filterzeros):
    '''Merge outputs from multiple featurecounts output files.'''

    assert(len(raw_count_filenames)==len(summary_filenames))

    #Concatenate featurecounts files
    list_of_dfs = []

    for idx,filename in enumerate(raw_count_filenames):
        df_temp = pd.read_csv(filename,sep='\t',index_col=0,header=1)
        if idx==0:
            #Feature dataframe should be the same across all samples
            fDF=df_temp.ix[:,:5]
            #Basic check on file format
            assert(len(df_temp.columns)==6)
            #Select only the counts, which are in the last column (which has a variable column name)
            df = df_temp.ix[:,-1]
            list_of_dfs.append(df)

    eDF = pd.concat(list_of_dfs,axis=1)

    if filterzeros:
        eDF = eDF.loc[~(eDF==0).all(axis=1)]

    #Parse featurecounts summary file to get total numbers of aligned reads used in each case
    list_of_dfs = []
    for filename in summary_filenames:
        df = pd.read_csv(filename,sep='\t',index_col=0)
        list_of_dfs.append(df)
        summary_df = pd.concat(list_of_dfs,axis=1)

    summary_df.index = pd.Index(['_'+x.lower() for x in summary_df.index])

    eDF = pd.concat([eDF,summary_df])

    #Check that the indices line up correctly (i.e. that summary and count files match up)
    assert(len(summary_df.columns)==len(eDF.columns))
    assert(len(set(summary_df.columns)&set(eDF.columns))==len(eDF.columns))

    fDF.to_csv(output_prefix+'_features.tsv',sep='\t')
    eDF.to_csv(output_prefix+'_counts.tsv',sep='\t')

    return 0
