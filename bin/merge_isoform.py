'''
Take corrected reads and perform the following steps:
Step 1. get first pass isoforms
Step 2. filter and merge
'''

import pandas as pd
import numpy as np
import os
import pandas as pd
import seaborn as sns
from collections import defaultdict
from tqdm import tqdm


def get_1_pass_isoform(read_df):
    """
    Step 1. get first pass isoforms
        1.1 count how many reads have exactly same intron chain
        1.2 determine the tss tts for each isoform (mode of the read tss/tts)
    """
    tss = []
    tts = []
    junc = []
    reference_name = []
    transcript_strand = []
    count = []
    grp_iter = read_df.groupby(by=['reference_name','junc']).__iter__()

    for _,grp in grp_iter:
        start_mode = grp.trans_start[grp.intron_removed_start == 0].mode()
        if len(start_mode):
            tss.append(min(start_mode))
        else:
            tss.append(np.nan)

        end_mode = grp.trans_end[grp.intron_removed_end == 0].mode()
        if len(end_mode):
            tts.append(max(end_mode))
        else:
            tts.append(np.nan)
        junc.append(grp.junc.iloc[0])
        reference_name.append(grp.reference_name.iloc[0])
        transcript_strand.append(grp.transcript_strand.iloc[0])
        count.append(len(grp))

    isoform_df = pd.DataFrame({
        'reference_name':reference_name,
        'transcript_strand':transcript_strand,
        'tss':tss,
        'junc':junc,
        'tts':tts,
        'counts':count
    })
    return isoform_df


def merge_transcript(df, min_count_to_collapes,max_ts_dff=100, max_exbound_ins=10,merge_count=True):
    '''
    Step 2. filter and merge
    rule:
    1. the parent isoform must have >= `min_count_to_collapes`
    2. the tss and tts of the child isoform must not exceed  
        the corresponding exon boundary of the parant isoform
        more than `max_exbound_ins`
    3. the tss and tts of the child isoform must not exceed  
        the corresponding tss and tts of the parant isoform
        more than `max_ts_dff`
        '''

    def merge(row, merged_d, max_ts_dff, max_exbound_ins, 
                append_if_not_exist=True, merge_count=True, check_count = False):
        """
        check if an intron_chain in row is partailly contained in merged d
            if no: create a new entry in merged_d if no
            if yes: add the read count to the corresponding isoform in merge_d 
        """
        
        if not len(merged_d):
            merged_d = pd.DataFrame([row]).set_index('Index')
            return merged_d
        
        chain_len = len(row.junc)
        merge_into = []
        for x in merged_d[(merged_d.reference_name==row.reference_name) & 
                        (merged_d.transcript_strand==row.transcript_strand)].itertuples():
            if len(x.junc) <= len(row.junc):
                continue
            
            if check_count and x.counts < 0.2*row.counts: # do not merge if the longer transcript has lower count
                continue
            
            for i in range(len(x.junc)-chain_len+1):
                if x.junc[i:i+chain_len] == row.junc:
                    if ~np.isnan(row.tss) and i==0 and abs(row.tss-x.tss)>max_ts_dff:
                        continue
                    if ~np.isnan(row.tss) and i!=0 and x.junc[i-1][1] - row.tss > max_exbound_ins:
                        continue
                    if ~np.isnan(row.tts) and i==len(x.junc)-chain_len and abs(row.tts-x.tts)>max_ts_dff:
                        continue
                    if ~np.isnan(row.tts) and i!=len(x.junc)-chain_len and row.tts - x.junc[i+chain_len][0] > max_exbound_ins:
                        continue
                    # found merge: update count
                    merge_into.append(x)
        
        if len(merge_into) and merge_count:
            parent_sum = sum([x.counts for x in merge_into])
            for x in merge_into:
                merged_d.at[x.Index, 'counts'] += row.counts*x.counts/parent_sum
        if not len(merge_into):
            if append_if_not_exist and ~np.isnan(row.tss) and ~np.isnan(row.tts):
                merged_d = pd.concat([merged_d, pd.DataFrame([row]).set_index('Index')])
        
        return merged_d
    
    df['junc_count'] = df.junc.apply(len)
    df.sort_values(['reference_name','transcript_strand', 'junc_count'], 
                                ascending = [True,True,False], inplace = True)    
    
    # dict for confirmed long transcript e.g. {(ch1, '+'): row in df}
    merged_d = pd.DataFrame()
    for row in tqdm(df[df.counts >= min_count_to_collapes].itertuples()):
        merged_d = merge(row, merged_d, max_ts_dff=max_ts_dff, 
            max_exbound_ins=max_exbound_ins, append_if_not_exist=True, 
            merge_count=merge_count)
    
    if merge_count: # merge count of 1 pass isoforms that do not reach thres
        for row in tqdm(df[df.counts < min_count_to_collapes].itertuples()):
            merged_d = merge(row, merged_d, max_ts_dff=max_ts_dff, max_exbound_ins=max_exbound_ins, append_if_not_exist=False)
    
    merged_d.drop(columns = ['junc_count'], inplace =True)
    return merged_d



