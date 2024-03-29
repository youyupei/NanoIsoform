import argparse
import textwrap
import pandas as pd
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm
import numpy as np
import logging


import helper
from helper import add_summary
from arg_parser import *

logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def group_reads(all_reads, max_diff = GROUP_ARG['max_diff']):
    """group reads based on their junctions
    Args:
        all_reads (pandas.DataFrame): The dataframe output by restructure_per_jwr_datafrome
        function.
        max_diff (INT): maximum different allowed when grouping closeby junction

    Output:
        add "subgroup" columns to the dataframe, each subgroup of reads share 
        exactly same junctions
    """
    def split_groupby(df, max_diff):
        org_d = df.copy()

        # get junction np.array: [[],[],[],...]
        temp_d = df.index.unique().to_frame().copy()
        junc_array = temp_d.junc.apply(lambda x: np.array(x).flatten())
        junc_array = np.vstack(junc_array)

        # junc_array: a columns is a splice junc
        group_by_sites = np.vstack([get_groups(col, max_diff) for col in junc_array.T]).T
        group = np.unique(group_by_sites, axis=0, return_inverse=True)[1]
        temp_d['sub_group'] = group
        merge_d = org_d.merge(temp_d, 'left', left_index=True, right_index=True)
        merge_d = merge_d.drop(columns = ['reference_name',
                                        'non_overlap_group', \
                                        'junc_count', 
                                        'junc'])
        return merge_d
    
    def get_groups(array, max_diff = max_diff):
        # add index
        array_indexed = np.vstack([np.arange(len(array)), array]).T
        # sort 
        array_indexed = array_indexed[array_indexed[:, 1].argsort()]
        # get group
        rst_group = __get_groups_from_sorted_array(array_indexed[:,1])
        # revert oder
        return rst_group[array_indexed[:, 0].argsort()]

    def __get_groups_from_sorted_array(array, max_diff=max_diff):
            # note that array has to be sorted
            # padding
            array_pad = np.hstack([array[0], array])
            return np.cumsum(np.abs(array_pad[1:] - array_pad[:-1]) > max_diff)
  
    def add_non_overlap_group(df):
        """Output a copy of the df with a new columns "non_overlap_group"
        The reads in the df should mapped to a some chromosom.
        """
        def __check_non_overlap(junc_col):
            """Take a pd.Series (formatted the same as columns 'junc') for same 
            reference name and output a non-overlapping group label with same length
            """
            region_span = list(junc_col.apply(lambda x: (x[0][0], x[-1][1])))
            region_span = np.array(region_span)
            # structure: [[index, site1, site2],...]
            region_span = \
                np.hstack(
                    [np.arange(region_span.shape[0]).reshape(-1,1), 
                    region_span])
            
            # sort by junction start
            region_span = region_span[region_span[:, 1].argsort()]

            group_label = -1
            group_end_site = -1
            groups = []
            for i in region_span:
                if i[1] >= group_end_site:
                    group_end_site = i[2]
                    group_label += 1
                    groups.append(group_label)
                else:
                    groups.append(group_label)
                    group_end_site = max(group_end_site, i[2])
            return np.array(groups)[region_span[:, 0].argsort()]

        out_d = df.copy()
        out_d['non_overlap_group'] = out_d.junc.pipe(__check_non_overlap)
        return out_d

    #################     
    all_reads = all_reads.groupby(by=['reference_name'
                                    ]).apply(add_non_overlap_group)
    all_reads.reset_index(drop=True, inplace=True)
    all_reads.set_index(['reference_name',
                        'non_overlap_group',
                        'junc_count', 
                        'junc'], inplace=True) 
    
    df_list = []
    for k, d in all_reads.groupby(by=['reference_name',
                                        'transcript_strand',
                                            'non_overlap_group',
                                            'junc_count']).__iter__():
        df_list.append(split_groupby(d, max_diff))
    del all_reads
    all_reads = pd.concat(df_list)
    
    return all_reads
   ############################3 

    # all_reads = all_reads.set_index(['reference_name','junc_count', 'junc'])   

    # # all reads group by chr/number of junctions/ junctions()
    # df_list = []
    # for k, d in all_reads.groupby(level=[0,1]).__iter__():
    #     df_list.append(split_groupby(d, max_diff))
    # del all_reads
    
    # all_reads = pd.concat(df_list)

    # if all_reads.index.nlevels > 3:
    #     # sometime, after merging, there are duplicated columns. Not sure about
    #     # the reason
    #     all_reads = all_reads.droplevel([0,1])

    # return all_reads

'''
Correct junction for each junction in each reads group
For each junction:
    methods 1 (majority_vote): Look for nearby mapped junctions with corrected identification,
                and correct uncorrected junction based on majoirity vote
    methods 2 (probability): Look for nearby mapped junctions with corrected identification,
                and correct uncorrected junction based on probability 
'''
def correct_junction_per_group(all_reads, methods):
    """Correct Junction in each group generated in group_reads
    Args:
        all_reads: pd.DataFrame (will reset index inside the input df)
        method: choose from 'majority_vote' and 'probability'
    """
    
    # work on group_single, which is a copy 
    def correct_line(df_line, correct_dict, methods=methods):
        """Read a single line of the dataframe
        return a currected list of junction for each read
        """
        if methods == 'majority_vote':
            return \
            tuple([x if y else correct_dict[x] for x,y in zip(
                df_line.junc, df_line.corrected)])
        if methods == 'probability':
            return \
            tuple([x if y else next(correct_dict[x]) for x,y in zip(
                df_line.junc, df_line.corrected)])
               
    def generate_df_per_junction(df, num_of_junc):
        """
        Process each df and key in pandas.groupby.__iter__() output

        Yields:
            Counter of junctions (count the junction HCJWR support)
            Unique list of uncorrect junction
        """
        df = df[['junc', 'corrected','high_conf']].copy()
        for i in range(num_of_junc):
            junc = df['junc'].apply(lambda x: x[i])
            corrected = df['corrected'].apply(lambda x: x[i])

            high_conf = df['high_conf'].apply(lambda x: x[i])
            print(df['high_conf'])
            exit()
            # yield count of corrected junc and a list of uncorrected junc
            yield Counter(junc[high_conf]), junc[~high_conf].unique()

    def correct_junction_majority_vote(counter, junc, 
                                    max_diff, min_prop, min_correct_read):
        """Correct uncorrected JWRs based on nearby (with maxdiff) corrected JWRs.
        Each JWR will be corret to the major junction supported by the nearby corrected
        JWRs

        Args:
            counter (collection.Counter): Number of corrected JWRs supporting each 
                                            nearby junction
            junc (pd.Series): Unique uncorrected junctions 
            max_diff (int): Maximal distance to define "nearby junctions"
            min_prop (float): Minimum proportion of the major junction to trigger a correction
            min_correct_read (int): Minimum number of corrected reads for the major junction
        Return:
            corrected location for junction (`junc`) or 
            None if no near by junction reach the min_prop.
        """
        def check_max_diff(junc1, junc2, max_diff=max_diff):
            return np.all(np.abs([x-y for x,y in zip(junc1, junc2)]) <= max_diff)
        keys = [k for k in counter.keys() if check_max_diff(k, junc)]

        # return None if no nearby corrected junction
        if not len(keys):
            return None

        values = np.array([counter[k] for k in keys])
        prop = values/sum(values)
        if prop.max() >= min_prop and values.max() >= min_correct_read:
            return keys[np.random.choice(np.array(range(len(keys)))[prop == prop.max()])]
        else:
            return None

    def correct_junction_probability(
                            counter, junc, 
                            max_diff=CORRECTION_ARG['dist'],
                            min_prop=CORRECTION_ARG['prob_samp_min_prop'], 
                            min_correct_read=CORRECTION_ARG['prob_samp_min_count']):
        """Correct uncorrected JWRs based on nearby (with maxdiff) corrected JWRs.
        Each JWR will be corret to the nearby junctions based on probabilities (cal-
        culated based on the proportion of nearby JWRs supporting each junction)
        Args:
            counter (collection.Counter): Number of corrected JWRs supporting each 
                                            nearby junction
            junc (pd.Series): Unique uncorrected junctions 
            max_diff (int): Maximal distance to define "nearby junctions"
            min_prop (float): Minimum proportion of the each junction to be a candidate
            min_correct_read (int): Minimum number of corrected reads for the each junction
        Return:
            corrected location for junction (`junc`) or 
            None if no near by junction reach the min_prop.
        """
        def check_max_diff(junc1, junc2, max_diff=max_diff):
            return np.all(np.abs([x-y for x,y in zip(junc1, junc2)]) <= max_diff)
        
        keys = [k for k in counter.keys() if check_max_diff(k, junc)]

        # return None if no nearby corrected junction
        if not len(keys):
            while True:
                yield None
        values = np.array([counter[k] for k in keys])
        prop = values/sum(values)
        pass_key_idx = (prop >= min_prop) & (values >= min_correct_read)
        prop = prop[pass_key_idx]
        values = values[pass_key_idx]
        keys = [x for x,y in zip(keys, pass_key_idx) if y]

        if len(prop):
            while True:
                yield keys[np.random.choice(range(len(keys)), p = prop/sum(prop))]
        else:
            while True:
                yield None


    def correct_junction_HCJWR(
                            counter, junc, 
                            max_diff=CORRECTION_ARG['dist'],
                            min_hcjwr=HCJWR['min_count']):
        """Re-correct all JWRs, including ones with NanoSplicer output based on 
        nearby (with maxdiff) High-confidence JWRs (HCJWRs).
        
        Args:
            counter (collection.Counter): Number of HCJWRs supporting each 
                                            nearby junction
            junc (pd.Series): Unique junctions fron non-HC JWRs
            max_diff (int): Maximal distance to define "nearby junctions"
            min_hcjwr (int): Minimum number of HCJWR for the each junction
        Return:
            corrected location for junction (`junc`) or 
            None if no near by junction reach the min_prop.
        """
        def check_max_diff(junc1, junc2, max_diff=max_diff):
            return np.all(np.abs([x-y for x,y in zip(junc1, junc2)]) <= max_diff)
        
        keys = [k for k in counter.keys() if check_max_diff(k, junc)]

        # return None if no nearby HC junction
        if not len(keys):
            while True:
                yield None
        values = np.array([counter[k] for k in keys])
        prop = values/sum(values)
        pass_key_idx = (prop >= min_prop) & (values >= min_correct_read)
        prop = prop[pass_key_idx]
        values = values[pass_key_idx]
        keys = [x for x,y in zip(keys, pass_key_idx) if y]

        if len(prop):
            while True:
                yield keys[np.random.choice(range(len(keys)), p = prop/sum(prop))]
        else:
            while True:
                yield None

    def small_groups_jwr_recover(group_single, all_junc_counter, methods):
        """Recovering the uncorrected JWRs using cross-group information
        Args:
            df (pd.DataFrame): DataFrame containing reads corrected by
                * NanoSplicer
                * JAQ-based approach
                * Nearby-JWR-based correct in group with sufficient reads
            corrected_counter (collections.Counter): Counts of unique junctions 
                Supported by:
                    * JWRs corrected by NanoSplicer
                    * JAQ-based approach
        Return Df
            """
        correct_dict = {}
        junc_count = len(list(group_single.junc)[0])
        
        for counter, junc in generate_df_per_junction(group_single, junc_count):
            all_junc_counter += counter

            for j in junc:
                if methods == 'majority_vote':
                    corrected_junc = correct_junction_majority_vote(
                        counter, j, 
                        max_diff=CORRECTION_ARG['dist'],
                        min_prop=CORRECTION_ARG['maj_vot_min_prop'], 
                        min_correct_read=CORRECTION_ARG['maj_vot_min_count'])
                    
                    correct_dict[j] = corrected_junc

                elif methods == 'probability':
                    corrected_junc = correct_junction_probability(
                        counter, j, 
                        max_diff=CORRECTION_ARG['dist'],
                        min_prop=CORRECTION_ARG['prob_samp_min_prop'], 
                        min_correct_read=CORRECTION_ARG['prob_samp_min_count'])   
                    correct_dict[j] = corrected_junc

        group_single['corrected_junction'] = \
            group_single.apply(correct_line, axis=1, correct_dict= correct_dict) 
        group_single['all_corrected'] = \
            group_single['corrected_junction'].apply(all)
        
        return group_single



    all_reads.reset_index(drop=False, inplace = True)

    # count read in group (JAQ is just a random colname)
    # all_reads['group_count'] = all_reads.groupby(
    #                                         ['reference_name',
    #                                         'non_overlap_group',
    #                                         'junc_count',
    #                                         'sub_group']).JAQ.transform('count')
    
    groups_gb_region = all_reads.groupby(['reference_name', 
                                        'non_overlap_group'])
    
    corrected_group = []
    
    # iter over non-overlap region
    for _, d_region in tqdm(groups_gb_region.__iter__(), 
                            total = len(groups_gb_region.groups.keys()),
                            desc='Genome region'):

        all_junc_counter = Counter()
        uncorrected_group = []
        sub_groups_gb = d_region.groupby(['junc_count', 'sub_group'])
        
        # iter over sub_group (same number of junctions with similar location)
        for k,d in tqdm(sub_groups_gb.__iter__(), 
                        leave=False,
                        total = len(sub_groups_gb.groups.keys())):
            # k: (junc_count,sub_group)
            # d:'index', 'reference_name', 'non_overlap_group', 'junc_count', 'junc',
                # 'read_id', 'transcript_strand', 'junc_start', 'junc_end', 'corrected',
                # 'JAQ', 'sub_group'
            group_single = d.copy()

            # dict for correcting uncorrected jwrs
            correct_dict = {}
            junc_count = k[0]
            for counter, junc in generate_df_per_junction(group_single, junc_count):
                all_junc_counter += counter

                for j in junc:
                    corrected_junc = correct_junction_HCJWR(counter, j)
                    correct_dict[j] = corrected_junc


            group_single['corrected_junction'] = \
                group_single.apply(correct_line, axis=1, correct_dict= correct_dict) 
            
            group_single['corrected'] = \
                group_single.corrected_junction.apply(
                    lambda y: tuple([bool(x) for x in y])) 
            
            group_single['all_corrected'] = \
                group_single['corrected_junction'].apply(all)

            if all(group_single['all_corrected']):
                corrected_group.append(group_single)
            else:    
                uncorrected_group.append((group_single, counter))

###########################
        # deal with uncorrected_group all_junc_counter here and add it to corrected_group
        for group_single, counter in tqdm(uncorrected_group):
            corrected_group.append(small_groups_jwr_recover(group_single, 
                                    all_junc_counter,
                                    methods=methods))
###############

    corrected_d = pd.concat(corrected_group)
    corrected_d.drop(columns=['junc_start', 'junc_end', 'corrected'], inplace=True)
    
    # output columns:
        # index 
        # reference_name  
        # junc_count 
        # junc 
        # read_id 
        # transcript_strand 
        # junc_start* (dropped in corrected_d) 
        # junc_end*    
        # corrected*  
        # JAQ
        # sub_group  
        # group_count 
        # corrected_junction  
        # all_corrected  
    # small_groups_jwr_recover(uncorrected_group[0], all_junc_counter)
    # exit()
    return corrected_d





    # # dic structure: group_key:Counter

    # for key in tqdm(groups_gb_obj.groups.keys()): # can potentially do multiprocessing
    #     # dictionary for correction
    #     group_single = groups_gb_obj.get_group(key).copy()
        
    #     # dict for correcting uncorrected jwrs
    #     correct_dict = {}
    #     for counter, junc in generate_df_per_junction(group_single, key):#key:num_of_junc
    #         # save counts of certain subgroup
    #         all_junc_counter += counter

    #         for j in junc:
    #             if methods == 'majority_vote':
    #                 corrected_junc = correct_junction_majority_vote(
    #                     counter, j, 
    #                     max_diff=CORRECTION_ARG['dist'],
    #                     min_prop=CORRECTION_ARG['maj_vot_min_prop'], 
    #                     min_correct_read=CORRECTION_ARG['maj_vot_min_count'])
                    
    #                 correct_dict[j] = corrected_junc

    #             elif methods == 'probability':
    #                 corrected_junc = correct_junction_probability(
    #                     counter, j, 
    #                     max_diff=CORRECTION_ARG['dist'],
    #                     min_prop=CORRECTION_ARG['prob_samp_min_prop'], 
    #                     min_correct_read=CORRECTION_ARG['prob_samp_min_count'])   
    #                 correct_dict[j] = corrected_junc

    #     group_single['corrected_junction'] = \
    #         group_single.apply(correct_line, axis=1, correct_dict= correct_dict) 
    #     group_single['all_corrected'] = \
    #         group_single['corrected_junction'].apply(all)
    #     corrected_group.append(group_single)
    #     if all(group_single['all_corrected']):
    #         corrected_group.append(group_single)
    #     else:    
    #         uncorrected_group.append(group_single)


    # corrected_d = pd.concat(corrected_group)
    # corrected_d.drop(columns=['junc_start', 'junc_end', 'corrected'], inplace=True)
    
    # # output columns:
    #     # index 
    #     # reference_name  
    #     # junc_count 
    #     # junc 
    #     # read_id 
    #     # transcript_strand 
    #     # junc_start* (dropped in corrected_d) 
    #     # junc_end*    
    #     # corrected*  
    #     # JAQ
    #     # sub_group  
    #     # group_count 
    #     # corrected_junction  
    #     # all_corrected  
    # # small_groups_jwr_recover(uncorrected_group[0], all_junc_counter)
    # # exit()
    # return corrected_d, uncorrected_group, all_junc_counter

