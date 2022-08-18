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
from config import *

logging.basicConfig(format=LOG_FORMAT, datefmt='%a, %d %b %Y %H:%M:%S')
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)



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
    def correct_line(df_line, correct_dict):
        """Read a single line of the dataframe
        return a currected list of junction for each read
        """
        return \
        tuple([x if y else correct_dict[x] for x,y in zip(
            df_line.junc, df_line.corrected)])

    all_reads.reset_index(drop=False, inplace = True)
    # count read in group

    all_reads['group_count'] = all_reads.groupby(
        ['reference_name','junc_count','sub_group']).JAQ.transform('count')
    

    groups_gb_obj = all_reads.groupby(['reference_name', 'junc_count','sub_group', 'group_count'])
    # determine the order in group. Group with more reads comes first.
    keys_ordered = sorted(groups_gb_obj.groups.keys(), key = lambda x: x[3], reverse =True)

    # dic structure: group_key:Counter
    all_junc_counter = defaultdict(Counter)
    corrected_group = []
    for key in tqdm(keys_ordered): # can potentially do multiprocessing
        # dictionary for correction
        group_single = groups_gb_obj.get_group(key).copy()
        
        # dict for correcting uncorrected jwrs
        correct_dict = {}
        for counter, junc in generate_df_per_junction(group_single, key):
            # save counts of certain subgroup
            all_junc_counter[key] += counter

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
        corrected_group.append(group_single)
    corrected_d = pd.concat(corrected_group)
    corrected_d.drop(columns=['junc_start', 'junc_end', 'corrected'], inplace=True)
    return corrected_d

def generate_df_per_junction(df, key):
    """
    Process each df and key in pandas.groupby.__iter__() output

    Yields:
        Counter of junctions (count the corrected junction support)
        Unique list of uncorrect junction
    """
    df = df[['junc', 'corrected']].copy()
    for i in range(key[1]):
        junc = df['junc'].apply(lambda x: x[i])
        corrected = df['corrected'].apply(lambda x: x[i])
        # yield count of corrected junc and a list of uncorrected junc
        yield Counter(junc[corrected]), junc[~corrected].unique()

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
        return None

    values = np.array([counter[k] for k in keys])
    prop = values/sum(values)
    pass_key_idx = (prop >= min_prop) & (values >= min_correct_read)
    prop = prop[pass_key_idx]
    values = values[pass_key_idx]
    keys = [x for x,y in zip(keys, pass_key_idx) if y]

    if len(prop):
        return keys[np.random.choice(range(len(keys)), p = prop/sum(prop))]
    else:
        return None