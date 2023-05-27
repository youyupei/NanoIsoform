"""
JWR level analysis:
Take as input the NanoSplicer result and JWR_checker result:
1. Filter NanoSplicer result to remove low SIQ and low porbability JWR
2. Identify a list of HCJWR and count
3. Add NanoSplicer output to JWR checker result (some of the JWR do not have NanoSplicer result will be NA)
4. group junction based on initial mapping
5. correct each junction to HCJWR nearby
"""

import argparse
import textwrap
import pandas as pd
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm
import numpy as np
import logging
import swifter
import itertools

import helper
from helper import add_summary
# import configs
from arg_parser import *

logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def parse_nanosplicer_prob_table(args):
    """Parse the prob table file output by NanoSplicer as 
    Args:
        tb_fn (str): NanoSplicer file name
        correct_to_hcjwr (bool): whether or not correct all junction to HCJWR
    return:
        pandas data from
    """
    def format_initial_junction(s):
        x1,x2 = s.strip('[()]').split(',')
        return (int(x1), int(x2))

    def format_candidates(s):
        s = [x.strip('[()]').split(',') for x in s.split('),(')]
        return  tuple([(int(x1), int(x2)) for x1, x2 in s])
    def format_tuple_string(s):
        return tuple([float(x) for x in s.split(',')])
    d = pd.read_csv(args.prob_table, sep='\t')
    d['initial_junction'] = d.initial_junction.swifter.apply(format_initial_junction)
    d['candidates'] = d.candidates.swifter.apply(format_candidates)
    
    if args.uniform_prior:
        d['prob'] = d.prob_uniform_prior.swifter.apply(format_tuple_string)
    else:
        d['prob'] = d.prob_seq_pattern_prior.swifter.apply(format_tuple_string)
    d['best_prob'] = d.prob.apply(max)
    
    # filter based on SIQ and strongest probabiltiy
    d = d[(d.SIQ >= args.SIQ_thres) & (d.best_prob >= args.prob_thres)]
    d['NanoSplicer_junction'] = list(d.swifter.apply(
        lambda x: x.candidates[np.argmax(x.prob)], axis = 1))

    d.set_index(
        ['reference_name', 'read_id',  'initial_junction'], inplace=True)
    
    # select which columns in NanoSplicer output to keep here
    d.drop(columns = d.columns.difference(['is_hcjwr', 'NanoSplicer_junction', 'candidates', 'prob', 'SIQ', 'best_prob']), inplace = True)
    return d
    
def parse_nanosplicer_jwr_h5(args):
    """Read the h5 file from the jwr checker
    """
    all_jwr = pd.read_hdf(args.jwr_check_h5, 'data')
    all_jwr = all_jwr.rename(columns={'id':'read_id', 'loc':'initial_junction', 'chrID':'reference_name'})
    logger.info(f'Indexing jwr file ...')
    all_jwr.set_index(['reference_name', 'read_id',  'initial_junction'], inplace=True)
    all_jwr.drop(columns = all_jwr.columns.difference(['transcript_strand', 'JAQ']))
    return all_jwr

def get_groups(junc_serie, max_diff):
    """Group together the jwr with nearby location. In each group, it is not
    guarenteed that any two JWRs will be within the specified "max_diff", but for
    each JWR, there will be at lease one JWR i is in the range. 

    Args:
        junc_serie pd.Series: list of junction
        max_diff (int, optional): Defaults to 20.

    Returns:
        np.array: An array of group index match the order of input junction
    """
    def __get_groups_from_sorted_array(array, max_diff):
        # note that array has to be sorted
        # padding
        array_pad = np.hstack([array[0], array])
        return np.cumsum(np.abs(array_pad[1:] - array_pad[:-1]) > max_diff)

    def __grouping(array, max_diff):
        # add index
        array_indexed = np.vstack([np.arange(len(array)), array]).T
        # sort 
        array_indexed = array_indexed[array_indexed[:, 1].argsort()]
        # get group
        rst_group = __get_groups_from_sorted_array(array_indexed[:,1], max_diff)
        # revert oder
        return rst_group[array_indexed[:, 0].argsort()]

    junc_serie = junc_serie.apply(lambda x: np.array(x).flatten())
    junc_serie = np.vstack(junc_serie)

    # junc_serie: a columns is a splice junc
    group_by_sites = np.vstack([__grouping(col, max_diff) for col in junc_serie.T]).T
    group = np.unique(group_by_sites, axis=0, return_inverse=True)[1]  
    return group

def grp_generator(d, group_index, group_by=''):    
    """Generator to create copy of dataframe per group

    Args:
        d (pd.DataFrame): original frame to group
        group_index (list), index of group index, length of the list is the rows in d 
        group_by (str, optional): Description for tqdm bar.
    """
    # get group of jwrs in
    d = d.copy()
    d['group_idx'] = group_index
    d_grb = d.groupby('group_idx')
    for _, d_single in  tqdm(d_grb.__iter__(), 
                        total = len(d_grb.groups.keys()),
                        desc=group_by):
        yield d_single.copy()

def hcjwr_identification(d, args):
    """Identify hcjwr

    Args:
        d (DataFrame): DataFrame with combined NanoSplicer and minimap2 result
    return:
        Datafrome of unique HC junction with count
    """
    d = d.copy()
    thres_pass = (d.JAQ >= args.hcjwr_JAQ_thres) & \
                    (d.SIQ >= args.hcjwr_SIQ_thres) & \
                    (d.best_prob >= args.hcjwr_prob_thres)
    
    if args.hcjwr_consistence:
        thres_pass = thres_pass & (d.initial_junction == d.NanoSplicer_junction)
    
    hwjwr = d[thres_pass].initial_junction
    hcjwr_count_df = hwjwr.value_counts().to_frame('counts').reset_index()
    hcjwr_count_df.rename(columns = {'index':'initial_junction'}, inplace=True)
    # apply minimum count
    hcjwr_count_df = hcjwr_count_df[hcjwr_count_df.counts > args.hcjwr_min_count]
    return hcjwr_count_df, thres_pass

def single_jwr_group_correction(d_grp, args):
    """Function for correcting JWRs in a single junction group containing 
    nearby junctions.

    Step 1. Get hc junction list and count
    Step 2. Correct non-hc JWR

    Args:
        d_grp 
    Returns:
        prob table from NanoSplicer with updated columns:
            corrected (bool)
            HC_candidates (tuples)
            HC_prob (tuples)
    """
    # identify hc junction and index of which JWR are hc
    hc_junc, is_hcjwr = hcjwr_identification(d_grp, args)

    # when not hc junction find, mark junction as (-1, -1)
    if not len(hc_junc):
        d_grp['corrected_junction'] = [(-1,-1)] * len(d_grp)
        d_grp['HC_candidates'] = [()] * len(d_grp)
        d_grp['HC_prob'] = [()] * len(d_grp)
        # need to save the result
        return d_grp
    else:
        pass

    d_grp_to_correct = d_grp.copy()
    d_grp_to_correct.reset_index(inplace=True)

    # get HC candidates for each JWR (overlap between NanoSplicer candidate and HC junction)
    uniq_candidates = d_grp_to_correct['candidates'].drop_duplicates().dropna()
    hc_junc_set = set(hc_junc.initial_junction)
    candidate_map_dict = {k:tuple(set(k) & hc_junc_set) for k in uniq_candidates}
    candidate_map_dict[np.nan]=()
    d_grp_to_correct['HC_candidates'] = d_grp_to_correct.candidates.map(candidate_map_dict)
    
    # perform corretion
        # no. of HC junction nearby == 0: mark jwr as (-1, -1)
        # == 1: correct to it
        # > 1: do not correct, mark jwr as np.nan adding the probability of each HC junc for each JWR
    corrected_junctions = []
    hcjwr_candidate_probs = []
    for row in d_grp_to_correct.itertuples():   
        if pd.isna(row.HC_candidates) or len(row.HC_candidates) == 0:
            # completely missed jwr
            corrected_junctions.append((-1,-1))
            hcjwr_candidate_probs.append(())
        elif len(row.HC_candidates) == 1:
            corrected_junctions.append(row.HC_candidates[0])
            hcjwr_candidate_probs.append(())
        else:
            candidate = list(row.candidates)
            hcjwr_candidate_probs.append(
                tuple([row.prob[candidate.index(x)] for x in row.HC_candidates]))
            corrected_junctions.append(np.nan)
    d_grp_to_correct['corrected_junction'] = corrected_junctions
    d_grp_to_correct['hcjwr_candidates_probs'] =  hcjwr_candidate_probs

    return d_grp_to_correct

def group_and_correct_uncorrected_jwr(args, d, max_diff):
    """
    1. group JWRs 
    2. For each group: 
        1. get HCJWR
        2. correct in group

    Note, all junction are assumed to be on a same chromosome in this function

    Args:
        d (pd.DataFrame): prob table from NanoSplicer

    Returns:
        prob table from NanoSplicer with updated columns:
            corrected (bool)
            HC_candidates (tuples)
            HC_prob (tuples)
    """
    try:
        assert len(d.reference_name.unique()) == 1
    except AssertionError as msg:
        print(msg)
        sys.exit(1)

    group_idx = get_groups(d.initial_junction, max_diff)

    logger.info(f'Correcting JWRs based on High-confidence JWRs ...')
    # get group of jwrs in
    d=d.copy()
    d['group_idx'] = group_idx
    tqdm.pandas(desc="Processing nearby JWR groups", position=10, leave = False)
    d_correct = d.groupby('group_idx').progress_apply(
        lambda x: single_jwr_group_correction_v2(x, args))
    return d_correct

def single_jwr_group_correction_v2(d_grp, args):
    """Function for correcting JWRs in a single junction group containing nearby junctions.

    Args:
        d_grp (pd.DataFrame)
        args

    Step 1. Get hc junction list and count
    Step 2. Correct non-hc JWR (keep JWR uncorrected when there are multiple hc junction)

    Returns:
        prob table from NanoSplicer with updated columns:
            corrected (bool)
            HC_candidates (tuples)
            HC_prob (tuples)
    """
    if np.all(~d_grp.corrected_junction.isna()):
        return d_grp

    hc_junc, is_hcjwr = hcjwr_identification(d_grp, args)

    if not len(hc_junc):
        d_grp['corrected_junction'] = [(-1,-1)] * len(d_grp)
        d_grp['HC_candidates'] = [()] * len(d_grp)
        d_grp['HC_prob'] = [()] * len(d_grp)
        # need to save the result
        return d_grp

    # d_grp_corrected = d_grp[is_hcjwr].copy()
    # d_grp_corrected['corrected_junction'] = d_grp_corrected.NanoSplicer_junction
    # d_grp_corrected['HC_candidates'] = [()] * len(d_grp_corrected)
    # d_grp_corrected['HC_prob'] = [()] * len(d_grp_corrected)
    
    # correcting JWRs that are not HC
    # temp change: org
    # d_grp_to_correct = d_grp[~is_hcjwr].copy()

    # temp
    d_grp_to_correct = d_grp.copy()
    d_grp_to_correct.reset_index(inplace=True)

    # get overlap between NanoSplicer candidate and HC junction
    uniq_candidates = d_grp_to_correct['candidates'].drop_duplicates().dropna()
    hc_junc_set = set(hc_junc.initial_junction)
    candidate_map_dict = {k:tuple(set(k) & hc_junc_set) for k in uniq_candidates}
    candidate_map_dict[np.nan]=()
    d_grp_to_correct['HC_candidates'] = d_grp_to_correct.candidates.map(candidate_map_dict)
    
    corrected_junctions = []
    hcjwr_candidate_probs = []
    for row in d_grp_to_correct.itertuples():   

        # perform corretion
        if pd.isna(row.HC_candidates) or len(row.HC_candidates) == 0:
            # completely missed jwr
            corrected_junctions.append((-1,-1))
            hcjwr_candidate_probs.append(())

        elif len(row.HC_candidates) == 1:
            corrected_junctions.append(row.HC_candidates[0])
            hcjwr_candidate_probs.append(())

        else:
            candidate = list(row.candidates)
            hcjwr_prob = tuple([row.prob[candidate.index(x)] for x in row.HC_candidates])
            hcjwr_candidate_probs.append(hcjwr_prob)
            if row.SIQ > args.SIQ_thres and max(hcjwr_prob) > args.prob_thres:
                corrected_junctions.append(row.HC_candidates[np.argmax(hcjwr_prob)]) 
            elif row.initial_junction in row.HC_candidates \
                and row.JAQ > args.JAQ_thres:
                corrected_junctions.append(row.initial_junction) 
            else:
                corrected_junctions.append(np.nan) 
            #corrected_junctions.append(np.nan)
    d_grp_to_correct['corrected_junction'] = corrected_junctions
    d_grp_to_correct['hcjwr_candidates_probs'] =  hcjwr_candidate_probs

        
    #return pd.concat([d_grp_to_correct,d_grp_corrected])
    return d_grp_to_correct

def group_and_correct(args, d, max_diff):
    """
    1. group JWRs based on max_diff
    2. For each group: 
        1. get HCJWR
        2. correct in group

    Note, all JWR in input dataframe (d) are assumed to be on a same chromosome in this function

    Args:
        d (pd.DataFrame): prob table from NanoSplicer

    Returns:
        prob table from NanoSplicer with updated columns:
            corrected (bool)
            HC_candidates (tuples)
            HC_prob (tuples)
    """
    # check assumption
    try:
        assert len(d.reference_name.unique()) == 1
    except AssertionError as msg:
        print(msg)
        sys.exit(1)

    # Assign group id in the order of input JWRs. Nearby JWRs will have same group id.
    group_idx = get_groups(d.initial_junction, max_diff)
    d['group_idx'] = group_idx

    # correct each single group
    tqdm.pandas(desc="Processing nearby JWR groups")
    d_correct = d.groupby('group_idx').progress_apply(
        lambda x: single_jwr_group_correction(x, args))
    
    return d_correct


def restructure_per_jwr_dataframe(all_jwr, tss_tts_d, 
                                  on='corrected_junction',
                                  output_uncorrected = True,
                                  rm_cpl_missed_from_ends = True,
                                  summary = True):
    """Restructure a DataFrame from JWR per row to reads per row.

    Args:
        all_jwr (df): DataFrame with JWR per row
        tss_tts_d (df): DataFrame containing transcript start and end site for each read
        on (str, optional): Which junction to use. Defaults to 'corrected_junction'.
        output_uncorrected (bool, optional): 
            Whether output a DataFrame with JWRs from uncorrected reads. Defaults to True.
        rm_cpl_missed_from_ends (bool, optional): 
            Whether or not to remove the completely missed junctions if they are
            the first/last junctions. Defaults to True.

    Returns:
        Output 1:
            Restructured table:
            Columns:
            'reference_name':reference/chromosom name 
            'read_id':read_id 
            'transcript_strand':transcript_strand
            'junc': list of corrected junctions
            'junc_count' number of junctions
            'junc_start': 5' splice site of each junction
            'junc_end': 3' splice site of each junction
        Note: if on='corrected_junction', this will only include reads with all JWR corrected.
        Output 2 (if `output_uncorrected=True`): 
            DataFrame containing more infomation of reads with uncorrected JWRs (This is a JWR table)
    """
    ref_names, read_ids, transcript_strands, num_junc, juncs,\
    junc_starts, junc_ends, intron_removed_ends, intron_removed_starts =\
            [],[],[],[],[],[],[],[],[]
    
    uncorrected_read = []
    uncorrected_read_count = 0
    discard_read_count = 0
    all_jwr_gb = all_jwr.groupby(['reference_name', 
                                    'read_id'])
    
    for key, df in tqdm(all_jwr_gb.__iter__(), 
                            total = len(all_jwr_gb.groups.keys()),
                            desc='Restructuring reads'):
        df.sort_values('initial_junction', inplace=True)
        df.reset_index(inplace=True, drop=True)
        # temp
        if any(df[on].isnull()):
            uncorrected_read.append(df)
            uncorrected_read_count += 1
            continue

        is_hc_junc = df['corrected_junction'] != (-1,-1)
        
        # remove completely missed junction in the middle
        if sum([k for k, g in itertools.groupby(is_hc_junc)])!=1:
            discard_read_count += 1
            continue

        
        ref_names.append(list(df.reference_name)[0])
        read_ids.append(df.read_id[0])
        transcript_strands.append(df.transcript_strand[0])
        
        intron_removed_start = 0
        intron_removed_end = 0
        if rm_cpl_missed_from_ends:
            # remove completely missed junction at begining or end
            rst_junc = tuple(df[on][is_hc_junc])
            for junc, i in zip(df['initial_junction'], is_hc_junc):
                if i == 1:
                    break
                intron_removed_start += np.abs(junc[1]-junc[0])

            for junc, i in zip(df['initial_junction'][::-1], is_hc_junc[::-1]):
                if i:
                    break
                intron_removed_end += np.abs(junc[1]-junc[0])
            
        else:
            rst_junc = tuple(df[on])
            
        
        juncs.append(tuple(rst_junc))
        junc_starts.append(tuple([x for x,y in rst_junc]))
        junc_ends.append(tuple([y for x,y in rst_junc]))
        intron_removed_ends.append(intron_removed_end)
        intron_removed_starts.append(intron_removed_start)
    all_read = pd.DataFrame({'reference_name':ref_names,  
                            'read_id':read_ids,  
                            'transcript_strand':transcript_strands, 
                            'junc': juncs,
                            'junc_start':junc_starts, 
                            'junc_end':junc_ends,
                            'intron_removed_start': intron_removed_starts,
                            'intron_removed_end':intron_removed_ends
    })

    all_read = all_read.merge(tss_tts_d, on = 'read_id')
    all_read['junc_count'] = all_read.junc.apply(len)

    
    if output_uncorrected:
        logger.info('Concating table...')
        uncorrected_read = pd.concat(uncorrected_read) if len(uncorrected_read) else None
        return all_read, uncorrected_read
    else:
        return all_read, None

# main function for correction
def correction_round1(args, tss_tts_d):
    # read input NanoSplicer result
    logger.info(f'Parsing jwr checker file ...')
    all_jwr = parse_nanosplicer_jwr_h5(args)
    logger.info(f'Parsing prob table file ...')
    prob_table = parse_nanosplicer_prob_table(args)
    
    # correct JWRs in all_jwrs using NanoSplicer output
    logger.info(f'Preparing data frame for correction...')
    all_jwr = all_jwr.merge(prob_table, 
                            how='left', 
                            #right_on = ['SIQ', 'best_prob'],
                            left_index=True,
                            right_index=True)
    logger.info(
        helper.green_msg(
            f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}',
            print_out=False))

    # proceed each chromosome strand
    logger.info(f'Correcting JWRs based on High-confidence JWRs ... (Round 1)')
    all_jwr.reset_index(inplace = True)
    all_jwr_chr_grps = all_jwr.groupby(by=['reference_name', 'transcript_strand'])
    corrected_all_jwr = []
    for key, all_jwr_single_chr in tqdm(all_jwr_chr_grps.__iter__(), 
                            total = len(all_jwr_chr_grps.groups.keys()),
                            desc='Chromosome (strand-aware)', position=0, leave = False):
        # ignore result without strand information
        if key[1] not in ['+', '-']:
            continue
        corrected_all_jwr.append(
            group_and_correct(args, all_jwr_single_chr, 
                              max_diff=CORRECTION_ARG['dist']))

    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')
    
    logger.info(f'Converting JWR table to read table...')
    corrected_all_jwr = pd.concat(corrected_all_jwr)
    
    return corrected_all_jwr
    #all_read, uncorrected_jwr = restructure_per_jwr_dataframe(corrected_all_jwr, tss_tts_d)
    # return all_read, uncorrected_jwr

    # add_summary(textwrap.dedent(
    #     f'''
    #     After correcting based on nearby JWRs:
    #         Number of reads with all JWRs corrected: {len(all_read)}
    #         Number of reads containing uncorrected JWRs: {uncorrected_read_count}
    #         Number of reads discarded: {discard_read_count}
    #     '''))
    