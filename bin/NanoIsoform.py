#!/usr/bin/env python3 

# goal correct JWR without NanoSplicer output
import argparse
import textwrap
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', 50)
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm
import numpy as np
import logging


import helper
from config import *

# ACTION REQUIRED FOR THE FINAL VERSION
#   remove levelname, filname, funcName, lineno for the final version
logging.basicConfig(format='[%(asctime)s] %(levelname)s [%(filename)s.%(funcName)s:%(lineno)d] %(message)s', datefmt='%a, %d %b %Y %H:%M:%S')
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.info('Start running...')

# store some useful stats during the run
summary_msg = ''
def add_summary(msg):
    global summary_msg
    summary_msg += msg + '\n'
    
def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        '''),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Required positional argument
    parser.add_argument('prob_table', type=str,
                        help='Filename of the probability table output from NanoSplicer')
    parser.add_argument('jwr_check_h5', type=str,
                        help='Filename of the HDF5 file output from NanoSplicer module'
                        'jwr_checker')

    # Optional argument
    parser.add_argument('--SIQ_thres', type=float, default=DEFAULT_INPUT['SIQ_thres'],
                        help= textwrap.dedent(
                            '''
                            SIQ threshold for high qualilty squiggle matching （JWRs with）
                            SIQ < threshold will be ignored from NanoSplicer output.
                            '''))
    parser.add_argument('--prob_thres', type=float, default=DEFAULT_INPUT['prob_thres'],
                        help= textwrap.dedent(
                            '''
                            The minimum probability of the NanoSplicer identified junction.
                            NanoSplicer identified junction with probability < threshold
                            will be ignored
                            '''))
    parser.add_argument('--JAQ_thres', type=float, default=DEFAULT_INPUT['JAQ_thres'],
                        help= textwrap.dedent(
                            '''
                            Fow JWRs that NanoSplicer does not give identification to. The initial mapped location 
                            will be used as "corrected junction" if the Junction Alignment Quality (JAQ) 
                            of the initial mapping >= this threshold. 
                            '''))
    
    parser.add_argument('--nearby_jwr_correction_mode', type=str, choices=['majority_vote', 'probability'], default=CORRECTION_ARG['method'],
                        help= textwrap.dedent(
                            '''
                            How to correct jwr based on nearby corrected jwrs.
                            'majority_vote': Corrected to most supported junction nearby
                            'probability': randomly choose from the nearby junctions with 
                                            probability based on the proportion of support. 
                            '''))

    parser.add_argument('--uniform_prior', action='store_true',
                        help='Output the whitelist in Cellranger style')
    parser.add_argument('--output_fn', type=str, default=DEFAULT_INPUT['output_fn'],
                        help='Output filename')

    # Developer only argument
    Dev_arg = parser.add_argument_group('For the developers only:')
    Dev_arg.add_argument('--test_mode', action='store_true',
                        help='Run in test mode.')
    Dev_arg.add_argument('--no_nearby_jwr_correction', action = 'store_true',
                        help= textwrap.dedent(
                            '''For test purpose only: Turn off the correction based on nearby jwr.
                            Note that reads with uncorrected jwr(s) will not appear in the 
                            output.
                            '''))
    args = parser.parse_args()

    # check file 
    helper.check_exist([args.prob_table, args.jwr_check_h5])
    return args

def parse_format_input_file(args):
    """Perform the following steps:
        1. Read two input files as pd.DataFrame
        2. Filter NanoSplicer prob table using SIQ and strongest assignment prob
        3. merge two filtered DataFrame (JWRs corrected to NanoSplicer identified ones)
        4. Further Correct JWR based on JAQ
        5. reformat the Dataframe (a jwr per row -> a read per row)
    """

    def parse_nanosplicer_prob_table(args):
        """Parse the prob table file output by NanoSplicer as 

        Args:
            tb_fn (str): NanoSplicer file name
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
        d['initial_junction'] = d.initial_junction.apply(format_initial_junction)
        d['candidates'] = d.candidates.apply(format_candidates)
        #d['candidate_sequence_motif_preference'] = d.candidate_sequence_motif_preference.apply(format_tuple_string)
        
        
        
        if args.uniform_prior:
            d['prob_uniform_prior'] = d.prob_uniform_prior.apply(format_tuple_string)
            d['best_prob'] = d.prob_uniform_prior.apply(max)
            d = d[(d.SIQ >= args.SIQ_thres) & (d.best_prob >= args.prob_thres)]
            d['corrected_junction'] = list(d.apply(
                lambda x: x.candidates[np.argmax(x.prob_uniform_prior)], axis = 1))
        else:
            d['prob_seq_pattern_prior'] = d.prob_seq_pattern_prior.apply(format_tuple_string)
            d = d[(d.SIQ >= args.SIQ_thres) & (d.best_prob >= args.prob_thres)]       
            d['corrected_junction'] = list(d.apply(
                lambda x: x.candidates[np.argmax(x.prob_seq_pattern_prior)], axis = 1))
            
        logger.info(f'Indexing prob table file ...')
        d.set_index(
            ['reference_name', 'read_id',  'initial_junction'], inplace=True)
        
        d.drop(columns = d.columns.difference(['corrected_junction']), inplace = True)
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

    def restructure_per_jwr_dataframe(all_jwr):
        """Restructured table:
        Columns:
        'reference_name':reference/chromosom name 
        'read_id':read_id 
        'transcript_strand':transcript_strand
        'junc_count': Number of junctions
        'junc_start': 5' splice site of each junction
        'junc_end': 3' splice site of each junction
        'corrected': whether the corresponding junction is corrected
        'JAQ': Junction alignment quality

        Note: In the junc start and junce end, the mapped splice sites were recorded
        if they have not been corrected.
        """
        ref_names,  read_ids,  transcript_strands, num_junc, juncs,\
        junc_starts, junc_ends, processed, JAQs= [],[],[],[],[],[],[],[],[]
        for key, df in all_jwr.groupby(level=1):
            ref_names.append(df.index.get_level_values(0)[0])
            read_ids.append(df.index.get_level_values(1)[0])
            JAQs.append(tuple(df.JAQ))
            num_junc.append(len(df.JAQ))
            transcript_strands.append(df.transcript_strand[0])
            
            identified = tuple(~df.corrected_junction.isnull())
            processed.append(identified)
            mapping = df.index.get_level_values('initial_junction') 
            corrected = df.corrected_junction 
            df_junc = \
                [corrected[i] if identified[i] else mapping[i] for i in range(len(corrected))]
            df_junc = sorted(df_junc)
            juncs.append(tuple(df_junc))
            junc_starts.append(tuple([x for x,y in df_junc]))
            junc_ends.append(tuple([y for x,y in df_junc]))

        all_read = pd.DataFrame({'reference_name':ref_names,  
                                'read_id':read_ids,  
                                'transcript_strand':transcript_strands, 
                                'junc': juncs,
                                'junc_count': num_junc,
                                'junc_start':junc_starts, 
                                'junc_end':junc_ends, 
                                'corrected': processed, 
                                'JAQ':JAQs
        })

        return all_read

    # read the input file (prob_table, all_jwr.h5)
    logger.info(f'Parsing prob table file ...')
    prob_table = parse_nanosplicer_prob_table(args)
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    logger.info(f'Parsing jwr file ...')
    all_jwr = parse_nanosplicer_jwr_h5(args)
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    # correct JWRs in all_jwrs using NanoSplicer output
    logger.info(f'Correcting JWRs using NanoSplicer output ...')
    all_jwr = all_jwr.merge(
        prob_table, 
        how='left', 
        # right_on = ['corrected_junction'],
        left_index=True, 
        right_index=True)
    del prob_table
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    # correct JWRs with JAQ > args.JAQ_thres
    logger.info(f'Correcting JWRs using JAQ threshold ...')
    JAQ_pass = (all_jwr.corrected_junction.isnull()) & (all_jwr.JAQ >= args.JAQ_thres)
    all_jwr.loc[JAQ_pass, 'corrected_junction'] =\
        all_jwr.loc[JAQ_pass,].index.get_level_values('initial_junction')  
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    # add some text summary
    add_summary(textwrap.dedent(
        f'''
        Total number of JWRs: {len(all_jwr)}
            Number of JWRs without NanoSplicer identification: {np.sum(all_jwr.corrected_junction.isnull())}
            Number of JWRs uncorrected after JAQ-based correction: {np.sum((all_jwr.corrected_junction.isnull()) & (all_jwr.JAQ < args.JAQ_thres))}
        '''
    ))
    
    # restructure table jwr per row -> read per raw
    logger.info('Restructuring table...')
    all_read = restructure_per_jwr_dataframe(all_jwr)
    del all_jwr
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    return all_read

def group_reads(all_reads, max_diff = GROUP_ARG['max_diff']):
    """group reads based on their junctions
    Args:
        all_reads (pandas.DataFrame): The dataframe output by restructure_per_jwr_datafrome
        function.
        max_diff (INT): maximum different allowed when grouping closeby junction

    Output:
        read_group (dataframe or list): indicate the membership of each group
    """
    def split_groupby(df, max_diff):
        org_d = df.copy()
        temp_d = df.index.unique().to_frame().copy()
        junc_array = temp_d.junc.apply(lambda x: np.array(x).flatten())
        junc_array = np.vstack(junc_array)
        test = np.vstack([get_groups(col, max_diff) for col in junc_array.T]).T
        group = np.unique(test, axis=0, return_inverse=True)[1]
        temp_d['sub_group'] = group
        merge_d = org_d.merge(temp_d, 'left', left_index=True, right_index=True)
        merge_d = merge_d.drop(columns = ['reference_name', 'junc_count', 'junc'])
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
        
    # the grouping 
    all_reads = all_reads.set_index(['reference_name',
                            'junc_count', 'junc'])   

    # all reads group by chr/number of junctions/ junctions()
    key, test_df = list(all_reads.groupby(level=[0,1]).__iter__())[1]
    #all_reads = all_reads.groupby(level=[0,1]).apply(split_groupby, max_diff)
    df_list = []
    for k, d in all_reads.groupby(level=[0,1]).__iter__():
        df_list.append(split_groupby(d, max_diff))
    del all_reads
    
    all_reads = pd.concat(df_list)

    if all_reads.index.nlevels > 3:
        # sometime, after merging, there are duplicated columns. Not sure about
        # the reason
        all_reads = all_reads.droplevel([0,1])

    return all_reads

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
                        min_prop=CORRECTION_ARG['maj_vot_min_prop'], 
                        min_correct_read=CORRECTION_ARG['maj_vot_min_count'])
                    
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
    if prop.max() >= min_prop and values.max() >= min_correct_read:
        return keys[np.random.choice(np.array(range(len(keys)))[prop == prop.max()])]
    else:
        return None



def output_h5_file_corrected_reads(d, filename, key='data'):
    output_d = d[d.all_corrected == True].copy()
    output_d['junc_start'] = output_d.corrected_junction.apply(lambda y: tuple([x[0] for x in y]))
    output_d['junc_end'] = output_d.corrected_junction.apply(lambda y: tuple([x[1] for x in y]))
    output_d.to_hdf(filename, key)
    output_d.to_csv(filename+'.csv')
    logger.info(helper.green_msg(f"Output saved as {args.output_fn}!"))

def main(args):
    print(args)
    logger.info("Formatting input data...")
    # get input data and reformat
    all_reads = parse_format_input_file(args)
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    # add some text summary
    add_summary(textwrap.dedent(
        f'''
        After correcting using NanoSplicer and JAQ >= {args.JAQ_thres}:
            Total number of reads: {len(all_reads)}
            Number of reads with all JWRs corrected: {np.sum(all_reads.corrected.apply(all))}
        '''))


    logger.info('Grouping reads...')
    all_reads = group_reads(all_reads,max_diff = GROUP_ARG['max_diff'])
    all_reads.reset_index(drop=False, inplace = True)
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    if args.no_nearby_jwr_correction:
        all_reads[all_reads.corrected.apply(all)].to_hdf(args.output_fn, key='data')
        all_reads[all_reads.corrected.apply(all)].to_csv(args.output_fn+'.csv')
        return None

    logger.info('Correcting reads in each group...')
    # get corrected junction for groups
    corrected_d = correct_junction_per_group(
                    all_reads, methods=args.nearby_jwr_correction_mode)
    print(corrected_d)
    print(corrected_d.junc.apply(len))
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    output_h5_file_corrected_reads(corrected_d, args.output_fn, key='data')
    
    # add some text summary
    add_summary(textwrap.dedent(
        f'''
        After correcting based on nearby JWRs:
            Number of reads with all JWRs corrected: {np.sum(corrected_d.all_corrected)}
        '''))

def test(args):
    # test setup
    cached_data = \
        '/home/ubuntu/data/github_repo/youyupei/NanoIsoform/test/large_set.h5'
    logger.info('Reading input dataset...')
    helper.check_memory_usage()
    all_reads = pd.read_hdf(cached_data, 'data')
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    # add some text summary
    add_summary(textwrap.dedent(
        f'''
        After correcting using NanoSplicer and JAQ >= {args.JAQ_thres}:
            Total number of reads: {len(all_reads)}
            Number of reads with all JWRs corrected: {np.sum(all_reads.corrected.apply(all))}
        '''))

    if args.no_nearby_jwr_correction:
        all_reads[all_reads.corrected.apply(all)].to_hdf(args.output_fn, key='data')
        all_reads[all_reads.corrected.apply(all)].to_csv(args.output_fn+'.csv')
        return None
        
    logger.info('Grouping reads...')
    all_reads = group_reads(all_reads,max_diff = GROUP_ARG['max_diff'])
    all_reads.reset_index(drop=False, inplace = True)
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    logger.info('Correcting reads in each group...')
    # get corrected junction for groups
    corrected_d = correct_junction_per_group(
                    all_reads, methods=args.nearby_jwr_correction_mode)
    
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')
    

    output_h5_file_corrected_reads(corrected_d, args.output_fn, key='data')

    # add some text summary
    add_summary(textwrap.dedent(
        f'''After correcting based on nearby JWRs:
            Number of reads with all JWRs corrected: {np.sum(corrected_d.all_corrected)}
        '''))


if __name__ == '__main__':
    args = parse_arg()
    if args.test_mode:
        # test command line input
        print(args)
        test(args)
    else:
        main(args)
    
    print('Summary:\n',summary_msg)