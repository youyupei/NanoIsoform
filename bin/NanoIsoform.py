#!/usr/bin/env python3 

# goal correct JWR without NanoSplicer output
import argparse
import textwrap
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', 20)
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm
import numpy as np


import helper
from config import *

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

    # required name argment
    # requiredNamed = parser.add_argument_group('Either one of these argument is required')
    # requiredNamed.add_argument('--expect-cells',type=int, help='<INT>:  Expected number of cells.')
    # requiredNamed.add_argument('--count-threshold', type=int,
    #                     help='Output the whitelist in Cellranger style')

    # Optional positional argument
    #parser.add_argument('opt_pos_arg', type=int, nargs='?',help)

    # Optional argument
    parser.add_argument('--SIQ_thres', type=float, default='-0.8',
                        help= textwrap.dedent(
                            '''
                            SIQ threshold for high qualilty squiggle matching （JWRs with）
                            SIQ < threshold will be ignored from NanoSplicer output.
                            '''))
    parser.add_argument('--prob_thres', type=float, default='0.8',
                        help= textwrap.dedent(
                            '''
                            The minimum probability of the NanoSplicer identified junction.
                            NanoSplicer identified junction with probability < threshold
                            will be ignored
                            '''))
    parser.add_argument('--JAQ_thres', type=float, default='0.95',
                        help= textwrap.dedent(
                            '''
                            Fow JWRs that NanoSplicer does not give identification to. The initial mapped location 
                            will be used as "corrected junction" if the Junction Alignment Quality (JAQ) 
                            of the initial mapping >= this threshold. 
                            '''))
    parser.add_argument('--uniform_prior', action='store_true',
                        help='Output the whitelist in Cellranger style')
    


    args = parser.parse_args()

    # check file 
    helper.check_exist([args.prob_table, args.jwr_check_h5])
    return args



def parse_format_input_file(args):

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
            d['corrected_junction'] = d.apply(
                lambda x: x.candidates[np.argmax(x.prob_uniform_prior)], axis = 1)
        else:

            d['prob_seq_pattern_prior'] = d.prob_seq_pattern_prior.apply(format_tuple_string)
            d = d[(d.SIQ >= args.SIQ_thres) & (d.best_prob >= args.prob_thres)]
            d['corrected_junction'] = d.apply(
                lambda x: x.candidates[np.argmax(x.prob_seq_pattern_prior)], axis = 1)
        print(f'Indexing prob table file ...')
        d.set_index(
            ['reference_name', 'read_id',  'initial_junction'], inplace=True)
        return d[['corrected_junction']]

    def parse_nanosplicer_jwr_h5(args):
        """Read the h5 file from the jwr checker
        """
        all_jwr = pd.read_hdf(args.jwr_check_h5, 'data')
        all_jwr = all_jwr.rename(columns={'id':'read_id', 'loc':'initial_junction', 'chrID':'reference_name'})
        print(f'Indexing jwr file ...')
        all_jwr.set_index(['reference_name', 'read_id',  'initial_junction'], inplace=True)
        return all_jwr[['transcript_strand', 'JAQ']]


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
    print(f'Parsing prob table file ...')
    prob_table = parse_nanosplicer_prob_table(args)
    print(f'Parsing jwr file ...')
    all_jwr = parse_nanosplicer_jwr_h5(args)
    
    # correct JWRs in all_jwrs using NanoSplicer output
    print(f'Correcting JWRs using NanoSplicer output ...')
    all_jwr = all_jwr.merge(
        prob_table, 
        how='left', 
        # right_on = ['corrected_junction'],
        left_index=True, 
        right_index=True)
    del prob_table

    # correct JWRs with JAQ > args.JAQ_thres
    print(f'Correcting JWRs using JAQ threshold ...')
    JAQ_pass = (all_jwr.corrected_junction.isnull()) & (all_jwr.JAQ >= args.JAQ_thres)
    all_jwr.loc[JAQ_pass, 'corrected_junction'] =\
        all_jwr.loc[JAQ_pass,].index.get_level_values('initial_junction')  

    # restructure table jwr per row -> read per raw
    all_read = restructure_per_jwr_dataframe(all_jwr)
    del all_jwr

    return all_read

def group_reads(all_reads, max_diff = MAX_DIFF):
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
        return merge_d.drop(columns = ['reference_name', 'junc_count', 'junc'])
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
    all_reads = all_reads.groupby(level=[0,1]).apply(split_groupby, max_diff)     
    return all_reads

def main(args):
    print("parsing skipped for testing purpose..")
    # # get input data and reformat
    # helper.check_memory_usage()
    # all_reads = parse_format_input_file(args)
    # all_reads.to_hdf("large_set.h5", 'data')

    # group reads 
    helper.check_memory_usage()
    print("parsing large_set.h5")
    all_reads = pd.read_hdf("large_set.h5", 'data')
    all_reads = group_reads(all_reads, max_diff=MAX_DIFF).sort_index()
    print(all_reads)

    # filter 

if __name__ == '__main__':
    args = parse_arg()
    # test command line input
    print(args)
    main(args)
    helper.check_memory_usage()
