#!/usr/bin/env python3 

# goal correct JWR without NanoSplicer output




import argparse
import textwrap
import pandas as pd
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm

import helper

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
    parser.add_argument('--uniform_prior', action='store_true',
                        help='Output the whitelist in Cellranger style')
    
    args = parser.parse_args()

    # check file 
    helper.check_exist([args.prob_table, args.jwr_check_h5])
    return args

def parse_nanosplicer_prob_table(args):
    """Parse the prob table file output by NanoSplicer as 

    Args:
        tb_fn (str): NanoSplicer file name
    return:
        pandas data from
    """
    def format_inital_junction(s):
        x1,x2 = s.strip('[()]').split(',')
        return (int(x1), int(x2))

    def format_candidates(s):
        s = [x.strip('[()]').split(',') for x in s.split('),(')]
        return  tuple([(int(x1), int(x2)) for x1, x2 in s])
    def format_tuple_string(s):
        return tuple([float(x) for x in s.split(',')])
    
    d = pd.read_csv(args.prob_table, sep='\t')
    d.inital_junction = d.inital_junction.apply(format_inital_junction)
    d.candidates = d.candidates.apply(format_candidates)
    d.candidate_sequence_motif_preference = d.candidate_sequence_motif_preference.apply(format_tuple_string)
    d.prob_uniform_prior = d.prob_uniform_prior.apply(format_tuple_string)
    d.prob_seq_pattern_prior = d.prob_seq_pattern_prior.apply(format_tuple_string)
    
    d = d[(d.SIQ >= args.SIQ_thres) & (d.best_prob >= args.prob_thres)]
    if args.uniform_prior:
        d.best_prob = d.prob_uniform_prior.apply(max)
    
    return d

def main(args):
    # read the input file (prob_table, all_jwr.h5)
    prob_table = parse_nanosplicer_prob_table(args)
    print(prob_table.columns)
    all_jwr = pd.read_hdf(args.jwr_check_h5, 'data')
    all_jwr = all_jwr.rename(columns={'id':'read_id', 'loc':'inital_junction', 'chrID':'reference_name'})
    print(all_jwr)
    #pd.merge(all_jwr, prob_table, by=['id',])
    # correct JWRs in all_jwrs 
        # if a jwr present in the prob table with good SIQ and strongest prob
        # correct the junction in all_jwr table and mark it as "corrected"
        # elif the JAQ > T
        # leave it to initial mapping but marked as "corrected"

    # restructure table jwr per row -> read per raw

    # group reads 
    pass


    # filter 

if __name__ == '__main__':
    args = parse_arg()
    # test command line input
    print(args)
    main(args)
