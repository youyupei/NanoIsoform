#!/usr/bin/env python3 

# goal correct JWR without NanoSplicer output
import argparse
import importlib
import warnings
import textwrap
import pandas as pd
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm
import numpy as np
import logging

from config import *
import helper
from helper import add_summary

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
    parser.add_argument('--cfg_fn', type=str, default='',
                        help='Filename of customised config file.')

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
    # update config if provided
    if args.cfg_fn:
        #mdl = importlib.import_module(args.cfg_fn)
        spec = importlib.util.spec_from_file_location(
                os.path.basename(args.cfg_fn).split('.')[0], args.cfg_fn)
        mdl = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mdl)
        if "__all__" in mdl.__dict__:
            names = mdl.__dict__["__all__"]
        else:
            names = [x for x in mdl.__dict__ if not x.startswith("_")]
        globals().update({k: getattr(mdl, k) for k in names})
        args = parser.parse_args()

    # check file 
    helper.check_exist([args.prob_table, args.jwr_check_h5])
    return args

args = parse_arg()
from __restructure_input import *
from __group_and_correct import *
# ACTION REQUIRED FOR THE FINAL VERSION
#   remove levelname, filname, funcName, lineno for the final version
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# set up pandas
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', 50)


# output 
def output_h5_file_corrected_reads(d, filename, key='data', csv_output=OUTPUT_CSV):
    """Output h5 file containing the reads with all JWRs corrected. Those without
    all JWRs corrected will not be present in this output file.
    """
    output_d = d[d.all_corrected == True].copy()
    output_d['junc_start'] = output_d.corrected_junction.apply(lambda y: tuple([x[0] for x in y]))
    output_d['junc_end'] = output_d.corrected_junction.apply(lambda y: tuple([x[1] for x in y]))
    output_d.to_hdf(filename, key)
    if csv_output:
        output_d.to_csv(filename+'.csv')
    logger.info(helper.green_msg(f"Output saved as {args.output_fn}!"))


def main(args):
    logger.info("Formatting input data...")
    # get input data and reformat
    all_reads = parse_format_input_file(args)
    logger.info(helper.mem_time_msg())
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
    logger.info(helper.mem_time_msg())
    if args.no_nearby_jwr_correction:
        all_reads[all_reads.corrected.apply(all)].to_hdf(args.output_fn, key='data')
        all_reads[all_reads.corrected.apply(all)].to_csv(args.output_fn+'.csv')
        return None
    else:
        pass

    logger.info('Correcting reads in each group...')
    # get corrected junction for groups
    corrected_d = correct_junction_per_group(
                    all_reads, methods=args.nearby_jwr_correction_mode)
    logger.info(helper.mem_time_msg())
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
    logger.info(helper.mem_time_msg())
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
    logger.info(helper.mem_time_msg())
    logger.info('Correcting reads in each group...')
    # get corrected junction for groups
    corrected_d = correct_junction_per_group(
                    all_reads, methods=args.nearby_jwr_correction_mode)
    
    logger.info(helper.mem_time_msg())
    output_h5_file_corrected_reads(corrected_d, args.output_fn, key='data')

    # add some text summary
    add_summary(textwrap.dedent(
        f'''After correcting based on nearby JWRs:
            Number of reads with all JWRs corrected: {np.sum(corrected_d.all_corrected)}
        '''))

if __name__ == '__main__':
    print(args)
    if args.test_mode:
        test(args)
    else:
        main(args)

    print('\n\nSummary:\n', helper.summary_msg)