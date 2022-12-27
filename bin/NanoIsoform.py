#!/usr/bin/env python3 

# goal correct JWR without NanoSplicer output
import argparse
import pysam
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
    parser.add_argument('input_BAM', type=str,
                        help='Filename of the BAM file output from NanoSplicer module'
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
                        help='Use uniform prior probability for splice patterns. Recommended for synthetic RNA data (e.g. SIRV, Sequins)')
    parser.add_argument('--output_fn', type=str, default=DEFAULT_INPUT['output_fn'],
                        help='Output filename')
    parser.add_argument('--cfg_fn', type=str, default='',
                        help='Filename of customised config file.')
                

    # Developer only argument
    hcjwr_arg = parser.add_argument_group(textwrap.dedent(
        '''
        Definition of High-confidence JWR:
        Note that more stringent threshold than above should be specified. It will be overwritten when less stringet.
        '''))
    hcjwr_arg.add_argument('--hcjwr_SIQ_thres', type=float, default=HCJWR['SIQ_thres'],
                        help= textwrap.dedent(
                            '''
                            SIQ threshold for high qualilty squiggle matching 
                            '''))
    hcjwr_arg.add_argument('--hcjwr_prob_thres', type=float, default=HCJWR['prob_thres'],
                        help= textwrap.dedent(
                            '''
                            The minimum probability of the NanoSplicer identified junction.
                            '''))
    hcjwr_arg.add_argument('--hcjwr_JAQ_thres', type=float, default=HCJWR['JAQ_thres'],
                        help= textwrap.dedent(
                            '''
                            Junction Alignment Quality (JAQ) of the initial mapping >= this threshold. 
                            '''))
    hcjwr_arg.add_argument('--hcjwr_consistence', type=bool, default=HCJWR['consistence'],
                        help= textwrap.dedent(
                            '''
                            NanoSplicer junction is consistent with the minimap2 junction.
                            '''))
    hcjwr_arg.add_argument('--hcjwr_min_count', type=bool, default=HCJWR['min_count'],
                        help= textwrap.dedent(
                            '''
                            Minimum support from HCJWRs for each unique HC junction.
                            '''))


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
import jwr_correction

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


def get_mapped_start_end(fn):
    """finding the mapped start and end sites of each read. 
    Args:
        fn (str): BAM file name
        out_fn (str): output filename (CSV)
    """
    bam = pysam.AlignmentFile(fn, "rb")
    r_id = []
    r_start = []
    r_end = []

    for read in bam.fetch():
        # f_out.write(
        #     f"{read.query_name},{read.reference_start},{read.reference_end}\n")
        r_id.append(read.query_name)
        r_start.append(read.reference_start)
        r_end.append(read.reference_end)
    return pd.DataFrame({'read_id':r_id,
                         'trans_start':r_start,
                         'trans_end': r_end})
    

def main(args):
    logger.info("Formatting input data...")
    
    # get input data and reformat

    # get jwr from nanosplicer output
    all_reads = parse_format_input_file(args)
    logger.info(helper.mem_time_msg())
    
    # get TSS TTS
    logger.info("Getting TSS and TTS from BAM file")
    tss_tts_d = get_mapped_start_end(args.input_BAM)
    print(tss_tts_d)
    exit()

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
    # correct junction within large groups
    corrected_d = correct_junction_per_group(
                    all_reads, methods=args.nearby_jwr_correction_mode)
    logger.info('Recovering remaining reads using cross-group information...')
    logger.info(helper.mem_time_msg())
    output_h5_file_corrected_reads(corrected_d, args.output_fn, key='data')
    # add some text summary
    add_summary(textwrap.dedent(
        f'''
        After correcting based on nearby JWRs:
            Number of reads with all JWRs corrected: {np.sum(corrected_d.all_corrected)}/{len(corrected_d)}
        '''))

if __name__ == '__main__':
    # print(args)

    all_read, uncorrected_jwr = jwr_correction.main(args)

    all_read.to_hdf(args.output_fn, key='data')
    uncorrected_jwr.to_hdf(args.output_fn, key='uncorrected_jwr')

    # # read uncorrect jwr
    # uncorrected_jwr=pd.read_hdf(args.output_fn, key='uncorrected_jwr')
    # restructure it to read per row and group based on neaby junction
    uncorrected_read = jwr_correction.restructure_per_jwr_dataframe(
                                            uncorrected_jwr,
                                            on='initial_junction',
                                            output_uncorrected=False,
                                            rm_cpl_missed_from_ends=False)
    uncorrected_read.reset_index(drop=False, inplace = True)
    uncorrected_read.to_hdf(args.output_fn, key='uncorrected_read')



# temp input
    # uncorrected_jwr=pd.read_hdf(args.output_fn, key='uncorrected_jwr')
    # uncorrected_read = pd.read_hdf(args.output_fn, key='uncorrected_read')
    
    uncorrected_read['junc_count'] = uncorrected_read.junc.apply(len)
    uncorrected_read = group_reads(uncorrected_read,max_diff = GROUP_ARG['max_diff'])

    
    # correct JWR in each group
    read_id_grp = uncorrected_read.groupby(by=['reference_name',
                                            'non_overlap_group',
                                            'junc_count','sub_group'])['read_id'].apply(list)
    
    corrected_d_list = []
    uncorrected_d_list = []

    for i in read_id_grp:
        # correcting jwr_df per read group
        jwr_to_correct = uncorrected_jwr[uncorrected_jwr.read_id.isin(i)]
        corrected_jwr = jwr_correction.group_and_correct_uncorrected_jwr(args,
                jwr_to_correct, max_diff=CORRECTION_ARG['dist'])
        
        # restructure back to read_df
        corrected_d, uncorrected_d =\
             jwr_correction.restructure_per_jwr_dataframe(
                all_jwr=corrected_jwr,
                on='corrected_junction',
                output_uncorrected=True,
                rm_cpl_missed_from_ends=True,
                summary=False)
                                        
        corrected_d_list.append(corrected_d.drop(columns=['junc_count']))
        uncorrected_d_list.append(uncorrected_d)


    correct_r1 = pd.read_hdf(args.output_fn, key='data')
    correct_r2 = pd.concat(corrected_d_list)
    correct_r2.to_hdf(args.output_fn, key='data2')
    pd.concat(uncorrected_d_list).to_hdf(args.output_fn, key='uncorrected_read_new')
    pd.concat([correct_r1,correct_r2]).to_hdf(args.output_fn, key='corrected_read_new')
    all_read.to_hdf(args.output_fn, key='data')
    #group uncorrected jwr

    # keys in output h5
    # '/corrected_all_jwr', ALL JWR after HC junc correction (contain those with 0 or multiple HC junc)
    # '/uncorrected_jwr',  uncorrected JWR 
    # '/all_read', all reads restructure from all JWR (including all JWRs)
    # '/data', corrected read which all JWRs have only 1 HCJWRs nearby
    # '/uncorrected_read', uncorrected read which all JWRs have only 1 HCJWRs nearby
    # '/data2', corrected read use in second run
    # '/corrected_read_new', all corrected read using 1st and 2nd run
    # '/uncorrected_read_new' all uncorrected read using 1st and 2nd run

    print('\n\nSummary:\n', helper.summary_msg)