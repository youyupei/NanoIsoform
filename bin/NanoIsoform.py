#!/usr/bin/env python3 

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

# this will import variable args and update all global variable from config
from arg_parser import *
import helper
from helper import add_summary
from __restructure_input import * # group reads
from __group_and_correct import *
import jwr_correction

# ACTION REQUIRED FOR THE FINAL VERSION
#   remove levelname, filname, funcName, lineno for the final version
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# set up pandas (temp for debuging)
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

def main_old(args):
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

def main():
    # get TSS TTS for each read
    logger.info("Getting TSS and TTS from BAM file...")
    tss_tts_d = get_mapped_start_end(args.input_BAM)

    # get corrected jwrs (Round 1)
    logger.info("Getting corrected jwrs (Round 1)...")
    corrected_all_jwr = jwr_correction.correction_round1(args, tss_tts_d)
    corrected_all_jwr.to_hdf(args.output_fn, 'corrected_all_jwr')
    # get all reads after first round correction 
    corrected_read_r1, uncorrected_jwr = \
        jwr_correction.restructure_per_jwr_dataframe(
                                            corrected_all_jwr, 
                                            tss_tts_d,
                                            output_uncorrected = True)
    add_summary(f"Round 1: {len(corrected_read_r1)} reads successfully corrected")
    add_summary(f"""Round 1: {len(uncorrected_jwr.read_id.unique())} reads 
                    contains JWR with >1 HC junction nearby, remaining uncorrected""")
    logger.info(
        helper.green_msg(
            f'Round 1 correction is finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}',
            print_out=False))
    
    # add TSS TTS columns
    if True: # save the data
        corrected_read_r1.to_hdf(args.output_fn, key='corrected_read_r1')
    if uncorrected_jwr is not None:
        uncorrected_jwr.to_hdf(args.output_fn, key='uncorrected_jwr_r1')
    else: 
        return None
    # stop if no furthur correction required
    if args.skip_round2_correction:
        return None

    #####################################
    # Round 2: recover uncorrected reads
    #####################################

    # get reads
    uncorrected_read, _ = jwr_correction.restructure_per_jwr_dataframe(
                                            uncorrected_jwr,
                                            tss_tts_d,
                                            on='initial_junction',
                                            output_uncorrected=False,
                                            rm_cpl_missed_from_ends=True)
    uncorrected_read.reset_index(drop=False, inplace = True)
    uncorrected_read.to_hdf(args.output_fn, key='uncorrected_read_r1')
    uncorrected_read['junc_count'] = uncorrected_read.junc.apply(len)
    uncorrected_read = group_reads(uncorrected_read,max_diff = GROUP_ARG['max_diff'])
    
    # group together the reads with same number of junc
    read_id_grp = uncorrected_read.groupby(by=['reference_name',
                                            'non_overlap_group',
                                            'junc_count','sub_group'])['read_id'].apply(list)
    # correct JWR in each group
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
                                                corrected_jwr,
                                                tss_tts_d,
                                                on='corrected_junction',
                                                output_uncorrected=True,
                                                rm_cpl_missed_from_ends=True,
                                                summary=False)
        corrected_d_list.append(corrected_d.drop(columns=['junc_count']))
        uncorrected_d_list.append(uncorrected_d)
    corrected_read_r2 = pd.concat(corrected_d_list)
    corrected_read_r2.to_hdf(args.output_fn, key='corrected_read_r2')
    uncorrected_jwr_r2 = pd.concat(uncorrected_d_list)
    uncorrected_jwr_r2.to_hdf(args.output_fn, key='uncorrected_jwr_r2')

    pd.concat([corrected_read_r1,corrected_read_r2]).to_hdf(args.output_fn, key='corrected_read_all')
    
    add_summary(f"After Round 2: {len(corrected_read_r2) + len(corrected_read_r1)} "
                "reads successfully corrected")
    add_summary(f"After Round 2: {len(uncorrected_jwr_r2.read_id.unique())} "
                "reads contains JWR with >1 HC junction nearby, remaining uncorrected")

    # keys in output h5
    # '/corrected_all_jwr', ALL JWR after HC junc correction (contain those with 0 or multiple HC junc)
    # '/uncorrected_jwr',  uncorrected JWR 
    # '/corrected_read_r1', reads with all JWR corrected, restructured from all JWR
    # '/data', corrected read which all JWRs have only 1 HCJWRs nearby
    # '/uncorrected_read', uncorrected read which all JWRs have only 1 HCJWRs nearby
    # '/data2', corrected read use in second run
    # '/corrected_read_new', all corrected read using 1st and 2nd run
    # '/uncorrected_read_new' all uncorrected read using 1st and 2nd run
    if helper.summary_msg:
        print('\n\nSummary:\n', helper.summary_msg)


if __name__ == '__main__':
    if args.test_mode:
        test()
    else:
        main()

