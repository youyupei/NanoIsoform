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
    ref_names, read_ids, transcript_strands, num_junc, juncs,\
    junc_starts, junc_ends, processed, JAQs, high_conf, SIQ, best_prob =\
            [],[],[],[],[],[],[],[],[],[],[],[]
    for key, df in all_jwr.groupby(level=1):
        ref_names.append(df.index.get_level_values(0)[0])
        read_ids.append(df.index.get_level_values(1)[0])
        JAQs.append(tuple(df.JAQ))
        num_junc.append(len(df.JAQ))
        high_conf.append(tuple(df.is_hcjwr))
        transcript_strands.append(df.transcript_strand[0])
        
        identified = tuple(~df.corrected_junction.isnull())
        processed.append(identified)
        mapping = df.index.get_level_values('initial_junction') 
        corrected = df.corrected_junction 
        df_junc = \
            [corrected[i] if identified[i] else mapping[i] for i in range(len(corrected))]
        df_junc = sorted(df_junc)# I think this is a bug
        juncs.append(tuple(df_junc))
        junc_starts.append(tuple([x for x,y in df_junc]))
        junc_ends.append(tuple([y for x,y in df_junc]))
        SIQ.append(tuple(df.SIQ))
        best_prob.append(tuple(df.best_prob))

    all_read = pd.DataFrame({'reference_name':ref_names,  
                            'read_id':read_ids,  
                            'transcript_strand':transcript_strands, 
                            'junc': juncs,
                            'junc_count': num_junc,
                            'junc_start':junc_starts, 
                            'junc_end':junc_ends, 
                            'corrected': processed, 
                            'JAQ':JAQs,
                            'high_conf':high_conf,
                            'SIQ':SIQ,
                            'best_prob': best_prob
    })
    return all_read


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
        d['initial_junction'] = d.initial_junction.apply(format_initial_junction)
        d['candidates'] = d.candidates.apply(format_candidates)
        #d['candidate_sequence_motif_preference'] = d.candidate_sequence_motif_preference.apply(format_tuple_string)
        
        
        
        if args.uniform_prior:
            d['prob'] = d.prob_uniform_prior.apply(format_tuple_string)
        else:
            d['prob'] = d.prob_seq_pattern_prior.apply(format_tuple_string)
            
        d['best_prob'] = d.prob.apply(max)
        d = d[(d.SIQ >= args.SIQ_thres) & (d.best_prob >= args.prob_thres)]
        d['corrected_junction'] = list(d.apply(
            lambda x: x.candidates[np.argmax(x.prob)], axis = 1))

        logger.info(f'Indexing prob table file ...')


        # # corrected to HCJWR
        # corrected_junctions = []
        # hcjwr_candidate_columns = []
        # hcjwr_candidate_probs = []
        
        # # get unique HCJWR and count
        d['is_hcjwr'] = (d.initial_junction == d.corrected_junction) & \
                    (d.JAQ >= HCJWR['JAQ_thres']) & \
                        (d.SIQ >= HCJWR['SIQ_thres']) & \
                        (d.best_prob >= HCJWR['prob_thres'])
        # hcjwr_count_df = d[d['is_hcjwr']][['reference_name', 
        #                     'initial_junction']].value_counts().to_frame('counts').reset_index()
        # hcjwr_count_df = hcjwr_count_df[
        #     hcjwr_count_df.counts > HCJWR['min_count']]


        # # correct non-HC jwr
        # logger.info(f'Correcting JWRs based on High-confidence JWRs ...')
        # for row in tqdm(d.itertuples(), desc='JWRs'):
        #     if row.is_hcjwr:
        #         corrected_junctions.append(row.NanoSplicer_junction)
        #         hcjwr_candidate_columns.append(())
            
        #     else:
        #         # perform corretion
        #         hcjwr_candidate = set(hcjwr_count_df[
        #             hcjwr_count_df.reference_name==row.reference_name].initial_junction
        #             ) & set(row.candidates)

        #         if len(hcjwr_candidate) == 0:
        #             # completely missed jwr
        #             corrected_junctions.append((-1,-1))
        #             hcjwr_candidate_columns.append(())
        #             hcjwr_candidate_probs.append(())

                    
        #         elif len(hcjwr_candidate) == 1:
        #             corrected_junctions.append(hcjwr_candidate.pop())
        #             hcjwr_candidate_columns.append(())
        #             hcjwr_candidate_probs.append(())

        #         else:
        #             hcjwr_candidate = list(hcjwr_candidate)
        #             hcjwr_candidate_columns.append(tuple(hcjwr_candidate))
        #             hcjwr_candidate_probs.append(tuple([row.prob[row.candidates.index(x)] for x in hcjwr_candidate]))
        #             if row.NanoSplicer_junction in hcjwr_candidate:
        #                 corrected_junctions.append(row.NanoSplicer_junction) 
        #             elif row.initial_junction in hcjwr_candidate and \
        #                     row.JAQ > args.JAQ_thres:
        #                 corrected_junctions.append(row.initial_junction) 
        #             else:
        #                 corrected_junctions.append(np.nan) 

        # d['corrected_junction'] = corrected_junctions
        # d['hcjwr_candidates'] = hcjwr_candidate_columns
        # d['hcjwr_candidates_probs'] =  hcjwr_candidate_probs

        d.set_index(
            ['reference_name', 'read_id',  'initial_junction'], inplace=True)
        
        # select which columns in NanoSplicer output to keep here
        d.drop(columns = d.columns.difference(['is_hcjwr', 'corrected_junction', 'candidates', 'prob', 'SIQ', 'best_prob']), inplace = True)
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
        #right_on = ['SIQ', 'best_prob'],
        left_index=True,
        right_index=True)
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    # high-confident junction
    # all_jwr['high_confident_junction'] =\
    #      (all_jwr.index.get_level_values('initial_junction') == all_jwr.corrected_junction) & \
    #      (all_jwr.JAQ >= HCJWR['JAQ_thres']) & \
    #      (all_jwr.SIQ >= HCJWR['SIQ_thres']) & \
    #      (all_jwr.best_prob >= HCJWR['prob_thres'])


    # correct JWRs with JAQ > args.JAQ_thres
    logger.info(f'Correcting JWRs using JAQ threshold ...')
    JAQ_pass = (all_jwr.corrected_junction.isnull()) & (all_jwr.JAQ >= args.JAQ_thres)
    NanoSplicer_skipped = np.sum(all_jwr.corrected_junction.isnull())
    all_jwr.loc[JAQ_pass, 'corrected_junction'] =\
        all_jwr.loc[JAQ_pass,].index.get_level_values('initial_junction')  
    logger.info(f'Finished. Memory used: {helper.check_memory_usage()}, Total runtime:{helper.check_runtime()}')

    # add some text summary
    add_summary(textwrap.dedent(
        f'''
        Total number of JWRs: {len(all_jwr)}
            Number of JWRs without NanoSplicer identification: {NanoSplicer_skipped}
            Number of JWRs uncorrected after JAQ-based correction: {np.sum((all_jwr.corrected_junction.isnull()) & (all_jwr.JAQ < args.JAQ_thres))}
        '''
    ))
    
    # restructure table jwr per row -> read per raw
    logger.info('Restructuring table...')
    all_read = restructure_per_jwr_dataframe(all_jwr)
    print(all_read)
    del all_jwr
    return all_read
