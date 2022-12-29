'''
Parse command line argument and config file.
Note: this function also creats some global parameters.
'''
import argparse
import textwrap
import importlib
import os

from config import *
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
    Dev_arg.add_argument('--skip_round2_correction', action = 'store_true',
                        help= textwrap.dedent(
                            '''Skip second round. correction, which recover uncorrect JWRs in
                            the first round.
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
    helper.check_exist([args.prob_table, args.jwr_check_h5, args.input_BAM])
    return args

args = parse_arg()