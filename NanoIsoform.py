#!/usr/bin/env python


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
    parser.add_argument('input', type=str,
                        help='Filename of the probability table output from NanoSplicer'
                        'file output from blaze.py')

    # required name argment
    requiredNamed = parser.add_argument_group('Either one of these argument is required')
    requiredNamed.add_argument('--expect-cells',type=int, help='<INT>:  Expected number of cells.')
    requiredNamed.add_argument('--count-threshold', type=int,
                        help='Output the whitelist in Cellranger style')

    # Optional positional argument
    #parser.add_argument('opt_pos_arg', type=int, nargs='?',help)

    # Optional argument
    parser.add_argument('--kit-version', type=str, default='v3',
                        help= textwrap.dedent(
                            '''
                            Choose from v2 and v3 (for 10X Single Cell 3ʹ gene expression v2 or v3). 
                            '''))
    parser.add_argument('--minQ', type=int, default=15,
                        help= textwrap.dedent(
                            '''
                            <INT>: Minimum phred score for all bases in a putative BC. 
                            Reads whose putative BC contains one or more bases with 
                            Q<minQ is not counted in the "Putative BC rank plot".'''))

    parser.add_argument('--full-bc-whitelist', type=str, default=None,
                        help='''<path to file>: .txt file containing all the possible BCs. Users may provide
        their own whitelist. No need to specify this if users want to use the 10X whilelist. The correct version
        of 10X whilelist will be determined based on 10X kit version''')
    parser.add_argument('--out-bc-whitelist', type=str, default=DEFAULT_GRB_OUT_WHITELIST,
                        help='''<filename_prefix>: Output the whitelist identified from all the reads.''')
    parser.add_argument('--cr-style', type=bool, nargs='?',const=True, default=True,
                        help='Output the whitelist in Cellranger style')
    parser.add_argument('--chunk-size', type=int, default=1_000_000,
                        help='Chunksize when reading the input file. Please use'
                        'smaller number if memory is not sufficient.')
    
    args = parser.parse_args()

    if not args.expect_cells and not args.count_threshold:
        helper.err_msg("Missing argument --expect-cells or --count-threshold.") 
        sys.exit(1)
    if args.expect_cells and args.count_threshold:
        helper.warning_msg(textwrap.dedent(
                f'''
                Warning: You have specified both '--expect-cells' and '--count-threshold'. \
'--expect-cells' will be ignored.                
                '''))
    
    args.kit_version = args.kit_version.lower()
    if args.kit_version not in ['v2', 'v3']:
        helper.err_msg("Error: Invalid value of --kit-version, please choose from v3 or v2") 
        sys.exit()

    if args.full_bc_whitelist:
        helper.warning_msg(textwrap.dedent(
                f'You are using {os.path.basename(args.full_bc_whitelist)} as the full barcode'\
                'whitelist. Note that the barcodes not listed in the file will never be found.'))
    else:
        if args.kit_version == 'v3':
            args.full_bc_whitelist = DEFAULT_GRB_WHITELIST_V3
        elif args.kit_version == 'v2':
            args.full_bc_whitelist = DEFAULT_GRB_WHITELIST_V2

    # check file 
    helper.check_exist([args.full_bc_whitelist, args.putative_bc_csv])
    return args



if __name__ == '__main__':
    args = parse_arg()
    # test command line input
    print(args)
    # main(args)
