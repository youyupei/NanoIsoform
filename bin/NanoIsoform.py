import os
import logging

from read_correction import nanoisoform_correction_pipeline
import merge_isoform, read_correction, output_formator, helper
from config import *
from arg_parser import *

logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

import pandas as pd
def main():
    # get corrected reads
    corrected_read_all,un_corrected_read = \
        read_correction.nanoisoform_correction_pipeline(args.save_hdf)
    #corrected_read_all = pd.read_hdf('/data/projects/punim0614/yupei/project/NanoIsoform_proj/SIRV_data/NanoIsoform_run_SIRV/2022-12-30_v0.3.1-alpha/2022-1-21_44f5ae2_test1/barcode1_SAP0.95/NanoIsoform_out_run_c3.h5',key='corrected_read_all')
    
    # get 1 pass isoforms
    logger.info("Getting first pass isoforms...")
    first_pass_iso_d = merge_isoform.get_1_pass_isoform(corrected_read_all)

    # filter and merge
    logger.info("Getting final isoforms...")
    isoform_rst_d = merge_isoform.merge_transcript(first_pass_iso_d,
        min_count_to_collapes=MERGE_ARG['min_count'],
        max_ts_dff=MERGE_ARG['max_ts_dff'], 
        max_exbound_ins=MERGE_ARG['max_exbound_ins'])
    
    # output
    ext = os.path.splitext(args.output_fn)[-1].lower()
    if ext in ['.h5', '.hdf5']:
        isoform_rst_d.to_hdf(args.output_fn, 'data')
        helper.green_msg(f"Finished! Output saved in {args.output_fn}!")
    elif ext == '.csv':
        isoform_rst_d.to_csv(args.output_fn)
        helper.green_msg(f"Finished! Output saved in {args.output_fn}!")
    elif ext == '.bed':
        output_formator.df2bed_convertor(isoform_rst_d, args.output_fn)
        helper.green_msg(f"Finished! Output saved in {args.output_fn}!")
    elif ext == '.gtf':
        output_formator.df2gtf_convertor(isoform_rst_d, args.output_fn)
        helper.green_msg(f"Finished! Output saved in {args.output_fn}!")
if __name__ == '__main__':
    main()

