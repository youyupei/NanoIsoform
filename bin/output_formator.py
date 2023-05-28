"""
Input: NanoIsoform isoform identification and quantification in pd.DataFrame
Output:GTF/BED file

Acknowledgement: 
    The DF to GTF conversion code is modified from a script provided by
    Xinyue Zhao, MSc student from University of Melbourne
"""
    
import argparse
import textwrap
import pickle
import pandas as pd
from collections import defaultdict, Counter
import sys
import os
from tqdm import tqdm
import numpy as np
import multiprocessing
import csv
        
##############################################
# DF to GTF conversion
##############################################

def df2gtf_row_convertor(row, genename, transcript_name):
    Blocksize = len(row.junc) + 1
    juncs = row.junc
    
    # convert from 0-based to 1-based
    isoform = [row.reference_name, 'NanoIsoform', 'transcript', 
            int(row.tss + 1), int(row.tts), 
            row.counts, row.transcript_strand, '.',
            genename, transcript_name]
    
    result = [isoform]
    for i in range(Blocksize):
        if i == 0:
            exon_start = row.tss + 1
            exon_end = juncs[i][0]
        elif i == Blocksize - 1:
            exon_start = juncs[i-1][1] + 1
            exon_end = row.tts
        else:
            exon_start = juncs[i-1][1] + 1
            exon_end = juncs[i][0]
        exonline = [row.reference_name, 'NanoIsoform', 'exon', 
                int(exon_start), int(exon_end), 
                row.counts, row.transcript_strand, '.',
                genename, transcript_name]
        result.append(exonline)
    
    return result


def df2gtf_convertor(df, out_fn, format=None):
    """Convert isoform data frame to GTF or BED file

        format (str): "gtf" or "bed"
    """
    # sort the isoform by tss
    df = df.sort_values(['reference_name', 'tss'], ignore_index = True)
    current_chr = None
    current_tss, current_tts = -1, -1
    gene_idx = 0

    gtf_out = []
    # get gene and transcript name
    for idx, row in enumerate(df.itertuples()):
        if current_chr != row.reference_name:
            current_chr = row.reference_name
            gene_idx = 1
            current_tss, current_tts = row.tss, row.tts
            trans_idx = 1
            genename = f'{current_chr}_G1'
            transcript_name = f'{current_chr}_T1.1'

        elif row.tss < current_tts: # isoform overlap -> same gene new transcript
            trans_idx += 1
            transcript_name = f'{current_chr}_T{gene_idx}.{trans_idx}'
            current_tts = max(row.tts, current_tts)
        else: # isoform do not overlap -> new gene
            gene_idx += 1
            trans_idx = 1
            current_tss, current_tts = row.tss, row.tts
            genename = f'{current_chr}_G{gene_idx}'
            transcript_name = f'{current_chr}_T{gene_idx}.{trans_idx}'

    # get gtf entries 
        result = df2gtf_row_convertor(row, genename, transcript_name)
        gtf_out = gtf_out + result
    inGTF = pd.DataFrame(gtf_out)
    inGTF.columns = ["seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    "gene_id","transcript_id"]
    
    cols=inGTF.columns.tolist()
    df=inGTF[cols[:8]].copy()
    df['attribute']=""
    for c in cols[8:]:
        if c == cols[len(cols)-1]:
            df.loc[:,'attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'";'
        else:
            df.loc[:,'attribute']=df['attribute']+c+' "'+inGTF[c].astype(str)+'"; '
            
    df.to_csv(out_fn, sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)


##############################################
# DF to BED conversion
##############################################
def df2bed_convertor(df, out_fn, color_rgb = '196,196,196'):
    '''
    df: input data frame
    out_fn: output filename
    color_rgb: color of output isoform
    '''
    df = df.sort_values(['reference_name', 'tss'], ignore_index = True)
    current_chr = None
    current_tss, current_tts = -1, -1
    gene_idx = 0
    
    rst = [[] for i in range(12)]
    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = rst
    for idx, row in enumerate(df.itertuples()):

        # find gene name and transcript name
        if current_chr != row.reference_name:
            current_chr = row.reference_name
            gene_idx = 1
            current_tss, current_tts = row.tss, row.tts
            trans_idx = 1
            genename = f'{current_chr}_G1'
            transcript_name = f'{current_chr}_T1.1'

        elif row.tss < current_tts: # isoform overlap -> same gene new transcript
            trans_idx += 1
            transcript_name = f'{current_chr}_T{gene_idx}.{trans_idx}'
            current_tts = max(row.tts, current_tts)
        else: # isoform do not overlap -> new gene
            gene_idx += 1
            trans_idx = 1
            current_tss, current_tts = row.tss, row.tts
            genename = f'{current_chr}_G{gene_idx}'
            transcript_name = f'{current_chr}_T{gene_idx}.{trans_idx}'

        f1.append(row.reference_name)
        f2.append(int(row.tss))
        f3.append(int(row.tts))
        f4.append(transcript_name)
        f5.append(row.counts)
        f6.append(row.transcript_strand)
        f7.append(0)
        f8.append(0)
        f9.append(color_rgb)
        f10.append(len(row.junc) +1)
        if len(row.junc) > 1:
            f11.append(f"{int(row.junc[0][0]-row.tss)},{','.join([str(int(row.junc[x+1][0])-int(row.junc[x][1])) for x in range(len(row.junc)-1)])},{int(row.tts-row.junc[-1][-1])}")
        else:
            f11.append(f"{int(row.junc[0][0]-row.tss)},{int(row.tts-row.junc[-1][-1])}")
        f12.append(f"{0},{','.join([str(int(x[1]-row.tss)) for x in row.junc])}")
    
    rst = pd.DataFrame(rst).transpose()
    # output isoform with count > 0
    rst[rst[4] > 0].to_csv(out_fn , header=False, index=False, sep='\t')
    return rst[rst[4] > 0]