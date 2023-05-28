# NanoIsoform: Accurate identification of RNA isoforms using Oxford Nanopore Sequencing
This tool is to perform the downstream analysis on output from [NanoSplicer](https://github.com/shimlab/NanoSplicer), which is software I for accurate identification of splice junction using nanopore data. 

## Keywords:
Oxford Nanopore sequencing, Transcriptomics, Isofrom identification and quantification.

# Installation

```
git clone https://github.com/youyupei/NanoIsoform.git
```
## Conda environment
All dependencies can be installed using "conda", the configuration file can be found at `conda_env/`:

```
conda env create -f conda_env/environment.yml
```

# Running NanoIsoform

## Input file: 
* Probability table from NanoSplicer (TSV)
* NanoSplicer `jwr_checker` output (HDF5)
* BAM file (same as what has been input to NanoSplicer)

```
python3 bin/NanoIsoform.py  < TSV from NanoSplicer>  <HDF5 NanoSplicer> <BAM file>
```

## Output:
BED file. Each entry is a unique isoform and the "score" field (columns 5) represents the read counts for each identified isoforms. 

**For more details:**
```
python3 bin/NanoIsoform.py -h 
```
**Arguments:**
```
usage: NanoIsoform.py [-h] [--SIQ_thres SIQ_THRES] [--prob_thres PROB_THRES]
                      [--JAQ_thres JAQ_THRES]
                      [--nearby_jwr_correction_mode {majority_vote,probability}]
                      [--uniform_prior] [--output_fn OUTPUT_FN]
                      [--cfg_fn CFG_FN] [--hcjwr_SIQ_thres HCJWR_SIQ_THRES]
                      [--hcjwr_prob_thres HCJWR_PROB_THRES]
                      [--hcjwr_JAQ_thres HCJWR_JAQ_THRES]
                      [--hcjwr_consistence HCJWR_CONSISTENCE]
                      [--hcjwr_min_count HCJWR_MIN_COUNT] [--save_hdf]
                      [--no_nearby_jwr_correction] [--skip_round2_correction]
                      prob_table jwr_check_h5 input_BAM

positional arguments:
  prob_table            Filename of the probability table output from
                        NanoSplicer
  jwr_check_h5          Filename of the HDF5 file output from NanoSplicer
                        modulejwr_checker
  input_BAM             Filename of the input BAM (same as what has been input
                        to NanoSplicer module jwr_checker)

optional arguments:
  -h, --help            show this help message and exit
  --SIQ_thres SIQ_THRES
                        SIQ threshold for high qualilty squiggle matching
                        （JWRs with） SIQ < threshold will be ignored from
                        NanoSplicer output. (default: -0.8)
  --prob_thres PROB_THRES
                        The minimum probability of the NanoSplicer identified
                        junction. NanoSplicer identified junction with
                        probability < threshold will be ignored (default: 0.8)
  --JAQ_thres JAQ_THRES
                        Fow JWRs that NanoSplicer does not give identification
                        to. The initial mapped location will be used as
                        "corrected junction" if the Junction Alignment Quality
                        (JAQ) of the initial mapping >= this threshold.
                        (default: 0.95)
  --nearby_jwr_correction_mode {majority_vote,probability}
                        How to correct jwr based on nearby corrected jwrs.
                        'majority_vote': Corrected to most supported junction
                        nearby 'probability': randomly choose from the nearby
                        junctions with probability based on the proportion of
                        support. (default: majority_vote)
  --uniform_prior       Use uniform prior probability for splice patterns.
                        Recommended for synthetic RNA data (e.g. SIRV,
                        Sequins) (default: False)
  --output_fn OUTPUT_FN
                        Output filename, please specify valid extension: .bed
                        (for BED file) or .gtf (for GTF) file, .csv (for CSV
                        file) (default: NanoIsoform_out.bed)
  --cfg_fn CFG_FN       Filename of customised config file. (default: )

Definition of High-confidence JWR:
Note that more stringent threshold than above should be specified. It will be overwritten when less stringet.
:
  --hcjwr_SIQ_thres HCJWR_SIQ_THRES
                        SIQ threshold for high qualilty squiggle matching
                        (default: -0.4)
  --hcjwr_prob_thres HCJWR_PROB_THRES
                        The minimum probability of the NanoSplicer identified
                        junction. (default: 0.95)
  --hcjwr_JAQ_thres HCJWR_JAQ_THRES
                        Junction Alignment Quality (JAQ) of the initial
                        mapping >= this threshold. (default: 0.95)
  --hcjwr_consistence HCJWR_CONSISTENCE
                        NanoSplicer junction is consistent with the minimap2
                        junction. (default: True)
  --hcjwr_min_count HCJWR_MIN_COUNT
                        Minimum support from HCJWRs for each unique HC
                        junction. (default: 5)

For the developers only::
  --save_hdf            Save intermediate dataframe into hdf5 file for
                        debugging. (default: False)
  --no_nearby_jwr_correction
                        For test purpose only: Turn off the correction based
                        on nearby jwr. Note that reads with uncorrected jwr(s)
                        will not appear in the output. (default: False)
  --skip_round2_correction
                        Skip second round. correction, which recover uncorrect
                        JWRs in the first round. (default: False)
```
# Furture improvement: 
* In the current version, reads without strand information is discarded. This can potentially kept just for quantification
