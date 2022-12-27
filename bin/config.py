#####################################################
# default command line input
#####################################################
DEFAULT_INPUT = {
    'JAQ_thres': 0.95, 
    'SIQ_thres': -0.8, 
    'prob_thres': 0.8,
    'output_fn': 'NanoIsoform_out.h5'
}

#####################################################
# HCJWR definition
#####################################################
HCJWR = {
    'JAQ_thres': 0.95, 
    'SIQ_thres': -0.4, 
    'prob_thres': 0.95,
    'min_count':5,
    'consistence': True
    # minimap2 result should also agree with NanoSplicer, 
    # which is hard coded.
}



#####################################################
# grouping parameters
#####################################################
GROUP_ARG = {
    # max distence allowed when grouping the nearby junction
    'max_diff': 10
}

#####################################################
# correst junction based on nearby corrected jwrs
#####################################################
CORRECTION_ARG = {
    # max distence allowed when grouping the nearby junction
    'dist': 10,
    # how to correct jwrs for each junction
    'method': 'majority_vote', # choose between 'majority_vote' and 'probability'
    # majority voting mode
        # minimum proportion of corrected jwrs required for the major junction
        'maj_vot_min_prop':0.8,
        # minimum count of corrected jwrs required for the major junction
        'maj_vot_min_count': 100,
    # probability sampleing mode
        # minimum proportion of corrected jwrs required for any minor junction 
        # to be considered
        'prob_samp_min_prop': 0.1,
        'prob_samp_min_count': 20
}


#####################################################
# output format
#####################################################

# format for the logging (stdout)
# LOG_FORMAT = \
# '[%(asctime)s] %(levelname)s [%(filename)s.%(funcName)s:%(lineno)d] %(message)s'
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'

OUTPUT_CSV = True