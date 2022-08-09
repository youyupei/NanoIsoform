
# default command line input
DEFAULT_INPUT = {
    'JAQ_thres': 0.95, 
    'SIQ_thres': -0.8, 
    'prob_thres': 0.8
}

# grouping parameters
GROUP_ARG = {
    # max distence allowed when grouping the nearby junction
    'max_diff': 30
}


# correst junction based on nearby corrected jwrs
CORRECTION_ARG = {
    # max distence allowed when grouping the nearby junction
    'dist': 20,
    'method': 'majority_vote', # or 'probability'
    # majority voting mode
        # minimum proportion of corrected jwrs required for the major junction
        'maj_vot_min_prop':0.8,
        # minimum count of corrected jwrs required for the major junction
        'maj_vot_min_count': 10
}