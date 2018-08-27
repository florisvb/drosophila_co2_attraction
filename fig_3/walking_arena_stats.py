import make_summary_figure
import multicat_analysis as mcat
import numpy as np

from co2_paper_locations import data_locations
ffname_to_directory = data_locations.walking_arena_activity

import flystat

from multicat_analysis import plot_nflies_vs_speed_scatter


def get_preference_index(ffname, flowrate, time_of_day, average_within_cohorts, use_first_half_of_odor_presentation = True):
    '''
    measure = mean, max, or min
    '''
    directory = ffname_to_directory[ffname]
    response, speed, indices_to_sort, sp, re, experiment_type, t = plot_nflies_vs_speed_scatter.get_new_attraction_index( [directory], 
                                                                                                                          [flowrate], 
                                                                                                                          time_of_day, 
                                                                                                                          average_within_cohorts,
                                                                                                                          use_first_half_of_odor_presentation)
    return re

def full_comparison_of_two_datasets(ffname1, ffname2, flowrate1, flowrate2, time_of_day1, time_of_day2,
                                    average_within_cohorts, use_first_half_of_odor_presentation = True):
    pind1 = get_preference_index(ffname1, flowrate1, time_of_day1, average_within_cohorts, use_first_half_of_odor_presentation)
    pind2 = get_preference_index(ffname2, flowrate2, time_of_day2, average_within_cohorts, use_first_half_of_odor_presentation)

    pvals = []
    for i in range(100):
        r = flystat.resampling.calc_statistical_significance_through_resampling(pind1, pind2, iterations=1000, analysis='mean', pval_for_power_calc=0.05, ax=None)
        pval = r[0]
        power = np.mean(r[1:])
        pvals.append(pval)  
    
    cilo, cihi = flystat.resampling.bootstrap_confidence_intervals_from_data(np.array(pvals), iterations=1000, use='mean', confidence_lo=0.025, confidence_hi=0.975)

    print 'RESULTS:  '
    print
    print pind1
    print pind2
    print
    print
    print 'pval mean: ', np.mean(pvals)
    print '95 conf: ', cilo, cihi

    return np.mean(pvals)

def compare_singlepulse_results(use_first_half_of_odor_presentation = True):
    ffname1 = 'singlepulse_2hrs'
    ffname2 = 'singlepulse_2min'
    flowrate1 = 1
    flowrate2 = 1
    time_of_day1 = 'firstpulse'
    time_of_day2 = [14,36]
    average_within_cohorts = False # irrelevant here b/c one pulse per cohort

    pval = full_comparison_of_two_datasets( ffname1, ffname2, flowrate1, flowrate2, time_of_day1, time_of_day2, average_within_cohorts, 
                                            use_first_half_of_odor_presentation = use_first_half_of_odor_presentation)
    #'singlepulse_2hrs_nflies': 'firstpulse',
    #'singlepulse_2min_nflies': [14,36]
    
def compare_controls(use_first_half_of_odor_presentation = True):
    ffname1 = 'control_dusk'
    ffname2 = 'control_afternoon'
    flowrate1 = 1
    flowrate2 = 1
    time_of_day1 = 'dusk'
    time_of_day2 = 'afternoon'
    average_within_cohorts = False 

    pval = full_comparison_of_two_datasets( ffname1, ffname2, flowrate1, flowrate2, time_of_day1, time_of_day2, average_within_cohorts, 
                                            use_first_half_of_odor_presentation = use_first_half_of_odor_presentation)
    
def compare_control_to_co2_dusk(use_first_half_of_odor_presentation = True):
    ffname1 = 'control_dusk'
    ffname2 = 'low_flow_co2_dusk'
    flowrate1 = 1
    flowrate2 = 1
    time_of_day1 = 'dusk'
    time_of_day2 = 'dusk'
    average_within_cohorts = False 
    pval = full_comparison_of_two_datasets( ffname1, ffname2, flowrate1, flowrate2, time_of_day1, time_of_day2, average_within_cohorts, 
                                            use_first_half_of_odor_presentation = use_first_half_of_odor_presentation)


def compare_co2_to_ethanol_starved_dusk(use_first_half_of_odor_presentation = False):
    ffname1 = 'low_flow_co2_starved_dusk'
    ffname2 = 'ethanol_starved_dusk'
    flowrate1 = 1
    flowrate2 = 1
    time_of_day1 = 'dusk'
    time_of_day2 = 'dusk'
    average_within_cohorts = False 

    pval = full_comparison_of_two_datasets( ffname1, ffname2, flowrate1, flowrate2, time_of_day1, time_of_day2, average_within_cohorts, 
                                            use_first_half_of_odor_presentation = use_first_half_of_odor_presentation)

def compare_co2_to_ethanol_dusk(use_first_half_of_odor_presentation = False):
    ffname1 = 'low_flow_co2_dusk'
    ffname2 = 'ethanol_dusk'
    flowrate1 = 1
    flowrate2 = 1
    time_of_day1 = 'dusk'
    time_of_day2 = 'dusk'
    average_within_cohorts = False 

    pval = full_comparison_of_two_datasets( ffname1, ffname2, flowrate1, flowrate2, time_of_day1, time_of_day2, average_within_cohorts, 
                                            use_first_half_of_odor_presentation = use_first_half_of_odor_presentation)



