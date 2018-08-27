import fly_plot_lib.plot as fpl

import flystat

import os
import pickle

import multi_tracker_analysis as mta
import numpy as np
import flydra_pandas_dataset
import matplotlib.pyplot as plt


def calc_trajecs_in_volumes(pd):
    
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.06,.1], 
                                                                  [-.07, .07], 
                                                                  [0,.03], 
                                                                  region_name='downwind')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [-.03,.04], 
                                                                  [-.01, .066], 
                                                                  [-.003,.01], 
                                                                  region_name='onpad')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.3,.4], 
                                                                  [-.07, .07], 
                                                                  [-.15,-.1], 
                                                                  region_name='blob')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.31,.51], 
                                                                  [-.07, .07], 
                                                                  [-.03,.06], 
                                                                  region_name='control')

    return pd
    
def get_keys_in_volume(pd, volume):
    q = volume + ' == 1'
    pd_q = pd.query(q)
    return np.unique(pd_q.objid)
    
def calc_stats(directory):
    pd, config = flydra_pandas_dataset.load_pickled_pandas_dataframes_as_one_dataframe(directory)
    
    pd_odor_on = flydra_pandas_dataset.get_pd_for_odor_on(pd, config, 0, 30)
    pd_odor_on = calc_trajecs_in_volumes(pd_odor_on)
    __calc_keys_that_land(pd_odor_on, config, odor_state='on')
    __calc_keys_that_approach(pd_odor_on, config, odor_state='on')
    __calc_keys_near_blob_relative_to_control(pd_odor_on, config, odor_state='on')
    
    pd_odor_off = flydra_pandas_dataset.get_pd_for_odor_off(pd, config, 30, 60)
    pd_odor_off = calc_trajecs_in_volumes(pd_odor_off)
    __calc_keys_that_land(pd_odor_off, config, odor_state='off')
    __calc_keys_that_approach(pd_odor_off, config, odor_state='off')
    __calc_keys_near_blob_relative_to_control(pd_odor_off, config, odor_state='off')
    
def __calc_keys_that_land(pd, config, odor_state='on'):
    q = 'downwind' + ' == 1'
    pd_q = pd.query(q)
    keys_downwind = np.unique(pd_q.objid).tolist()
    
    q = 'onpad' + ' == 1'
    pd_q = pd.query(q)
    keys_onpad = np.unique(pd_q.objid).tolist()
    
    intersection = filter(lambda x: x in keys_downwind, keys_onpad)
    
    # note: checking order does not change anything really
    
    landed = []
    for key in keys_downwind:
        if key in intersection:
            landed.append(1)
        else:
            landed.append(0)
        
    print 100*np.mean(landed), '% landed | odor: ', odor_state
    print 'saving data'
    
    fname = os.path.join(config.path, 'results')
    s = 'landed_' + odor_state + '.pickle'
    fname = os.path.join(fname,  s)
    f = open(fname, 'w')
    pickle.dump(landed, f)
    f.close()
    
def __calc_keys_that_approach(pd, config, odor_state='on'):
    q = 'downwind' + ' == 1'
    pd_q = pd.query(q)
    keys_downwind = np.unique(pd_q.objid)
    
    q = 'control' + ' == 1'
    pd_q = pd.query(q)
    keys_control = np.unique(pd_q.objid)
    
    intersection = filter(lambda x: x in keys_downwind, keys_control)
    
    approached = []
    for key in keys_control:
        if key in intersection:
            approached.append(1)
        else:
            approached.append(0)
        
    print 100*np.mean(approached), '% approached | odor: ', odor_state
    print 'saving data'
    
    fname = os.path.join(config.path, 'results')
    s = 'approached_' + odor_state + '.pickle'
    fname = os.path.join(fname,  s)
    f = open(fname, 'w')
    pickle.dump(approached, f)
    f.close()
    
def __calc_keys_near_blob_relative_to_control(pd, config, odor_state='on'):
    q = 'blob' + ' == 1'
    pd_q = pd.query(q)
    keys_blob = np.unique(pd_q.objid)
    
    q = 'control' + ' == 1'
    pd_q = pd.query(q)
    keys_control = np.unique(pd_q.objid)
    
    intersection = filter(lambda x: x in keys_blob, keys_control)
    
    approached = []
    for key in keys_control:
        if key in intersection:
            approached.append(1)
        else:
            approached.append(0)
        
    print 100*np.mean(approached), '% approached blob | odor: ', odor_state
    print 'saving data'
    
    fname = os.path.join(config.path, 'results')
    s = 'blob_' + odor_state + '.pickle'
    fname = os.path.join(fname,  s)
    f = open(fname, 'w')
    pickle.dump(approached, f)
    f.close()
    
def load_stats(directory, analysis='landed', odor_state='on'):
    
    s = 'results/' + analysis + '_' + odor_state + '.pickle'
    fname = os.path.join(directory, s)
    f = open(fname)
    stats = pickle.load(f)
    f.close()
    
    return np.array(stats)


def make_figure():
    
    data_to_analyze = {'hcs_co2_15sccm' : '/media/caveman/Orchard/rotpad/15percentco2_dry',
                       'hcs_co2_60sccm' : '/media/caveman/Orchard/rotpad/60sccm_co2_dry_hcs',
                       'hcs_eth_15sccm' : '/media/caveman/Orchard/rotpad/hcs_15sccm_ethanol_dry',
                       'hcs_eth_60sccm' : '/media/caveman/Orchard/rotpad/hcs_60sccm_ethanol_dry',
                       'orco_co2_15sccm': '/media/Orchard2/rotpad2/orco_15percentco2_dry',
                       'orco_eth_15sccm': '/media/Orchard2/rotpad2/orco_15sccm_ethanol_dry',
                       'orco_eth_60sccm': '/media/Orchard2/rotpad2/orco_60sccm_ethanol_dry',
                       }
    order = ['hcs_co2_15sccm', 'hcs_co2_60sccm', 'hcs_eth_15sccm', 'hcs_eth_60sccm', 'orco_co2_15sccm', 'orco_eth_15sccm', 'orco_eth_60sccm']
    labels = [s.replace('_', '\n') for s in order]
    
    fig = plt.figure(figsize=(8,8))
    ax_landed = fig.add_subplot(311)
    ax_approached = fig.add_subplot(312)
    ax_blob = fig.add_subplot(313)

    def add_data_to_ax(directory, analysis, odor_state, ax, x):
        if odor_state == 'on':
            if 'co2' in directory:
                color = 'red'
            elif 'ethanol' in directory:
                color = 'blue'
            else:
                color = 'green'
        else:
            color = 'black'
        analysis_on = load_stats(directory, analysis=analysis, odor_state=odor_state)
        analysis_means_on = flystat.resampling.bootstrap_medians_from_data(analysis_on, iterations=1000, use='mean')
        confidence_interval_95 = flystat.resampling.bootstrap_confidence_intervals_for_medians(analysis_means_on, confidence_lo=0.025, confidence_hi=0.975)
        fpl.plot_confidence_interval(ax, x, np.mean(analysis_on), confidence_interval_95, confidence_interval_50=None, width=0.3, color=color, linewidth=3, alpha95=0.3, alpha50=0)
    
    x = 0
    xticks = []
    for label in order:
        directory = data_to_analyze[label]
        
        x += 1
        add_data_to_ax(directory, 'landed', 'on', ax_landed, x)
        add_data_to_ax(directory, 'approached', 'on', ax_approached, x)
        add_data_to_ax(directory, 'blob', 'on', ax_blob, x)
        
        x += 1
        add_data_to_ax(directory, 'landed', 'off', ax_landed, x)
        add_data_to_ax(directory, 'approached', 'off', ax_approached, x)
        add_data_to_ax(directory, 'blob', 'off', ax_blob, x)
        
        xticks.append(x-0.5)
    
    yticks = [0,.1,.2,.3]
    
    fpl.adjust_spines(ax_landed, ['left'], yticks=yticks)
    fpl.adjust_spines(ax_approached, ['left'], yticks=yticks)
    fpl.adjust_spines(ax_blob, ['left', 'bottom'], xticks=xticks, yticks=yticks)
    
    
    ax_blob.set_xticklabels(labels)
    
    ax_landed.set_ylabel('Frac.\nlanded')
    ax_approached.set_ylabel('Frac.\napproached source')
    ax_blob.set_ylabel('Frac.\napproached blob')
    
    fig.savefig( '/media/Orchard2/rotpad2/windtunnel_3d_approach_landing_summary.pdf', format='pdf')
      
'''

pd_odor_on = get_pd_for_odor_on(pd, config, start_minute=0, end_minute=None, odor_presentations=None)
            
            

pds = []
for i in range(-40,70,2):
    pd_subset = fpd.get_pd_for_odor_on(pd, config,i,i+2)
    pds.append(pd_subset)

pos = []
neg = []
blob = []
downwindpad = []
for pd_subset in pds:
    pd_zslice = pd_subset[(pd_subset.position_z>-.002) & (pd_subset.position_z<.01)]
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                              config.circular_region_pos['center'], 
                                                              config.circular_region_pos['radius'], 
                                                              region_name='pos')
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_neg['center'], 
                                                                  config.circular_region_neg['radius'], 
                                                                  region_name='neg')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd_subset, 
                                                              config.rectangular_region_blob['x_range'], 
                                                              config.rectangular_region_blob['y_range'], 
                                                              config.rectangular_region_blob['z_range'], 
                                                              region_name='blob')

'''
