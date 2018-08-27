from optparse import OptionParser

import numpy as np
import multi_tracker_analysis as mta
import os
import h5py
import scipy.signal
import matplotlib.pyplot as plt
import fly_plot_lib.plot as fpl
import fly_plot_lib.text as flytext
from matplotlib.backends.backend_pdf import PdfPages
import time
import matplotlib
import imp
import flydra_pandas_dataset as fpd
import heatmaps

import pickle

def count_total_flies(pd, config):
    pd_zslice = pd[(pd.position_z>-.002) & (pd.position_z<.01)]
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_pos['center'], 
                                                                  config.circular_region_pos['radius'], 
                                                                  region_name='pos')
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_neg['center'], 
                                                                  config.circular_region_neg['radius'], 
                                                                  region_name='neg')
                                                                  
                                                                  
    nflies_in_pos_region = pd_zslice['pos'].groupby(pd_zslice.index).sum()
    nflies_pos_region = np.sum(nflies_in_pos_region)
    nflies_in_neg_region = pd_zslice['neg'].groupby(pd_zslice.index).sum()
    nflies_neg_region = np.sum(nflies_in_neg_region)
    
    pd_9pm_to_7am = mta.data_slicing.get_data_in_epoch_timerange(pd, [0,10*3600.])
    nflies = pd_9pm_to_7am.shape[0]
    
    print 'pos: ',  nflies_pos_region
    print 'neg: ',  nflies_neg_region
    print 'total: ',  nflies
    
    
def subset_preference_index_plot(pd, config, save=''):
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
        mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd_subset, 
                                                                  [.058,.15], 
                                                                  [-.05,.05], 
                                                                  [.01,.06], 
                                                                  region_name='downwindpad')
                                                                      
                                                                      
        
        nflies_in_pos_region = np.sum(pd_zslice['pos'].groupby(pd_zslice.index).sum())
        nflies_in_neg_region = np.sum(pd_zslice['neg'].groupby(pd_zslice.index).sum())
        nflies_in_downwindpad_region = np.sum(pd_subset['downwindpad'].groupby(pd_subset.index).sum())
        nflies_in_blob_region = np.sum(pd_subset['blob'].groupby(pd_subset.index).sum())
        
        pos.append( nflies_in_pos_region )
        neg.append( nflies_in_neg_region )
        blob.append( nflies_in_blob_region )
        downwindpad.append( nflies_in_downwindpad_region )

    normalizer = np.max( np.sum( [pos, neg], axis=0) )
    pos /= normalizer
    neg /= normalizer
    blob /= normalizer
    downwindpad /= normalizer
    downwindpad *= 5
    blob *= 5

    fig = plt.figure(figsize=(8,8))
    ax_pos = fig.add_subplot(4,1,1)
    ax_neg = fig.add_subplot(4,1,2)
    ax_blob = fig.add_subplot(4,1,3)
    ax_downwindpad = fig.add_subplot(4,1,4)
    
    t = range(-40,70,2)
    ax_pos.plot(t, pos)
    ax_neg.plot(t, neg)
    ax_blob.plot(t, blob) 
    ax_downwindpad.plot(t, downwindpad) 
    
    ax_pos.set_ylim(0,1)
    ax_neg.set_ylim(0,1)
    ax_blob.set_ylim(0,1)
    ax_downwindpad.set_ylim(0,1)
    
    ax_pos.fill_between( [0,30], 0, 1, facecolor='red', edgecolor = 'none', zorder=-100, alpha=0.5)
    ax_neg.fill_between( [0,30], 0, 1, facecolor='red', edgecolor = 'none', zorder=-100, alpha=0.5)
    ax_blob.fill_between( [0,30], 0, 1, facecolor='red', edgecolor = 'none', zorder=-100, alpha=0.5)
    ax_downwindpad.fill_between( [0,30], 0, 1, facecolor='red', edgecolor = 'none', zorder=-100, alpha=0.5)
    
    fpl.adjust_spines(ax_pos, ['left', 'bottom'], yticks=[0,1], xticks=[-20,0,30,50])
    fpl.adjust_spines(ax_neg, ['left', 'bottom'], yticks=[0,1], xticks=[-20,0,30,50])
    fpl.adjust_spines(ax_blob, ['left', 'bottom'], yticks=[0,1], xticks=[-20,0,30,50])
    fpl.adjust_spines(ax_downwindpad, ['left', 'bottom'], yticks=[0,1], xticks=[-20,0,30,50])
    
    ax_downwindpad.set_xlabel('time, min')
    
    ax_pos.set_ylabel('N flies near odor\nnormalized')
    ax_neg.set_ylabel('N flies near control\nnormalized')
    ax_blob.set_ylabel('N flies near blob\nnormalized')
    ax_downwindpad.set_ylabel('N flies near pad\nnormalized')

    flytext.set_fontsize(fig, 8)
    
    if len(save) > 0:
        fig.savefig(save, format='pdf')
        
    return pos, neg, blob, downwindpad
    
def count_trajectories_that_landed_or_took_off(pd, config):
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.06,.16], 
                                                                  [-.05,.05], 
                                                                  [.01,.06], 
                                                                  region_name='downwindpad')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [-.05,.05], 
                                                                  [-.075,.075], 
                                                                  [-.001,.01], 
                                                                  region_name='onpad')
    
    onpad = pd.onpad.groupby(pd.objid).max()
    downwindpad = pd.downwindpad.groupby(pd.objid).max() 
    
    candidate_landing_indices = np.where( (onpad.values + downwindpad.values)>1 )[0]
    candidate_landing_keys = [onpad.keys()[i] for i in candidate_landing_indices]
    
    pd_subset = mta.data_slicing.get_pd_subset_from_keys(pd, candidate_landing_keys)
    
    dataset = mta.read_hdf5_file_to_pandas.Dataset(pd_subset)
    keys_landed = []
    for key in candidate_landing_keys:
        print key
        trajec = dataset.trajec(key)
        d = np.where(trajec.downwindpad>0)[0]
        o = np.where(trajec.onpad>0)[0]
        
        if np.max(o) > np.min(d):
            keys_landed.append(key)

    return len(keys_landed)
        
def stats_count_flies(pd, config, save=''):
    pd_zslice = pd[(pd.position_z>-.002) & (pd.position_z<.01)]
    
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_pos['center'], 
                                                                  config.circular_region_pos['radius'], 
                                                                  region_name='pos')
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_neg['center'], 
                                                                  config.circular_region_neg['radius'], 
                                                                  region_name='neg')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  config.rectangular_region_blob['x_range'], 
                                                                  config.rectangular_region_blob['y_range'], 
                                                                  config.rectangular_region_blob['z_range'], 
                                                                  region_name='blob')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.06,.16], 
                                                                  [-.05,.05], 
                                                                  [.01,.06], 
                                                                  region_name='downwindpad')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.5,.6], 
                                                                  [-.05,.05], 
                                                                  [.01,.06], 
                                                                  region_name='control')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [-.05,.05], 
                                                                  [-.075,.075], 
                                                                  [-.001,.01], 
                                                                  region_name='onpad')
    
    def get_keys_for_region(pd, region):                                                              
        control = pd[region].groupby(pd.objid).max()
        indices = np.where( control>0 )[0]
        control = [control.keys()[i] for i in indices]
        return control
        
    control = get_keys_for_region(pd, 'control')
    downwindpad = get_keys_for_region(pd, 'downwindpad')
    pos = get_keys_for_region(pd_zslice, 'pos')
    neg = get_keys_for_region(pd_zslice, 'neg')
    blob = get_keys_for_region(pd, 'blob')
    onpad = get_keys_for_region(pd, 'onpad')
    keylist = {'control': control,
               'downwindpad': downwindpad,
               'pos': pos,
               'neg': neg,
               'blob': blob,
               'onpad': onpad,
               }
    
    all_keys = pd.objid.unique()
    
    data = {'control': [], 
            'downwindpad': [],
            'pos': [],
            'neg': [],
            'blob': [],
            'onpad': [],
            }
    
    for i in range(iterations):
        print 'iteration: ', i
        indices = np.random.randint(0, len(all_keys), len(all_keys))
        keys = all_keys[indices]
        
        for region, regionkeylist in keylist.items():
            n = filter(lambda x: x in keys, regionkeylist)
            print len(n)
            data[region].append(len(n))
         
            
    
    
def count_flies_on_odor_pad_relative_to_total_trajecory_time(pd, config, save='', trajectory_based=True):
    pd_zslice = pd[(pd.position_z>-.002) & (pd.position_z<.01)]
    
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_pos['center'], 
                                                                  config.circular_region_pos['radius'], 
                                                                  region_name='pos')
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_neg['center'], 
                                                                  config.circular_region_neg['radius'], 
                                                                  region_name='neg')
    #mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
    #                                                              config.rectangular_region_blob['x_range'], 
    #                                                              config.rectangular_region_blob['y_range'], 
    #                                                              config.rectangular_region_blob['z_range'], 
    #                                                              region_name='blob')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.33, 0.37], 
                                                                  [-.05,.05], 
                                                                  [-.16,-.12], 
                                                                  region_name='blob')
    
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.065,.105], 
                                                                  [-.05,.05], 
                                                                  [0,.04], 
                                                                  region_name='downwindpad')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [-1,2], 
                                                                  [-.14,.14], 
                                                                  [-.14,.14], 
                                                                  region_name='control')
    #mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
    #                                                              [.06,.16], 
    #                                                              [-.05,.05], 
    #                                                              [.08,.13], 
    #                                                              region_name='control')
    
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [-.05,.05], 
                                                                  [-.075,.075], 
                                                                  [-.001,.01], 
                                                                  region_name='onpad')
                                                                  
    #nflies_in_pos_region = pd_zslice['pos'].groupby(pd_zslice.index).sum()
    #nflies_in_neg_region = pd_zslice['neg'].groupby(pd_zslice.index).sum()
    #nflies_in_downwindpad_region = pd['downwindpad'].groupby(pd.index).sum()
    #nflies_in_blob_region = pd['blob'].groupby(pd.index).sum()
    #nflies_total = len(pd)
    #nflies_in_control = np.sum(pd['control'].groupby(pd.index).sum())
    #nflies = pd.shape[0]
    
    
    if trajectory_based:
        control = np.sum(pd.control.groupby(pd.objid).max())
        downwindpad = np.sum(pd.downwindpad.groupby(pd.objid).max())
        pos = np.sum(pd_zslice.pos.groupby(pd_zslice.objid).max())
        neg = np.sum(pd_zslice.neg.groupby(pd_zslice.objid).max())
        blob = np.sum(pd.blob.groupby(pd.objid).max())
        landed = count_trajectories_that_landed_or_took_off(pd, config)
    else:
        control = np.sum(pd.control)
        downwindpad = np.sum(pd.downwindpad)
        pos = np.sum(pd_zslice.pos)
        neg = np.sum(pd_zslice.neg)
        blob = np.sum(pd.blob)
        landed = np.sum(pd.onpad)
    
    
    #print 'pos / total: ', np.sum(nflies_in_pos_region) / float(nflies_in_control)
    #print 'neg / total: ', np.sum(nflies_in_neg_region) / float(nflies_in_control)
    #print 'downwindpad / total: ', np.sum(nflies_in_downwindpad_region) / float(nflies_in_control)
    #print 'blob / total: ', np.sum(nflies_in_blob_region) / float(nflies_in_control)
    
    
    
    # trajectory based analysis
    #ans = pd['objid'].groupby(pd.downwindpad).unique()
    #nkeys_in_downwindpad = len(ans.values[1])
    #nkeys_not_in_downwindpad = len(ans.values[0])
    #ans = pd_zslice['objid'].groupby(pd_zslice.pos).unique()
    #nkeys_pos_pad = len(ans.values[1])
    #nkeys_not_pos_pad = len(ans.values[0])
    
    #n_trajectories_downwindpad = len(np.unique(pd.iloc[np.where(pd['downwindpad']==1)].objid))
    
    
    #nflies_in_downwindpad_region = pd['downwindpad'].groupby(pd.index).sum()
    
    if 0:
        if len(save) > 0:
            f = open(save, 'w')
            f.writelines('pos / control: ' + str(np.sum(nflies_in_pos_region) / float(np.sum(nflies_in_downwindpad_region))) + '\n')
            f.writelines('neg / control: ' + str(np.sum(nflies_in_neg_region) / float(np.sum(nflies_in_downwindpad_region))) + '\n')
            f.writelines('downwindpad / control: ' + str(np.sum(nflies_in_downwindpad_region) / float(nflies_in_control)) + '\n')
            f.writelines('blob / control: ' + str(np.sum(nflies_in_blob_region) / float(nflies_in_control)) + '\n')
            f.writelines('control / total: ' + str(np.sum(nflies_in_control) / float(nflies_total)) + '\n')
            f.writelines('landed / approached: ' + str(fraction_of_approaches_that_land) + '\n')
            f.close()
    else:
        if len(save) > 0:
            f = open(save, 'w')
            f.writelines('pos / control: ' + str(pos / float(control)) + '\n')
            f.writelines('neg / control: ' + str(neg / float(control)) + '\n')
            f.writelines('downwindpad / control: ' + str(downwindpad / float(control)) + '\n')
            f.writelines('blob / control: ' + str(blob / float(control)) + '\n')
            f.writelines('control: ' + str(control) + '\n')
            f.writelines('landed / approached: ' + str( landed / float(downwindpad) ) + '\n')
            f.close()
            
    counts = {'pos': pos,
              'neg': neg,
              'downwindpad': downwindpad,
              'blob': blob,
              'control': control,
              'landed': landed,
              }
    
    return counts
    
def get_count_statistics(pd, config, trajectory_based=True):
    
    fly_counts_on = None
    fly_counts_off = None
    for f in range(0, np.max(pd.dataset_code)+1):
        pd_f = pd[pd.dataset_code==f]
        
        for op in range(0,8):
            pd_odor_on = fpd.get_pd_for_odor_on(pd_f, config, odor_presentations=[op])
            pd_odor_off = fpd.get_pd_for_odor_off(pd_f, config, 30, 60, odor_presentations=[op])
            
            counts_on = count_flies_on_odor_pad_relative_to_total_trajecory_time(pd_odor_on, config, save='', trajectory_based=trajectory_based )
            counts_off = count_flies_on_odor_pad_relative_to_total_trajecory_time(pd_odor_off, config, save='', trajectory_based=trajectory_based )
            
            print 'DL: ', counts_on['downwindpad'], counts_on['landed']
        
            if fly_counts_on is None:
                fly_counts_on = counts_on
                for key, value in fly_counts_on.items():
                    fly_counts_on[key] = [value]
                fly_counts_on.setdefault('odor_presentation', [op])
                
                fly_counts_off = counts_off
                for key, value in fly_counts_off.items():
                    fly_counts_off[key] = [value]
                fly_counts_off.setdefault('odor_presentation', [op])
                
            else:
                for key, value in counts_on.items():
                    fly_counts_on[key].append(value)
                fly_counts_on['odor_presentation'].append(op)
                for key, value in counts_off.items():
                    fly_counts_off[key].append(value)
                fly_counts_off['odor_presentation'].append(op)
                
    return fly_counts_on, fly_counts_off
    

def calc_pvals(counts_on, counts_off, measurement, controlmeasurement='control', iterations=5000):
    
    value_on = np.array(counts_on[measurement]) / ( np.array(counts_on[controlmeasurement]) )
    value_off = np.array(counts_off[measurement]) / ( np.array(counts_off[controlmeasurement]) )
    
    print np.nanmean(value_on), np.nanmean(value_off)
    
    actual_difference = np.nanmean(value_on - value_off)
    
    all_values = []
    all_values.extend(value_on)
    all_values.extend(value_off)
    all_values = np.array(all_values)
    n = len(all_values)
    
    differences = []
    for i in range(iterations):
        sample1_indices = np.random.randint(0, n, len(value_on))
        sample2_indices = np.random.randint(0, n, len(value_off))
        sample1 = all_values[sample1_indices]
        sample2 = all_values[sample2_indices]
        difference = np.nanmean(sample2 - sample1)    
        differences.append(difference)
        
    differences.sort()
    differences = np.array(differences)
    index = np.argmin( np.abs(differences-actual_difference) )
    
    
    return 1 - np.abs((index - iterations/2.) / (iterations/2.)), value_on, value_off
    
def calc_all_pvals(pd, config, trajectory_based=False):
    data = {}
        
    counts_on, counts_off = get_count_statistics(pd, config, trajectory_based=trajectory_based)
    measurements = ['downwindpad', 'blob', 'pos', 'neg', 'landed']
    
    for measurement in measurements:
        if measurement == 'landed':
            controlmeasurement = 'downwindpad'
        else:
            controlmeasurement = 'control'
        pval, on, off = calc_pvals(counts_on, counts_off, measurement, controlmeasurement=controlmeasurement, iterations=5000)
    
        data.setdefault(measurement, {'pval': pval,
                                      'on_mean': np.nanmean(on),
                                      'on_std': np.nanstd(on),
                                      'on_n': len(on),
                                      'off_mean': np.nanmean(off),
                                      'off_std': np.nanstd(off),
                                      'off_n': len(off),
                                      })
    
    if trajectory_based:
        s = 'results/stats_trajectory_based.pickle'
    else:
        s = 'results/stats.pickle'
    fname = os.path.join(config.path, s)
    f = open(fname, 'w')
    pickle.dump(data, f)
    f.close()
    
    return data
    
def plot_windtunnel_statistics(stat, trajectory_based=False):
    
    paths = ['0percentco2_dry',
             '15percentco2_dry',
             '60sccm_co2_dry_hcs',
             'gr21a_15sccm_co2',
             'hcs_15sccm_ethanol_dry',
             'hcs_60sccm_ethanol_dry',
             'orco_15percentco2_dry',
             'orco_15sccm_ethanol_dry'   ]
    
    labels = {  '15percentco2_dry': '15 CO2',
                 '60sccm_co2_dry_hcs': '60 CO2',
                 'gr21a_15sccm_co2': 'gr21xtnt 15 CO2',
                 'hcs_15sccm_ethanol_dry': '15 EtOH',
                 'hcs_60sccm_ethanol_dry': '60 EtOH',
                 'orco_15percentco2_dry': '15 Orco-',
                 'orco_15sccm_ethanol_dry': '15 EtOH Orco-',
                 '0percentco2_dry': 'Cln Air',
                 }
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for p, path in enumerate(paths):
        fname = os.path.join('/media/caveman/Orchard/rotpad', path)
        if trajectory_based:
            s = 'results/stats_trajectory_based.pickle'
        else:
            s = 'results/stats.pickle'
        fname = os.path.join(fname, s)
        f = open(fname)
        stats = pickle.load(f)
        f.close()
        
        on_mean = stats[stat]['on_mean']
        on_std = stats[stat]['on_std']
        
        off_mean = stats[stat]['off_mean']
        off_std = stats[stat]['off_std']
        
        if 'ethanol' in path:
            color = 'orange'
        else:
            color = 'red'
        
        ax.hlines(off_mean, p-0.3, p, color='black', linewidth=3)
        ax.fill_between([p-0.3, p], off_mean+off_std, off_mean-off_std, facecolor='black', edgecolor='none', alpha=0.5) 
        
        ax.hlines(on_mean, p, p+0.3, color=color, linewidth=3)
        ax.fill_between([p, p+0.3], on_mean+on_std, on_mean-on_std, facecolor=color, edgecolor='none', alpha=0.5)
        
        if stats[stat]['pval'] < 0.01:
            ax.plot(p, 0, 'o', markerfacecolor='green', markeredgecolor='none', markersize=10)
    
    xlabels = [labels[p] for p in paths]
    xticks = range(len(paths))
    
    plt.xticks(rotation=70)
    fpl.adjust_spines(ax, ['left', 'bottom'], xticks=xticks)
    
    ax.set_xticklabels(xlabels)
    
    
def olfactometer_analysis(pd, config, save=''):
    
    pd_zslice = pd[(pd.position_z>-.002) & (pd.position_z<.01)]
    
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_pos['center'], 
                                                                  config.circular_region_pos['radius'], 
                                                                  region_name='pos')
    mta.data_slicing.calc_frames_with_object_in_circular_region(  pd_zslice, 
                                                                  config.circular_region_neg['center'], 
                                                                  config.circular_region_neg['radius'], 
                                                                  region_name='neg')
                                                                  
                                                                  
    nflies_in_pos_region = pd_zslice['pos'].groupby(pd_zslice.index).sum()
    nflies_in_neg_region = pd_zslice['neg'].groupby(pd_zslice.index).sum()
    t = pd_zslice['time_epoch'].groupby(pd_zslice.index).mean()
    indices_sorted = np.argsort(t.values)
    
    a,b = scipy.signal.filter_design.butter(3,0.1)
    
    t = t.iloc[indices_sorted]
    nflies_in_neg_region = nflies_in_neg_region.iloc[indices_sorted]
    nflies_in_neg_region = scipy.signal.filtfilt(a,b,nflies_in_neg_region)
    nflies_in_pos_region = nflies_in_pos_region.iloc[indices_sorted]
    nflies_in_pos_region = scipy.signal.filtfilt(a,b,nflies_in_pos_region)
    
    fig = plt.figure(figsize=(20,8))
    ax_pos = fig.add_subplot(3,1,1)
    ax_neg = fig.add_subplot(3,1,2)
    
    ax_neg.plot(t,nflies_in_neg_region)
    ax_pos.plot(t,nflies_in_pos_region)

    add_odor_color_gradients(ax_pos, config, [0,10])
    ax_pos.set_ylim(0,2)
    add_odor_color_gradients(ax_neg, config, [0,10])
    ax_neg.set_ylim(0,2)
    
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  config.rectangular_region_blob['x_range'], 
                                                                  config.rectangular_region_blob['y_range'], 
                                                                  config.rectangular_region_blob['z_range'], 
                                                                  region_name='blob')
    blob_t = pd['time_epoch'].groupby(pd.index).mean()           
    indices_sorted = np.argsort(blob_t.values)
    nflies_in_blob_region = pd['blob'].groupby(pd.index).sum()
    
    ax_blob = fig.add_subplot(3,1,3)
    
    t = blob_t.iloc[indices_sorted]
    nflies_in_blob_region = nflies_in_blob_region.iloc[indices_sorted]
    nflies_in_blob_region = scipy.signal.filtfilt(a,b,nflies_in_blob_region)
    ax_blob.plot(t,nflies_in_blob_region)

    add_odor_color_gradients(ax_blob, config, [0,10])    
    ax_blob.set_ylim(0,0.1)              
    
    if len(save) > 0:
        fig.savefig(save, format='pdf')                                            
    
###########################################################################################################################

def add_odor_color_gradients(ax, config, ylim):
    odor_control = config.odor_control
    alpha = 0.5
    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    for i, status in enumerate(odor_control[0:-1]):
        if status[1] is True or status[1] > 0:
            if status[2] > 0: # pos side odor
                im = np.linspace(1,-1,50)
                img = im.reshape(len(im), 1)
                ax.imshow(img, norm=norm, cmap=plt.get_cmap('bwr'), extent=[odor_control[i][0], odor_control[i+1][0], ylim[0], ylim[1]], alpha=alpha, zorder=-100, aspect='auto')
            elif status[2] < 0: # neg side odor
                im = np.linspace(-1,1,50)
                img = im.reshape(len(im), 1)
                ax.imshow(img, norm=norm, cmap=plt.get_cmap('bwr'), extent=[odor_control[i][0], odor_control[i+1][0], ylim[0], ylim[1]], alpha=alpha, zorder=-100, aspect='auto')
        else: # no odor
            im = np.linspace(0,-1,50)
            img = im.reshape(len(im), 1)
            midpt = np.mean(ylim)
            ax.imshow(img, norm=norm, cmap=plt.get_cmap('bwr'), extent=[odor_control[i][0], odor_control[i+1][0], ylim[0], midpt], alpha=alpha, zorder=-100, aspect='auto')
            im = np.linspace(-1,0,50)
            img = im.reshape(len(im), 1)
            midpt = np.mean(ylim)
            ax.imshow(img, norm=norm, cmap=plt.get_cmap('bwr'), extent=[odor_control[i][0], odor_control[i+1][0], midpt, ylim[1]], alpha=alpha, zorder=-100, aspect='auto')
        ax.vlines(odor_control[i][0], ylim[0], ylim[1], color='green', linewidth=0.5)



##################################################################################################################################



if __name__ == '__main__':
    ## Read data #############################################################
    parser = OptionParser()
    parser.add_option('--path', type=str, help="the .bag file")
    (options, args) = parser.parse_args()
    
    if options.path=='all':
        paths = ['0percentco2_dry',
                 '15percentco2_dry',
                 '60sccm_co2_dry_hcs',
                 'gr21a_15sccm_co2',
                 'hcs_15sccm_ethanol_dry',
                 'hcs_60sccm_ethanol_dry',
                 'orco_15percentco2_dry',
                 'orco_15sccm_ethanol_dry']
                 
        for path in paths:
            directory = os.path.join('/media/caveman/Orchard/rotpad', path)

            pd, config = fpd.load_pickled_pandas_dataframes_as_one_dataframe(directory)
            
            calc_all_pvals(pd, config, trajectory_based=False)
            calc_all_pvals(pd, config, trajectory_based=True)
    
    else:
    
        pd, config = fpd.load_pickled_pandas_dataframes_as_one_dataframe(options.path)
            
        calc_all_pvals(pd, config, trajectory_based=False)
        calc_all_pvals(pd, config, trajectory_based=True)
    
    if 0:
        pd_odor_on = fpd.get_pd_for_odor_on(pd, config)
        pd_odor_off = fpd.get_pd_for_odor_off(pd, config, 20, 50)
        
        savetodirectory = os.path.join(options.path, 'results')
        if not os.path.exists(savetodirectory):
            os.mkdir(savetodirectory)
        
        heatmaps.plot_xz_heatmap(pd_odor_off, save=os.path.join(savetodirectory, 'odor_off_xz.pdf') )
        heatmaps.plot_xy_heatmap(pd_odor_off, save=os.path.join(savetodirectory, 'odor_off_xy.pdf') )
        
        heatmaps.plot_xz_heatmap(pd_odor_on, save=os.path.join(savetodirectory, 'odor_on_xz.pdf') )
        heatmaps.plot_xy_heatmap(pd_odor_on, save=os.path.join(savetodirectory, 'odor_on_xy.pdf') )
            
        #count_flies_on_odor_pad_relative_to_total_trajecory_time(pd_odor_on, config, save=os.path.join(savetodirectory, 'count_flies_odor_on.txt') )
        #count_flies_on_odor_pad_relative_to_total_trajecory_time(pd_odor_off, config, save=os.path.join(savetodirectory, 'count_flies_odor_off.txt') )
            
            
        
    
    

