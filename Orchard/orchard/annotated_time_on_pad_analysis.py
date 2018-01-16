from optparse import OptionParser
import sys, os
import imp

import pickle

import time
import numpy as np
import scipy.stats

from multi_tracker_analysis import read_hdf5_file_to_pandas
from multi_tracker_analysis import data_slicing

import matplotlib.pyplot as plt

import multi_tracker_analysis as mta
import copy

import fly_plot_lib.plot as fpl
import flystat
from fly_plot_lib.colormaps import viridis as viridis

#import three_port_orchard_simplified

import pandas


import data_fit

import co2_paper_locations




def get_filename(path, contains):
    cmd = 'ls ' + path
    ls = os.popen(cmd).read()
    all_filelist = ls.split('\n')
    try:
        all_filelist.remove('')
    except:
        pass
    filelist = []
    for i, filename in enumerate(all_filelist):
        if contains in filename:
            return os.path.join(path, filename)
            
def load_data_and_config(path):
    data_filename = get_filename(path, 'trackedobjects.hdf5')
    config_filename = get_filename(path, 'config')
    
    print data_filename
    print config_filename
    
    data_filename_pickled = data_filename.split('.')[0] + '.pickle'
    try:
        pd = pandas.read_pickle(data_filename_pickled)
    except:
        pd = mta.read_hdf5_file_to_pandas.load_data_as_pandas_dataframe_from_hdf5_file(data_filename)
        pandas.to_pickle(pd, data_filename_pickled)
        
    identifiercode = config_filename.split('config_')[-1].split('.')[0]
    Config = imp.load_source('Config', config_filename)
    config = Config.Config(path, identifiercode)
    
    if config.odor_control[0][1] != 0:
        first_row = [config.odor_control[0][0]-600, 0, 0]
        config.odor_control = np.vstack((first_row, config.odor_control))
    
    return pd, config
    
def get_filelist(path, contains):
    cmd = 'ls ' + path
    ls = os.popen(cmd).read()
    all_filelist = ls.split('\n')
    try:
        all_filelist.remove('')
    except:
        pass
    filelist = []
    for i, filename in enumerate(all_filelist):
        if contains in filename and '~' not in filename:
            filelist.append(os.path.join(path, filename))
    return filelist

def calculate_times_on_pad_for_all_days(directory, include='day'):
    filelist = get_filelist(directory , include)
    print filelist
    for filename in filelist:
        calculate_times_on_pad_for_annotated_trajectories(filename, preprocess=True)

def calculate_times_on_pad_for_annotated_trajectories(path, preprocess=True):
    data_filename = get_filename(path, 'trackedobjects.hdf5')
    if preprocess:
        pd, config = mta.read_hdf5_file_to_pandas.load_and_preprocess_data(data_filename)
        dataset = mta.read_hdf5_file_to_pandas.Dataset(pd)
        annotations_file = open(get_filename(path, 'annotations'))
        annotations = pickle.load(annotations_file)
        annotations_file.close()
        lengths = []
        for key, notes in annotations.items():
            try:
                trajec = dataset.trajec(key)
                l = len(trajec.speed)
                if l > 0:
                    lengths.append( len(trajec.speed) )
            except:
                print 'Failed to load: ', key
        
    else:
        pd, config = load_data_and_config(path)
        pd = mta.read_hdf5_file_to_pandas.cull_short_trajectories(pd, 100)
        keys, lengths = mta.data_slicing.get_nframes_per_key(pd)
        lengths.remove(0)

    fname = os.path.join(path, 'trajectory_lengths.pickle')
    f = open(fname, 'w')
    pickle.dump(lengths, f)
    f.close()
    return lengths

def calculate_mean_speeds_for_annotated_trajectories(path, preprocess=True):
    data_filename = get_filename(path, 'trackedobjects.hdf5')
    pd, config = mta.read_hdf5_file_to_pandas.load_and_preprocess_data(data_filename)
    dataset = mta.read_hdf5_file_to_pandas.Dataset(pd)
    annotations_file = open(get_filename(path, 'annotations'))
    annotations = pickle.load(annotations_file)
    annotations_file.close()
    speeds = []
    for key, notes in annotations.items():
        try:
            trajec = dataset.trajec(key)
            l = len(trajec.speed)
            if l > 0:
                speeds.append( np.mean(trajec.speed) )
        except:
            print 'Failed to load: ', key
        
    fname = os.path.join(path, 'trajectory_speeds.pickle')
    f = open(fname, 'w')
    pickle.dump(speeds, f)
    f.close()
    return speeds

def calculate_speed_distributions_for_annotated_trajectories(path, preprocess=True):
    if 'day' not in path:
        filelist = get_filelist(path, 'day')
    else:
        filelist = [path]

    for path in filelist:
        data_filename = get_filename(path, 'trackedobjects.hdf5')
        pd, config = mta.read_hdf5_file_to_pandas.load_and_preprocess_data(data_filename)
        dataset = mta.read_hdf5_file_to_pandas.Dataset(pd)
        annotations_file = open(get_filename(path, 'annotations'))
        annotations = pickle.load(annotations_file)
        annotations_file.close()
        speeds = []
        for key, notes in annotations.items():
            try:
                trajec = dataset.trajec(key)
                l = len(trajec.speed)
                if l > 0:
                    speeds.extend( trajec.speed.tolist() )
            except:
                print 'Failed to load: ', key
            
        fname = os.path.join(path, 'trajectory_speed_distributions.pickle')
        f = open(fname, 'w')
        pickle.dump(speeds, f)
        f.close()
        

def calculate_number_of_odor_entries(path, preprocess=True):
    data_filename = get_filename(path, 'trackedobjects.hdf5')
    if 1:
        pd, config = mta.read_hdf5_file_to_pandas.load_and_preprocess_data(data_filename)
        
        # Circle Estimates
        # Center (x,y):  380.409250175 240.568325158
        # Radius:  74.6072355957
        if 'co2_ethanol_windy_pad_choice' in path:
            if 'right' in path:
                mta.data_slicing.calc_frames_with_object_in_circular_region(pd, (350, 240), 75, region_name='odor')
            if 'left' in path:
                mta.data_slicing.calc_frames_with_object_in_circular_region(pd, (371, 257), 75, region_name='odor')
        else:
            mta.data_slicing.calc_frames_with_object_in_circular_region(pd, (380.409250175, 240.568325158), 75, region_name='odor')
        dataset = mta.read_hdf5_file_to_pandas.Dataset(pd)
        annotations_file = open(get_filename(path, 'annotations'))
        annotations = pickle.load(annotations_file)
        annotations_file.close()
        entries = []
        for key, notes in annotations.items():
            try:
                trajec = dataset.trajec(key)
            except:
                print 'Failed to load: ', key
                continue
            d = np.diff(trajec.odor)
            entries.append( len(np.where(d==1)[0]) )
    else:
        pass

    fname = os.path.join(path, 'trajectory_odor_entries.pickle')
    f = open(fname, 'w')
    pickle.dump(entries, f)
    f.close()
    return entries

def calculate_number_of_odor_entries_for_all_days(path, include='day'):
    filelist = get_filelist(path , include)
    print filelist
    for filename in filelist:
        calculate_number_of_odor_entries(filename, preprocess=True)
        
def calculate_number_of_odor_entries_for_entire_set():
    paths = co2_paper_locations.data_locations.windtunnel_walking

    for path in paths:
        calculate_number_of_odor_entries_for_all_days(path, include='day')

def calculate_distances_travelled_on_pad_for_annotated_trajectories(path, preprocess=True):
    data_filename = get_filename(path, 'trackedobjects.hdf5')
    if 1:
        pd, config = mta.read_hdf5_file_to_pandas.load_and_preprocess_data(data_filename)
        dataset = mta.read_hdf5_file_to_pandas.Dataset(pd)
        annotations_file = open(get_filename(path, 'annotations'))
        annotations = pickle.load(annotations_file)
        annotations_file.close()
        distances = []
        for key, notes in annotations.items():
            try:
                trajec = dataset.trajec(key)
                dt = 1/30.
                d = np.sum(trajec.speed*dt)
                if d > 0:
                    distances.append( d )
                # single pad: 47.2 pixels / cm
                # left and right: 45.1 pixels / cm
            except:
                print 'Failed to load: ', key
        
    else:
        pass

    fname = os.path.join(path, 'trajectory_distances.pickle')
    f = open(fname, 'w')
    pickle.dump(distances, f)
    f.close()
    return distances
    
def calculate_distances_for_all_days(path, include='day'):
    filelist = get_filelist(path , include)
    print filelist
    for filename in filelist:
        calculate_distances_travelled_on_pad_for_annotated_trajectories(filename, preprocess=True)
        
def calculate_distances(path, include='sccm'):
    filelist = get_filelist(path , include)
    print filelist
    for filename in filelist:
        calculate_distances_for_all_days(filename, include='day')
        
'''
def load_data_for_ethanol_co2_choice_assay(base_path='none', include='day', load_data='lengths'):
    if base_path == 'none':
        base_path = '/media/Orchard3/co2_ethanol_windy_pad_choice'
    filelist = get_filelist(base_path , include)
    print filelist
    lengths = []
    
    eth_only = []
    eth_and_co2 = []
    
    for filename in filelist:
        fname = 'trajectory_' + load_data + '.pickle'
        fname = os.path.join(filename, fname)
        try:
            f = open(fname)
            l = pickle.load(f)
            f.close()
        except:
            calculate_times_on_pad_for_annotated_trajectories(filename)
            calculate_number_of_odor_entries(filename)
            calculate_distances_travelled_on_pad_for_annotated_trajectories(filename)
            f = open(fname)
            l = pickle.load(f)
            f.close()
        
        lengths.extend(l)
        
        config = mta.read_hdf5_file_to_pandas.load_config_from_path(filename)
        a = config.odor_control.value['action']
        idx = np.where(a!='off')
        co2_odor_on_side = a[idx[0]][0]
        
        if '_N1' in config.identifiercode:
            if co2_odor_on_side == 'right':
                eth_and_co2.extend(lengths)
            else:
                eth_only.extend(lengths)
        if '_N2' in config.identifiercode:
            if co2_odor_on_side == 'left':
                eth_and_co2.extend(lengths)
            else:
                eth_only.extend(lengths)
    
    return eth_only, eth_and_co2
'''
    
def load_lengths(include='day', load_data='lengths'):
    def load(path):
        filelist = get_filelist(path , include)
        print filelist
        lengths = []
        for filename in filelist:
            fname = 'trajectory_' + load_data + '.pickle'
            fname = os.path.join(filename, fname)
            print fname
            
            try:
                f = open(fname)
                l = pickle.load(f)
                f.close()
            except:
                functions = {'trajectory_speeds.pickle': calculate_mean_speeds_for_annotated_trajectories,
                            'trajectory_lengths.pickle': calculate_times_on_pad_for_annotated_trajectories,
                            'trajectory_odor_entries.pickle': calculate_number_of_odor_entries,
                            'trajectory_distances.pickle': calculate_distances_travelled_on_pad_for_annotated_trajectories}
                functions[os.path.basename(fname)](filename)
                f = open(fname)
                l = pickle.load(f)
                f.close()
            
            lengths.extend(l)
        return lengths
    
    paths = co2_paper_locations.data_locations.windtunnel_walking
    
    data = {}
    for label, path in paths.items():
        lengths = load(path)
        data.setdefault(label, lengths)
        
    #eth_only, eth_and_co2 = load_data_for_ethanol_co2_choice_assay(load_data=load_data)
    #data['choice-eth-only_60sccm'] = eth_only
    #data['choice-eth-and-co2_60sccm_15sccm'] = eth_and_co2        

    print
    print '*******************'
    print 'Done loading data'
    print '*******************'
    
    return data

'''  
def calculate_nflies_on_pad_histogram():
    data_dict = load_lengths()
    
    fig = plt.figure()
    ax_co2_scatter = fig.add_subplot(221)
    ax_eth_scatter = fig.add_subplot(222)
    ax_co2 = fig.add_subplot(223)
    ax_eth = fig.add_subplot(224)

    for ax in [ax_co2, ax_eth, ax_co2_scatter, ax_eth_scatter]:
        ax.vlines(0, -10, 10, linestyle=':', color='black', linewidth=1)
    
    co2_times, m_co2_mean = three_port_orchard_simplified.plot_max_in_regions_from_directory_experiment_routine_analysis('/media/Orchard2/blacktop_delay_0_leds', portion_of_day='dusk', side=[1,-1], experiment_routine=3, include='day', ax=ax_co2, plot='flies', analysis='max', only_odor_mean=True, normalize_to_prior=True)
    co2_index_peak = np.argmax(m_co2_mean)
    co2_time_peak = co2_times[co2_index_peak]
    
    eth_times, m_eth_mean = three_port_orchard_simplified.plot_max_in_regions_from_directory_experiment_routine_analysis('/media/Orchard2/eth_concentrationsweep_0_5_20', portion_of_day='dusk', side=[1,-1], experiment_routine=5, include='day', ax=ax_eth, plot='flies', analysis='max', only_odor_mean=True, normalize_to_prior=True)
    eth_index_peak = np.argmax(m_eth_mean)
    eth_time_peak = eth_times[eth_index_peak]
    
    def get_bootstrapped_nflies(data_dict, label):
        co2_data = np.array(data_dict[label])
        co2_data.sort()
        # bootstrap confidence intervals
        co2_n_bootstrapped = []
        co2_t_bootstrapped = []
        for i in range(1000):
            indices = np.random.randint(0,len(co2_data), len(co2_data))
            indices.sort()
            co2_n = [len(indices)]
            co2_t = [0]
            for trajec in co2_data[indices]:
                co2_n.append(co2_n[-1]-1)
                co2_t.append(trajec/30./60.)
            co2_n = np.array(co2_n)
            co2_n = co2_n / float(np.max(co2_n))
            co2_t = co2_t + co2_time_peak
            co2_n_bootstrapped.append(co2_n)
            co2_t_bootstrapped.append(co2_t)
        co2_n_bootstrapped_interpolated = []
        co2_t_bootstrapped_interpolated = np.arange(0, int(co2_data[-1]/30./60.)+1,1)
        for i, c in enumerate(co2_n_bootstrapped):
            c_interpolated = np.interp(co2_t_bootstrapped_interpolated, co2_t_bootstrapped[i], co2_n_bootstrapped[i])
            co2_n_bootstrapped_interpolated.append(c_interpolated)
        co2_n_bootstrapped_interpolated = np.array(co2_n_bootstrapped_interpolated)
        co2_n_bootstrapped_interpolated = np.sort(co2_n_bootstrapped_interpolated, axis=0)
        print int(co2_n_bootstrapped_interpolated.shape[0]*0.025), int(co2_n_bootstrapped_interpolated.shape[0]*0.975)
        co2_conf_low = co2_n_bootstrapped_interpolated[ int(co2_n_bootstrapped_interpolated.shape[0]*0.025) ]
        co2_conf_high = co2_n_bootstrapped_interpolated[ int(co2_n_bootstrapped_interpolated.shape[0]*0.975) ]
        return co2_t_bootstrapped_interpolated, co2_n_bootstrapped_interpolated, co2_conf_low, co2_conf_high
        
    co2_t_bootstrapped_interpolated, co2_n_bootstrapped_interpolated, co2_conf_low, co2_conf_high = get_bootstrapped_nflies(data_dict, 'co2_60sccm')
    eth_t_bootstrapped_interpolated, eth_n_bootstrapped_interpolated, eth_conf_low, eth_conf_high = get_bootstrapped_nflies(data_dict, 'eth_60sccm')
    
    ax_co2.plot(co2_t_bootstrapped_interpolated, np.mean(co2_n_bootstrapped_interpolated, axis=0), color='red', linewidth=3)
    ax_co2.fill_between(co2_t_bootstrapped_interpolated, co2_conf_low, co2_conf_high, edgecolor='none', facecolor='red', alpha=0.3, zorder=-100)
    
    ax_eth.plot(eth_t_bootstrapped_interpolated, np.mean(eth_n_bootstrapped_interpolated, axis=0), color='blue', linewidth=3)
    ax_eth.fill_between(eth_t_bootstrapped_interpolated, eth_conf_low, eth_conf_high, edgecolor='none', facecolor='blue', alpha=0.3, zorder=-100)
    
    
    ax_co2.set_xlim(-10,30)
    ax_eth.set_xlim(-10, 30)
    
    ax_co2.set_ylim(0,1.1)
    ax_eth.set_ylim(0, 1.1)
    
    
    ##
    
    color_dict = {'co2_60sccm': 'red',
             'eth_60sccm': 'blue'}
     
    fpl.scatter_box(ax_co2_scatter, 0, (np.array(data_dict['co2_60sccm']) / 30. / 60.).tolist(), color=color_dict['co2_60sccm'], flipxy=True)
    fpl.scatter_box(ax_eth_scatter, 0, (np.array(data_dict['eth_60sccm']) / 30. / 60.).tolist(), color=color_dict['eth_60sccm'], flipxy=True)
    
    ax_co2_scatter.set_xlim(-10,30)
    ax_co2_scatter.set_ylim(-.5,.5)
    ax_eth_scatter.set_xlim(-10,30)
    ax_eth_scatter.set_ylim(-.5,.5)
    
    xticks = [-10,0,10,30]
    
    fpl.adjust_spines(ax_co2_scatter, [], xticks=xticks)
    fpl.adjust_spines(ax_eth_scatter, [], xticks=xticks)
    
    fpl.adjust_spines(ax_co2, ['left', 'bottom'], xticks=xticks, yticks=[0,1])
    fpl.adjust_spines(ax_eth, ['left', 'bottom'], xticks=xticks, yticks=[0,1])
    
    ax_co2.set_xlabel('Time, min')
    ax_co2.set_ylabel('Normalized time on pad\nAnd time near CO2 source')
    
    fig.savefig('/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_co2_mag/time_on_pad_compared_to_walking_data.pdf')
'''

'''
def plot_lengths():
    data_dict = load_lengths()
    color_dict = {'orco_co260': 'black',
                    'co2_5sccm': 'red',
                'co2_15sccm': 'red',
                'co2_60sccm': 'red',
                'co2_hungry': 'red',
             'eth_60sccm': 'magenta',
             'co2_200sccm': 'red',
             'eth_200sccm': 'magenta',
             'eth_15sccm': 'magenta',
             'co2_eth_60sccm': 'red',
             'vinegar_15sccm': 'brown',
             'vinegar_60sccm': 'brown',
             'vinegar_200sccm': 'brown',
             'vinegar_60sccm_15co2': 'red',
             'h2o_60sccm': 'cyan',
             'ethacet_15_50': 'green',
             'ethyl_acetate_15sccm': 'green',
             'ethyl_acetate_60sccm': 'green'}
     
    labels = ['co2_5sccm', 'co2_15sccm', 'co2_60sccm', 'co2_200sccm', 'eth_15sccm', 'eth_60sccm', 'eth_200sccm', 'co2_eth_60sccm', 'vinegar_15sccm', 'vinegar_60sccm', 'vinegar_200sccm', 'vinegar_60sccm_15co2', 'ethacet_15_50', 'ethyl_acetate_15sccm', 'ethyl_acetate_60sccm']
    
    #'vinegar_60sccm_15co2',  'h2o_60sccm', 'orco_co260', 
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    xticks = []
    for x, label in enumerate(labels):
        length = np.array(data_dict[label]) / 30.
        length = length.tolist()
        fpl.scatter_box(ax, x, np.array(length), color=color_dict[label], shading='95conf', markersize=1, edgecolor=color_dict[label])
        xticks.append(x)
    
    stat_analysis = 'median'
    print 'Stats for CO2 60 vs Eth 60: '
    print flystat.resampling.calc_statistical_significance_through_resampling(data_dict['co2_60sccm'], data_dict['eth_60sccm'], iterations=10000, analysis=stat_analysis)
    print
    print 'Stats for CO2 60 vs CO2 200: '
    print flystat.resampling.calc_statistical_significance_through_resampling(data_dict['co2_60sccm'], data_dict['co2_200sccm'], iterations=10000, analysis=stat_analysis)
    
    print 'Stats for Eth 60 vs Eth 200: '
    print flystat.resampling.calc_statistical_significance_through_resampling(data_dict['eth_60sccm'], data_dict['eth_200sccm'], iterations=10000, analysis=stat_analysis)
    
    
    ax.set_ylim(0,3000)
    yticks = [0,600,1200,1800,2400,3000]
    fpl.adjust_spines(ax, ['left', 'bottom'], xticks=xticks, yticks=yticks)
    
    labels = [label.replace('_', '\n') for label in labels]
    ax.set_xticklabels(labels)
    ax.set_ylabel('Time on pad, sec')
    
    fig.savefig('/media/Orchard2/singlepad_windtunnel_blacktop/time_on_pad.pdf')
'''
    
def save_annotated_dataset(path):
    data_filename = get_filename(path, 'trackedobjects.hdf5')
    pd, config = mta.read_hdf5_file_to_pandas.load_and_preprocess_data(data_filename)
    dataset = mta.read_hdf5_file_to_pandas.Dataset(pd)
    annotations_file = open(get_filename(path, 'annotations'))
    annotations = pickle.load(annotations_file)
    annotations_file.close()
    keys = []
    for key, notes in annotations.items():
        try:
            trajec = dataset.trajec(key)
            if len(trajec.speed) > 5:
                keys.append(key)
        except:
            pass
    dataset.keys = keys
    dataset.copy_trajectory_objects_to_dataset()
    fname = os.path.join(path, 'dataset.pickle')
    f = open(fname, 'w')
    pickle.dump(dataset, f)
    f.close()
    
def load_landings_and_takeoffs(path):
    try:
        f = open(get_filename(path, 'dataset.pickle'))
        dataset = pickle.load(f)
        f.close()
    except:
        save_annotated_dataset(path)
        f = open(get_filename(path, 'dataset.pickle'))
        dataset = pickle.load(f)
        f.close()
    landings_x = []
    landings_y = []
    takeoffs_x = []
    takeoffs_y = []
    for key in dataset.keys:
        trajec = dataset.trajecs[key]
        landings_x.append(trajec.position_x[0])
        landings_y.append(trajec.position_y[0])
        
        takeoffs_x.append(trajec.position_x[-1])
        takeoffs_y.append(trajec.position_y[-1])
    return [landings_x, landings_y], [takeoffs_x, takeoffs_y]
            
def plot_landings_and_takeoffs(directory, include='day'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    filelist = get_filelist(directory , include)
    print filelist
    for filename in filelist:
        landings, takeoffs = load_landings_and_takeoffs(filename)
        
        ax.plot(landings[0], landings[1], 'o', color='green', markersize=5)
        ax.plot(takeoffs[0], takeoffs[1], 'o', color='red', markersize=5)
        
    img = plt.imread(get_filename(filename, 'bgimg'))
    ax.imshow(img, cmap='gray', zorder=-100)
    
    ax.set_aspect('equal')
    fpl.adjust_spines(ax, [])
    ax.set_frame_on(False)
    
    fname = os.path.join(directory, 'landings_and_takeoffs.pdf')
    fig.savefig(fname, format='pdf')
    
def get_speed_leading_up_to_takeoff(path, analysis='takeoff'):
    try:
        f = open(get_filename(path, 'dataset.pickle'))
        dataset = pickle.load(f)
        f.close()
    except:
        save_annotated_dataset(path)
        f = open(get_filename(path, 'dataset.pickle'))
        dataset = pickle.load(f)
        f.close()
        
    def is_in_circle(x, y):
        if np.sqrt( (x-384)**2 + (y-238)**2 ) < 200:
            return True
        else:
            return False
    
    mint = 3*60*30    
    speeds = []
    ts = []
    for key in dataset.keys:
        trajec = dataset.trajecs[key]
        if is_in_circle(trajec.position_x[-1], trajec.position_y[-1]):
            t = trajec.time_epoch
            if len(t) > mint:
                
                if analysis == 'takeoff':
                    t -= t[-1]
                    t = t[-1*mint:]
                    s = trajec.speed[-1*mint:]
                    speeds.append(s)
                    ts.append(t)
                elif analysis == 'landing':
                    t -= t[0]
                    t = t[0:mint]
                    s = trajec.speed[0:mint]
                    speeds.append(s)
                    ts.append(t)
        
    
    return ts, speeds
            
def plot_speed_leading_up_to_takeoff(directory, include='day', analysis='takeoff'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    filelist = get_filelist(directory , include)
    print filelist
    ts, speeds = [], []
    resttimes = []
    for filename in filelist:
        t, speed = get_speed_leading_up_to_takeoff(filename, analysis)
        ts.extend(t)
        speeds.extend(speed)
        
        for s in speed:
            indices = np.where(s[0:-6] > 0.2)[0]
            resttimes.append(len(s[0:-6])-indices[-1])
        
    for i, s in enumerate(speeds):
        ax.plot(ts[i], s, color='gray', linewidth=0.5)
        
    s = np.mean(speeds, axis=0)
    t = np.mean(ts, axis=0)
    ax.plot(t, s, color='black', linewidth=3)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(resttimes)
            
            
            
def plot_trajectories(paths, keys, key_lengths=None):
    # dead flies in 200 sccm co2 day 4: 5932, 4394; use lengths: {5932: 9000, 4394: 3000}
    # ethanol trajectory: day2, trajectory 13437 '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_ethanol/day2'
    # h20 trajectory: day2, trajec 299  '/media/Orchard2/singlepad_windtunnel_blacktop/h20_60sccm/day2'
    if type(paths) is not list:
        paths = [paths]
        keys = [keys]
        key_lengths = [key_lengths]
        
    path = paths[0]
    bgimg_filename = get_filename(path, 'bgimg_N1.png')
    bgimg = plt.imread(bgimg_filename)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    ax.imshow(bgimg, cmap='gray')
    
    for p, path in enumerate(paths):
        dataset_filename = os.path.join(path, 'dataset.pickle')
        if not os.path.exists(dataset_filename):
            save_annotated_dataset(path)
        f = open(dataset_filename)
        dataset = pickle.load(f)
        f.close()
        
        for key in keys[p]: 
            trajec = dataset.trajec(key)
            
            if key_lengths is not None:
                l = key_lengths[p][key]
            else:
                l = len(trajec.speed)
            
            interpolated_indices = np.where(trajec.interpolated==1)[0]
            print len(interpolated_indices)
            r = np.arange(0, len(interpolated_indices), 5)
            interpolated_indices = interpolated_indices[r]
            print len(interpolated_indices)
            
            trajec.position_x[interpolated_indices] = np.nan
            trajec.position_y[interpolated_indices] = np.nan
            
            fpl.colorline(ax, trajec.position_x[0:l], trajec.position_y[0:l], trajec.time_epoch[0:l]-trajec.time_epoch[0:l][0], linewidth=2, colormap='none', norm=None, zorder=1, alpha=1, linestyle='solid', cmap=viridis)
    
    fpl.adjust_spines(ax, [])
    ax.set_frame_on(False)
    
    fname = os.path.join(path, 'trajectory_plot.pdf')
    fig.savefig(fname, format='pdf')
        
        
def plot_heatmaps(directory, include='day'):
    # dead flies in 200 sccm co2 day 4: 5932, 4394; use lengths: {5932: 9000, 4394: 3000}
        
        
    filelist = get_filelist(directory , include)
    print filelist
    pds = []
    for path in filelist:
    
        pd, config = mta.read_hdf5_file_to_pandas.load_and_preprocess_data(path)
    
        annotations_file = open(get_filename(path, 'annotations'))
        annotations = pickle.load(annotations_file)
        annotations_file.close()
        
        pd_subset = pd[pd.objid.isin(annotations.keys())]

        if 'interpolated' in pd_subset.keys():
            pd_subset = pd_subset.query('interpolated != 1')
    
        pds.append(pd_subset)
        
    pd = pandas.concat(pds)
    path = filelist[0]
    bgimg_filename = get_filename(path, 'bgimg_N1.png')
    binsx, binsy = mta.plot.get_bins_from_backgroundimage(bgimg_filename)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
        
    vmax = pd.shape[0]*.00003
    mta.plot.plot_heatmap(pd, binsx, binsy, ax=ax, vmin=0, vmax=vmax, logcolorscale=False, position_x='position_x', position_y='position_y')
    
    fpl.adjust_spines(ax, [])
    ax.set_frame_on(False)
    
    fname = os.path.join(directory, 'heatmap.pdf')
    fig.savefig(fname, format='pdf')


def get_cumulative_histogram_of_time_on_pad(lengths):
    t = [0]
    n = [0]
    lengths.sort()
    for l in lengths:
        t.append(l)
        n.append(n[-1]+1)
    
    return np.array(t), np.array(n)
    
def plot_cumulative_histograms(lengths):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    odors = ['co2_60sccm', 'vinegar_60sccm', 'eth_60sccm']
    colors = {'co2_60sccm': 'red', 
              'eth_60sccm': 'blue', 
              'vinegar_60sccm': 'brown',
              }
              
    models = {}
    medians = {}
    
    for o, odor in enumerate(odors):
        t, n = get_cumulative_histogram_of_time_on_pad(lengths[odor])
        t /= 30.
        t_interp = np.linspace(0, 4000, 5000)
        n_interp = np.interp(t_interp, t, n)
        integral = np.sum(n_interp*np.mean(np.diff(t_interp)))
        n = n / float(integral)
        n *= 5000
        print integral, np.max(n)
        
        model = data_fit.models.ExponentialDecay()
        model.parameters['assymptote'] = 0.9*np.max(n)
        model.parameters['gain'] = 0.0001
        model.fit(n, t)              
        models[odor] = model
        medians[odor] = np.median(lengths[odor])
        
        ax.plot(t, n, 'o', color=colors[odor])
        model_val = model.get_val(t_interp)
        ax.plot(t_interp, model_val, color=colors[odor])
        
        fpl.scatter_box(ax, o*0.50-1.30, np.array(lengths[odor])/30., xwidth=.30, ywidth=.7, color=colors[odor], flipxy=True, shading='95conf')
        ax.vlines(np.median(lengths[odor])/30., -200, 200, linestyle=':', color=colors[odor])
    ax.set_ylim(-1.50, 1.50)


    def get_slope(odor, v):
        a = models[odor].parameters['assymptote']
        g = models[odor].parameters['gain']
        return a*g*np.exp(-g*v)
    
    def get_a_for_slope_of_1(odor, v):
        a = models[odor].parameters['assymptote']
        g = models[odor].parameters['gain']
        new_a = 1.18568896993e-11/(g*np.exp(-g*v))
        return new_a
        
    # reference_slope = 1
    slopes = {odor: get_slope(odor, medians[odor]) for odor in odors}
    
    for key, model in models.items():
        print key
        print 'median: ', medians[key]
        print 'old slope: ', get_slope(key, medians[key])
        print 'new a: ', get_a_for_slope_of_1(key, medians[key])
        models[key].parameters['assymptote'] = get_a_for_slope_of_1(key, medians[key])
        models[key].parameters['xoffset'] = 0
        models[key].parameters['yoffset'] = 0
        print 'new slope: ', get_slope(key, medians[key])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    t_interp = np.linspace(0, 4000, 5000)
    
    for key, model in models.items():
        model_val = model.get_val(t_interp)
        ax.plot(t_interp, model_val, color=colors[key])
        
def plot_probability(ax, lengths, o, color):
    lengths.sort()
    lengths = np.array(lengths) / 30.
    
    fpl.scatter_box(ax, o, lengths, xwidth=.15, ywidth=1, color=color, flipxy=True, shading='95conf')

    
    
    #return lengths, y

    model = data_fit.models.LogDecay(fixed_parameters={'b': 100})
    #model.fit(y, lengths)
    model.parameters['a'] = np.median(lengths) / 5387#slope_at_median*(model.fixed_parameters['b']*np.median(lengths)+1) / model.fixed_parameters['b']
    model.parameters['xoffset'] = 0
    model.parameters['yoffset'] = 0
    
    print 'Gaussian stats: ', model.get_val(np.median(lengths)), model.get_val(np.median(lengths))/10.#np.median(lengths)
    print 'Max val: ', model.get_val(3000)
    gnorm = scipy.stats.norm( model.get_val(np.median(lengths)), model.get_val(np.median(lengths))/10.) 
    y = np.array([gnorm.rvs() for i in range(len(lengths))])
    y = y
    y.sort()
    fpl.scatter_box(ax, o*500, np.array(y), xwidth=100, ywidth=1, color=color, flipxy=False, shading='95conf')
    ax.plot(lengths, y, 'o', color=color)
    
    print model.parameters['a'], model.fixed_parameters['b']
    
    def get_slope(model, v):
        a = model.parameters['a']
        b = model.fixed_parameters['b']
        slope = a*b / (b*v + 1)
        return slope
    print 'Slope at median t: ', get_slope(model, np.median(lengths))
    print 'Value at median t: ', model.get_val(np.median(lengths))
    
    t_interp = np.linspace(0, lengths[-1], 1000)
    model_values = model.get_val(t_interp)
    ax.plot(t_interp, model_values, color=color)
    
    ax.vlines(np.median(lengths), -2000, 2000, linestyle=':', color=color)
    ax.hlines(np.median(y), -2000, 2000, linestyle=':', color=color)
    
    print model.get_val(1), model.get_val(0)
    
def plot_all_probabilities():
    lengths = load_lengths()
    
    fig = plt.figure()
    ax = fig.add_subplot(111) 
    
    print 'co2'
    plot_probability(ax, lengths['co2_60sccm'], -.3, 'red')
    
    print
    print 'vinegar'
    plot_probability(ax, lengths['vinegar_60sccm'], -.60, 'brown')
    
    print
    print 'ethanol'
    plot_probability(ax, lengths['eth_60sccm'], -.90, 'blue')
    
    ax.set_ylim(-1,1.5)
    ax.set_xlim(-600, 3000)
    fpl.adjust_spines(ax, ['bottom'])
    fig.savefig('all_odors_optimal_searching.pdf', format='pdf')

def plot_simple_probability():
    lengths = load_lengths()['co2_60sccm']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    color = 'red'
    o = -1
    
    lengths.sort()
    lengths = np.array(lengths) / 30.
    
    fpl.scatter_box(ax, 0.3, lengths, xwidth=0.3, ywidth=1, color=color, flipxy=True, shading='95conf')

    gnorm = scipy.stats.norm( 1, 0.1) 
    y = np.array([gnorm.rvs() for i in range(len(lengths))])
    y = y
    y.sort()
    fpl.scatter_box(ax, -100, np.array(y), xwidth=100, ywidth=1, color=color, flipxy=False, shading='95conf')
    ax.plot(lengths, y, 'o', color=color)
    
    model = data_fit.models.LogDecay()
    model.fit(y, lengths)
    
    t_interp = np.linspace(0, lengths[-1], 1000)
    model_values = model.get_val(t_interp)
    ax.plot(t_interp, model_values, color=color)
    
    ax.vlines(np.median(lengths), -2, 2, linestyle=':', color=color)
    ax.hlines(np.median(y), -100, np.max(lengths), linestyle=':', color=color)
    
    ax.set_ylim(0,1.5)
    fpl.adjust_spines(ax, ['bottom'])
    fig.savefig('co2_logcurve.pdf', format='pdf')
    
    
    
    # now first principles
    lengths = load_lengths()['co2_60sccm']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_probability(ax, lengths, -0.3, 'red')
    ax.set_ylim(-1,1.5)
    ax.set_xlim(-600, 3000)
    fpl.adjust_spines(ax, ['bottom'])
    fig.savefig('co2_logcurve_1stprinciples.pdf', format='pdf')
    
def plot_simple_probability_with_1st_principles_fit():
    lengths = load_lengths()['co2_60sccm']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    color = 'red'
    o = -1
    
    lengths.sort()
    lengths = np.array(lengths) / 30.
    
    fpl.scatter_box(ax, 0.3, lengths, xwidth=0.3, ywidth=1, color=color, flipxy=True, shading='95conf')

    gnorm = scipy.stats.norm( 1, 0.1) 
    y = np.array([gnorm.rvs() for i in range(len(lengths))])
    y = y
    y.sort()
    fpl.scatter_box(ax, -100, np.array(y), xwidth=100, ywidth=1, color=color, flipxy=False, shading='95conf')
    ax.plot(lengths, y, 'o', color=color)
    
    model = data_fit.models.LogDecay(fixed_parameters={'b': 100})
    model.parameters['a'] = 1
    
    t_interp = np.linspace(0, lengths[-1], 1000)
    model_values = model.get_val(t_interp)
    ax.plot(t_interp, model_values, color=color)
    
    ax.vlines(np.median(lengths), -2, 2, linestyle=':', color=color)
    ax.hlines(np.median(y), -100, np.max(lengths), linestyle=':', color=color)
    
    ax.set_ylim(0,1.5)
    fpl.adjust_spines(ax, ['bottom'])
    fig.savefig('co2_logcurve_1stprinciples.pdf', format='pdf')
    
    
    
'''

To make figure:

(0) (optional for length plot) run "save_annotated_dataset" on each path
(1) run "calculate_times_on_pad_for_annotated_trajectories"
(2) run "plot_lengths" to plot lengths for all the directories listed in the function definition



'''
            
            
            
            
            
            
            
            
            
    
