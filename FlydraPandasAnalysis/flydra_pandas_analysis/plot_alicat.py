import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

import multi_tracker_analysis as mta
import h5py
import time

import fly_plot_lib.plot as fpl 

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
            

def load_bag_as_hdf5(bag, skip_messages=[]):
    output_fname = bag.split('.')[0] + '.hdf5'
    print output_fname
    if not os.path.exists(output_fname):
        mta.bag2hdf5.bag2hdf5(   bag,
                                 output_fname,
                                 max_strlen=200,
                                 skip_messages=skip_messages)    
    metadata = h5py.File(output_fname, 'r')
    return metadata
    
def plot_alicat_data(filename, ax):
    alicat_data = load_bag_as_hdf5(filename)
    
    if 'N1' in filename:
        topic = '/alicat_flow_rate'
    elif 'N2' in filename:
        topic = '/alicat_flow_rate_2'
    else:
        topic = '/alicat_flow_rate'
    
    indices_on = np.where(alicat_data[topic]['data'] > 0.9)[0]
    epoch_time = alicat_data[topic]['t_secs'][indices_on[0]]
    print time.localtime(epoch_time)
    
    t = alicat_data[topic]['t_secs'] - epoch_time
    ax.plot(t, alicat_data[topic]['data'])   
    #add_odor_color_gradients(ax1, config, [0,2])
    #ax1.set_xlim(0, np.max(config.alicat_data[topic]['t_secs']))
    ax.set_ylim(-1, 70)
    
    
