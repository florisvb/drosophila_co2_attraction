import os, sys
import imp
import numpy as np
import matplotlib.pyplot as plt
import pandas
from optparse import OptionParser
import multi_tracker_analysis as mta

import gc
import time

def get_filename(path, contains, does_not_contain=[]):
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
            fileok = True
            for nc in does_not_contain:
                if nc in filename:
                    fileok = False
            if fileok:
                return os.path.join(path, filename)

def extract_data_to_pd(config, pd, action, flowrate, localtimerange, frames_before, frames_after):
    odor_localtimes, odor_data_times, odor_n_flies = config.get_data_for_action_flowrate_localtimerange(action, flowrate, localtimerange, frames_before, frames_after, 'odor', interpolate_with_time_resolution=0.25, apply_function='sum')
    control_localtimes, control_data_times, control_n_flies = config.get_data_for_action_flowrate_localtimerange(action, flowrate, localtimerange, frames_before, frames_after, 'control', interpolate_with_time_resolution=0.25, apply_function='sum')
    speed_localtimes, speed_times, speed_data = config.get_data_for_action_flowrate_localtimerange(action, flowrate, localtimerange, frames_before, frames_after, 'speed', interpolate_with_time_resolution=0.25, apply_function='mean')
    
    if len(odor_localtimes) > 0:
        for i in range(len(odor_localtimes)):
            
            if config.legacy:
                d = {'action': action, 'flowrate': flowrate, 'n_flies_odor': control_n_flies[i], 'n_flies_control': odor_n_flies[i], 'localtime': odor_localtimes[i], 't': odor_data_times[i], 'identifiercode': config.identifiercode, 'speed': speed_data[i]/ config.pixels_per_mm * config.frames_per_second}
            else:
                if 'N1' in config.identifiercode:
                    d = {'action': action, 'flowrate': flowrate, 'n_flies_odor': odor_n_flies[i], 'n_flies_control': control_n_flies[i], 'localtime': odor_localtimes[i], 't': odor_data_times[i], 'identifiercode': config.identifiercode, 'speed': speed_data[i]/ config.pixels_per_mm * config.frames_per_second}
                elif 'N2' in config.identifiercode: # need to switch odor / control because camera reversed
                    d = {'action': action, 'flowrate': flowrate, 'n_flies_odor': control_n_flies[i], 'n_flies_control': odor_n_flies[i], 'localtime': odor_localtimes[i], 't': odor_data_times[i], 'identifiercode': config.identifiercode, 'speed': speed_data[i]/ config.pixels_per_mm * config.frames_per_second}
            pd = pd.append(d, ignore_index=True)

    return pd    
    
def analyze_n_flies(path, frames_before=18000, frames_after=36000, flowrate_index=0, config=None):
    if config is None:
        configuration_filename = get_filename(path, 'config_', does_not_contain=['~', '.pyc'])
        configuration = imp.load_source('configuration', configuration_filename)
        config = configuration.Config(path)
        config.load_and_process_data()
    
    actions = ['left', 'right']
    flowrates = np.unique(config.odor_control['flowrate'])
    
    localtimerange = [14, 24+11]

    flowrate = flowrates[flowrate_index]
    t = time.time()
    pd = pandas.DataFrame()
    for action in actions:
        print '***: ', action, flowrate
        pd = extract_data_to_pd(config, pd, action, flowrate, localtimerange, frames_before, frames_after)

    fname = config.identifiercode + '_pd_data_' + str(flowrate) + '.pickle'
    fname = os.path.join(config.path, fname)
    pd.to_pickle(fname)
    del(pd)
    print flowrate, time.time()-t
        
    # now merge pandas dataframes
    #pds = []
    #for flowrate in flowrates:
    #    fname = config.identifiercode + '_pd_data_' + str(flowrate) + '.pickle'
    #    fname = os.path.join(config.path, fname)
    #    pds.append(pandas.from_pickle(fname))
    #fname = config.identifiercode + '_pd_data.pickle'
    #fname = os.path.join(config.path, fname)
        
    
    del(config)
    gc.collect()
    
def analyze_circadian(path, config):
    circadian_speed_all = config.data['speed'].groupby(config.data.index).__getattribute__('mean')() 
    circadian_time_all = config.data['time_local'].groupby(config.data.index).__getattribute__('mean')() 
    circadian_time = np.arange(config.data['time_local'].min(), config.data['time_local'].max(), config.circadian_time_resolution) # in hours
    circadian_speed = np.interp(circadian_time, circadian_time_all, circadian_speed_all) / config.pixels_per_mm * config.frames_per_second
    
    pd = pandas.DataFrame()
    pd = pd.append({'circadian_time': circadian_time, 'circadian_speed': circadian_speed}, ignore_index=True)
    fname = config.identifiercode + '_pd_circadian.pickle'
    fname = os.path.join(config.path, fname)
    pd.to_pickle(fname)
    
    del(config)
    gc.collect()
    
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("--path", type="str", dest="path", default='',
                        help="path to data")
    (options, args) = parser.parse_args()  
   
    if 'day' not in options.path:
        paths = mta.read_hdf5_file_to_pandas.get_filenames(options.path, contains='day')
    else:
        paths = [options.path]
    
    for path in paths:
    
    
        configuration_filename = get_filename(path, 'config_', does_not_contain=['~', '.pyc'])
        print 'Path: ', path
        print 'Config filename: ', configuration_filename
        configuration = imp.load_source('configuration', configuration_filename)
        config = configuration.Config(path)
        config.load_and_process_data()
        
        flowrates = np.unique(config.odor_control['flowrate'])
    
        for n in range(len(flowrates)):
            if 'frames_before' in config.__dict__.keys() and 'frames_after' in config.__dict__.keys():
                analyze_n_flies(path, frames_before=config.frames_before, frames_after=config.frames_after, flowrate_index=n, config=config)
            else:
                analyze_n_flies(path, flowrate_index=n, config=config)
        
        # now merge pandas dataframes
        print 'merging'
        pds = []
        for flowrate in flowrates:
            fname = config.identifiercode + '_pd_data_' + str(flowrate) + '.pickle'
            fname = os.path.join(config.path, fname)
            pds.append(pandas.read_pickle(fname))
        pd_all = pandas.concat(pds)
        fname = config.identifiercode + '_pd_data.pickle'
        fname = os.path.join(config.path, fname)
        pd_all.to_pickle(fname)
    
        analyze_circadian(path, config)
    
