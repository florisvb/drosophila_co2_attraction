import os, sys
import imp
import numpy as np
import matplotlib.pyplot as plt
import pandas
from optparse import OptionParser
import multi_tracker_analysis as mta

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

def analyze_n_flies(path, frames_before=18000, frames_after=36000):
    configuration_filename = get_filename(path, 'config_', does_not_contain=['~', '.pyc'])
    configuration = imp.load_source('configuration', configuration_filename)
    config = configuration.Config(path)
    config.load_and_process_data()
    
    actions = ['left', 'right']
    flowrates = np.unique(config.odor_control['flowrate'])
    
    pd = pandas.DataFrame()
    
    localtimerange = [15, 24+11]
    for action in actions:
        for flowrate in flowrates:
            odor_localtimes, odor_data_times, odor_n_flies = config.get_data_for_action_flowrate_localtimerange(action, flowrate, localtimerange, frames_before, frames_after, 'odor', interpolate_with_time_resolution=0.25, apply_function='sum')
            control_localtimes, control_data_times, control_n_flies = config.get_data_for_action_flowrate_localtimerange(action, flowrate, localtimerange, frames_before, frames_after, 'control', interpolate_with_time_resolution=0.25, apply_function='sum')
            speed_localtimes, speed_times, speed_data = config.get_data_for_action_flowrate_localtimerange(action, flowrate, localtimerange, frames_before, frames_after, 'speed', interpolate_with_time_resolution=0.25, apply_function='mean')
            
            if len(odor_localtimes) > 0:
                for i in range(len(odor_localtimes)):
                    if 'N1' in config.identifiercode:
                        d = {'action': action, 'flowrate': flowrate, 'n_flies_odor': odor_n_flies[i], 'n_flies_control': control_n_flies[i], 'localtime': odor_localtimes[i], 't': odor_data_times[i], 'identifiercode': config.identifiercode, 'speed': speed_data[1]}
                    elif 'N2' in config.identifiercode: # need to switch odor / control because camera reversed
                        d = {'action': action, 'flowrate': flowrate, 'n_flies_odor': control_n_flies[i], 'n_flies_control': odor_n_flies[i], 'localtime': odor_localtimes[i], 't': odor_data_times[i], 'identifiercode': config.identifiercode, 'speed': speed_data[1]}
                    pd = pd.append(d, ignore_index=True)
    
    fname = config.identifiercode + '_pd_data.pickle'
    fname = os.path.join(config.path, fname)
    pd.to_pickle(fname)
    
    # circadian
    
    circadian_speed_all = config.data['speed'].groupby(config.data.index).__getattribute__('mean')() 
    circadian_time_all = config.data['time_local'].groupby(config.data.index).__getattribute__('mean')() 
    circadian_time = np.arange(config.data['time_local'].min(), config.data['time_local'].max(), 1/60.) # in hours
    circadian_speed = np.interp(circadian_time, circadian_time_all, circadian_speed_all) / config.pixels_per_mm * config.frames_per_second
    
    pd = pandas.DataFrame()
    pd = pd.append({'circadian_time': circadian_time, 'circadian_speed': circadian_speed}, ignore_index=True)
    fname = config.identifiercode + '_pd_circadian.pickle'
    fname = os.path.join(config.path, fname)
    pd.to_pickle(fname)
    
    del(config)
    
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
        analyze_n_flies(path)
   
    
    
    
    
    
