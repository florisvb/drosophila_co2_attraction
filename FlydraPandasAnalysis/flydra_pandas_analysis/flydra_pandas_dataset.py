import sys, os
from optparse import OptionParser

import numpy as np
import pickle, imp
import time
import copy
import pandas
import h5py

import matplotlib.pyplot as plt

try:
    import flydra.a2.core_analysis as core_analysis
    import flydra.analysis.result_utils as result_utils
except:
    print 'cannot load flydra.. may not be an issue if it is not needed'
import multi_tracker_analysis as mta

def get_localtime(t):
    lt = time.localtime(t)
    lt_hr = lt.tm_hour + lt.tm_min/60. + lt.tm_sec/3600.
    return lt_hr

def get_filename(path, contains, does_not_contain='none'):
    if type(contains) is not list:
        contains = [contains]
    if type(does_not_contain) is not list:
        does_not_contain = [does_not_contain]
    cmd = 'ls ' + path
    ls = os.popen(cmd).read()
    all_filelist = ls.split('\n')
    try:
        all_filelist.remove('')
    except:
        pass
    filelist = []
    for i, filename in enumerate(all_filelist):
        fok = True
        for c in contains:
            if c not in filename:
                fok = False
        for d in does_not_contain:
            if d in filename:
                fok = False
        if fok:
            filelist.append(os.path.join(path, filename))
    return filelist

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

def load_data(filename, kalman_smoothing=True, dynamic_model=None, fps=None, save_covariance=False):
    # set up analyzer
    ca = core_analysis.get_global_CachingAnalyzer()
    (obj_ids, use_obj_ids, is_mat_file, data_file, extra) = ca.initial_file_load(filename)
    data_file.flush()
    
    # data set defaults
    if fps is None:
        fps = result_utils.get_fps(data_file)
    if dynamic_model is None:
        try:
            dyn_model = extra['dynamic_model_name']
        except:
            print 'cannot find dynamic model'
            print 'using EKF mamarama, units: mm'
            dyn_model = 'EKF mamarama, units: mm'
    if dynamic_model is not None:
        dyn_model = dynamic_model

    # if kalman smoothing is on, then we cannot use the EKF model - remove that from the model name
    print '** Kalman Smoothing is: ', kalman_smoothing, ' **'
    if kalman_smoothing is True:
        dyn_model = dyn_model[4:]
    print 'using dynamic model: ', dyn_model
    print 'framerate: ', fps
    print 'loading data.... '
    
    # load object id's and save as pandas rows
    attributes = {       'objid'                : 'obj_id',
                         'time_epoch'           : 'timestamp',
                         'position_x'           : 'x',
                         'position_y'           : 'y',
                         'position_z'           : 'z',
                         'velocity_x'           : 'xvel',
                         'velocity_y'           : 'yvel',
                         'velocity_z'           : 'zvel',
                         'frames'               : 'frame'}
    pd = None
    first_time = None
    for obj_id in use_obj_ids:
        try:
            print 'processing: ', obj_id
            kalman_rows = ca.load_data( obj_id, data_file,
                             dynamic_model_name = dyn_model,
                             use_kalman_smoothing= kalman_smoothing,
                             frames_per_second= fps)
            
            d = {}
            for attribute, name in attributes.items():
                d.setdefault(attribute, kalman_rows[name].flat)
            index = d['frames']
            
            t = extra['time_model'].framestamp2timestamp(d['frames'][0])
            if first_time is None:
                first_time = t
            else:
                if t - first_time > 24*3600.:
                    'time exceeds 24 hrs'
                    break
            
            pd_tmp = pandas.DataFrame(d, index=index)
            if pd is None:
                pd = pd_tmp
            else:
                pd = pandas.concat([pd, pd_tmp])
        except:
            print 'not enough data: ', obj_id
            
    pd['time_epoch'] = extra['time_model'].framestamp2timestamp(pd.frames)
        
    return pd
    
def process_odor_control_data(path, identifiercode):
    odor_control_bag_filename_contains = identifiercode + '_odor_control.bag'
    print odor_control_bag_filename_contains
    odor_control_bag_filename = get_filename(path, odor_control_bag_filename_contains)[0]

    odor_bag_data = load_bag_as_hdf5(odor_control_bag_filename)
    t = odor_bag_data['alicat_flow_control'].value['t']
    state = odor_bag_data['alicat_flow_control'].value['data']
    
    odor_control = np.vstack((t, state, np.zeros_like(state))).T
    
    # [1,0]: left (neg) side on
    # [0,1]: right (pos) side on 
    
    indices = np.where( odor_bag_data['alicat_flow_control']['data'] != 0 )[0]
    odor_control[indices, 2] = -1
    
    return odor_control
    
def load_and_align_timestamps_for_pickled_pandas_dataframe(filename, identifiercode, odor=True):
    pd = pandas.read_pickle(filename)
    path = os.path.dirname(filename)
    
    if odor:
        odor_control = process_odor_control_data(path, identifiercode)
    else:
        odor_control = None
        
    # get alignment timestamp
    if odor:
        first_odor_row = np.where(odor_control[:,1]>0)[0][0]
        alignment_time = odor_control[first_odor_row][0]
        odor_control[:,0] -= alignment_time
    else:
        first_time = np.min(pd.time_epoch)
        first_time_local = get_localtime(first_time)
        d = 21 - first_time_local
        alignment_time = first_time + d*3600.
    
    # get alignment frame
    dataset = mta.read_hdf5_file_to_pandas.Dataset(pd)
    alignment_frame = dataset.timestamp_to_framestamp(alignment_time)
    
    pd.time_epoch -= alignment_time
    pd.index = pd.frames - alignment_frame
    
    print 'updating objids, may take some time'
    date = identifiercode.split('_')[0]
    pd['objid_int'] = (date + pd.objid.astype('str')).astype('int')
    pd.objid = identifiercode + '_' + pd.objid.astype('str')
    
    
    return pd, odor_control
    
def load_pickled_pandas_dataframes_as_one_dataframe(path, odor=True, days=[]):
    contains = ['trackedobjects.pickle']
    contains.extend(days)
    filenames = get_filename(path, contains)
    
    print filenames
    
    config_filename = get_filename(path, 'config.py')[0]  
    Config = imp.load_source('Config', config_filename)
    config = Config.Config(path)
    
    pds = []
    odor_controls = []
    for f, filename in enumerate(filenames):
        identifiercode = os.path.basename(filename).split('_trackedobjects')[0]
        print identifiercode
        pd, odor_control = load_and_align_timestamps_for_pickled_pandas_dataframe(filename, identifiercode, odor)
        side = config.sides[identifiercode]
        if side == 'left':
            print identifiercode, 'left', 'flipped'
            pd.position_y = pd.position_y*-1
        pd['dataset_code'] = f
        pds.append(pd)
        odor_controls.append(odor_control)
    pd = pandas.concat(pds)
    
    
    config.odor_control = odor_controls[0]
    
    return pd, config

def load_config(path, days=[]):
    contains = ['trackedobjects.pickle']
    contains.extend(days)
    filenames = get_filename(path, contains)
    
    print filenames
    
    config_filename = get_filename(path, 'config.py')[0]  
    Config = imp.load_source('Config', config_filename)
    config = Config.Config(path)
    
    pds = []
    odor_controls = []
    for f, filename in enumerate(filenames):
        identifiercode = os.path.basename(filename).split('_trackedobjects')[0]
        print identifiercode
        odor_control = process_odor_control_data(path, identifiercode)
            
        # get alignment timestamp
        first_odor_row = np.where(odor_control[:,1]>0)[0][0]
        alignment_time = odor_control[first_odor_row][0]
        odor_control[:,0] -= alignment_time
            
        side = config.sides[identifiercode]
        odor_controls.append(odor_control)
    
    config.odor_control = odor_controls[0]
    
    return config
    
'''
def load_pickled_pandas_and_config(path):
    pd = load_pickled_pandas_dataframes_as_one_dataframe(path)
    config_filename = get_filename(path, 'config.py')[0]  
    Config = imp.load_source('Config', config_filename)
    config = Config.Config(path)
    return pd, config
    
'''
    
def get_pd_for_odor_on(pd, config, start_minute=0, end_minute=None, odor_presentations=None):
    odor_on_indices = np.where(config.odor_control[:,1]!=0)[0]
    if odor_presentations is None:
        odor_presentations = np.arange(0, len(odor_on_indices))
    print odor_presentations
        
    times_when_odor_is_on = []
    
    for op in odor_presentations:
        row = config.odor_control[odor_on_indices][op]
        i = odor_on_indices[op]
        if row[1] > 0:
            if end_minute is not None:
                times_when_odor_is_on.append( [config.odor_control[i][0]+start_minute*60, config.odor_control[i][0]+end_minute*60] )
            else:
                print 'all odor on times'
                times_when_odor_is_on.append( [config.odor_control[i][0], config.odor_control[i+1][0]] )
    
    pds = []
    for timerange in times_when_odor_is_on:    
        pd_tmp = mta.data_slicing.get_data_in_epoch_timerange(pd, timerange)
        pds.append(pd_tmp)
    
    pd_odor_on = pandas.concat(pds)
    
    return pd_odor_on

def get_pd_for_odor_off(pd, config, start_minute=0, end_minute=None, odor_presentations=None):
    if config.odor_control[0][1] != 0:
        t = np.min(pd.time_epoch)
        row = np.array([t, 0, 0])
        config.odor_control = np.vstack((row, config.odor_control))
    odor_off_indices = np.where(config.odor_control[:,1]==0)[0]
    
    if odor_presentations is None:
        odor_presentations = np.arange(0, len(odor_off_indices)-1)
    print odor_presentations
        
    times_when_odor_is_off = []
    
    '''
    for i, row in enumerate(config.odor_control[0:-1]):
        if config.odor_control[i][1] == 0:
            print 'odor off'
            times_when_odor_is_off.append( [config.odor_control[i][0], config.odor_control[i+1][0]] )
    '''
    
    for op in odor_presentations:
        row = config.odor_control[odor_off_indices][op]
        i = odor_off_indices[op]
        if row[1] == 0:
            if end_minute is not None:
                times_when_odor_is_off.append( [config.odor_control[i][0]+start_minute*60, config.odor_control[i][0]+end_minute*60] )
            else:
                print 'all odor on times'
                times_when_odor_is_off.append( [config.odor_control[i][0], config.odor_control[i+1][0]] )
                
    pds = []
    for timerange in times_when_odor_is_off:    
        pd_tmp = mta.data_slicing.get_data_in_epoch_timerange(pd, timerange)
        pds.append(pd_tmp)
    
    pd_odor_off = pandas.concat(pds)
    
    return pd_odor_off

def get_pd_for_2_hrs_prior_to_odor_onset(pd, config):
    pd_tmp = mta.data_slicing.get_data_in_epoch_timerange(pd, [-6400,0])
    return pd_tmp
            
            
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--path", type="str", dest="path", default='',
                        help="path to empty data folder, where you have a configuration file")
    parser.add_option("--contains", type="str", dest="contains", default='trackedobjects.h5',
                        help="contains")
    parser.add_option("--does-not-contain", type="str", dest="does_not_contain", default='none',
                        help="does not contain")
    (options, args) = parser.parse_args()
    
    path = options.path    
        
    filenames = get_filename(path, options.contains, options.does_not_contain)
    print 'Found files: ', filenames
    
    for filename in filenames:
        picklename = filename.split('.')[0] + '.pickle'
        if os.path.exists(picklename):
            continue
        else:
            pd = load_data(filename, kalman_smoothing=True, dynamic_model=None, fps=None, save_covariance=False)
            pd.to_pickle(picklename)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
