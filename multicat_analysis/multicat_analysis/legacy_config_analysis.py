import imp
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import itertools

import multi_tracker_analysis as mta
import h5py

import cv2, time

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

class Config(object):
    def __init__(self, path, identifiercode=None):
        if '.py' in path:
            self.path = os.path.dirname(path)
        else:
            self.path = path
        if identifiercode is None:
            config_filename = os.path.basename( get_filename(self.path, 'config_', ['~', '.pyc']) )
            self.identifiercode = config_filename.split('config_')[1].split('.py')[0]
        else:
            self.identifiercode = identifiercode
        self.legacy = True
        self.preprocess_data_function = self.preprocess_data
        self.postprocess_data_function = self.postprocess_data
        self.nflies = 10
        self.xlim_nflies = [-600, 1200]
        
        self.pixels_per_mm = 331/77.
        self.frames_per_second = 30.
        self.center_xpixel = 335
        
        self.circular_region_pos = {'center': [558,245], 'radius': 80}
        self.circular_region_neg = {'center': [116,245], 'radius': 80}
        self.circular_regions = {'pos': self.circular_region_pos,
                                 'neg': self.circular_region_neg,
                                 }
        
        self.lights_on = 7 # localtime
        self.daylength = 16 # hours
        self.circadian_time_resolution = 1/60.
        
        self.portion_of_day_to_local_time_range = {'afternoon': [15, 20.6],
                                                   'dusk': [20.6,25],
                                                   'night': [25,24+7+.4],
                                                   'morning': [24+7+.4, 24+11],
                                                   }
        
        self.skip_odor_control_messages=[]
        self.skip_alicat_messages=[]
        
    def load_and_process_data(self, preprocess=True, postprocess=True):
        contains = self.identifiercode + '_trackedobjects.hdf5'
        tracked_objects_filename = get_filename(self.path, contains)
        self.data = mta.read_hdf5_file_to_pandas.load_data_as_pandas_dataframe_from_hdf5_file(tracked_objects_filename)
        if preprocess:
            self.data = self.preprocess_data_function(self.data)
        if postprocess:
            self.data = self.postprocess_data_function(self.data)
    
    def get_keyframe_times(self, action, flowrate):
        indices_action = np.where(self.odor_control['action']==action)
        indices_flowrate = np.where(self.odor_control['flowrate']==flowrate)
        indices_ok = np.intersect1d(indices_flowrate[0], indices_action[0])
        keyframe_times = [self.odor_control['t'][i] for i in indices_ok]
        return keyframe_times
        
    def get_alicat_for_action_flowrate_localtimerange(self, unit, action, flowrate, seconds_before, seconds_after):
        if '_N1' in self.identifiercode:
            unit = 1
        elif '_N2' in self.identifiercode:
            unit = 2  
        key = 'alicat_flow_rate_' + unit
        keyframe_times = self.get_keyframe_times(action, flowrate)
        
        flowdata = []
        for keyframe_time in keyframe_times:
            first_index = np.argmin( np.abs( self.alicat_data[key]['t'] - (keyframe_time-seconds_before) ) )
            last_index = np.argmin( np.abs( self.alicat_data[key]['t'] - (keyframe_time+seconds_after) ) )
            d = self.alicat_data[key]['data']
            flowdata.append(d[first_index:last_index])
        t = self.alicat_data[key]['t'][first_index:last_index]-keyframe_time
    
        return t, flowdata
        
    
    def get_data_for_action_flowrate_localtimerange(self, action, flowrate, localtimerange, frames_before, frames_after, key, interpolate_with_time_resolution=None, apply_function='sum'):
        keyframe_times = self.get_keyframe_times(action, flowrate)
        
        data_for_key_at_times = []
        data_times = []
        localtimes = []
        
        if key == 'odor':
            if action=='right':
                key = 'pos'
            elif action=='left':
                key = 'neg'
        elif key == 'control':
            if action=='right':
                key = 'neg'
            elif action=='left':
                key = 'pos'
                
        length_should_be = frames_after+frames_before-1
        length_in_seconds_should_be = (frames_after + frames_before)/self.frames_per_second
        
        localtimerange_query = 'time_local > ' + str(localtimerange[0]) + ' and time_local < ' + str(localtimerange[1])
        data_for_localtimerange = self.data.query(localtimerange_query)
        for keyframe_time in keyframe_times:
            pandas_index = np.argmin( np.abs(data_for_localtimerange.time_epoch - keyframe_time) ) 
            try:
                frame = data_for_localtimerange.frames.loc[pandas_index].median()
            except AttributeError:
                frame = data_for_localtimerange.frames.loc[pandas_index]
                print frame
            #    frame = data_for_localtimerange.frames.loc[pandas_index]
            
            query = 'frames > ' + str(frame-frames_before) + ' and frames < ' + str(frame+frames_after)
            pandas_slice = data_for_localtimerange.query(query)
            data_for_key_at_time = pandas_slice[key].groupby(pandas_slice.index).__getattribute__(apply_function)()
            data_time = pandas_slice['time_epoch'].groupby(pandas_slice.index).mean()
            
            #if len(data_time) < length_should_be or len(data_for_key_at_time) < length_should_be:
            #    print len(data_time), len(data_for_key_at_time), length_should_be
            #    continue # must have picked a stretch that saddles the localtimerange
            
            if interpolate_with_time_resolution is not None:
                new_time_base = np.arange(data_time.min(), data_time.min()+length_in_seconds_should_be, interpolate_with_time_resolution) # length_in_seconds_should_be term used to fix single fly data, which sometimes has a missing data point
                data_for_key_at_time = np.interp(new_time_base, data_time, data_for_key_at_time.values)
                data_times.append(new_time_base-keyframe_time)
                data_for_key_at_times.append(data_for_key_at_time)
            else:
                data_times.append(data_time-keyframe_time)
                data_for_key_at_times.append(data_for_key_at_time.values)
                
            if (data_time-keyframe_time).min() < -1*(frames_before/self.frames_per_second):
                continue # something wrong
            
            else:
                localtimes.append( pandas_slice['time_local'].groupby(pandas_slice.index).mean().min() )
            
        return localtimes, data_times, data_for_key_at_times
    
    def process_odor_control_data(self):
        odor_control_bag_filename_contains = self.identifiercode + '_led_and_odor_control.bag'
        odor_control_bag_filename = get_filename(self.path, odor_control_bag_filename_contains)
        self.odor_bag_data = load_bag_as_hdf5(odor_control_bag_filename, skip_messages=self.skip_odor_control_messages)
        if '_N1' in self.identifiercode:
            rig = 'rig1'
        elif '_N2' in self.identifiercode:
            rig = 'rig2'  
        self.odor_control = self.odor_bag_data['alicat_flow_control']
    
    def calculate_actions_and_flowrates(self):
        # calculate actions
        actions = []
        flowrates = []
        seq = []
        for r, row in enumerate(self.odor_bag_data['phidgets_daq/digital_output'].value):
        
            # check for shift
            shift = 0
            if (self.odor_bag_data['phidgets_daq/digital_output'][0][2] == True) or (self.odor_bag_data['phidgets_daq/digital_output'][0][3] == True):
                if self.odor_bag_data['alicat_flow_control'].value[0][0] == 0:
                    shift = 1
        
            if 'flowrate_override' in self.__dict__.keys():
                if self.odor_bag_data['alicat_flow_control'].value[r+shift][0] != 0:
                    flowrates.append( self.flowrate_override )
                else:
                    flowrates.append( 0 )
            else:
                flowrates.append( self.odor_bag_data['alicat_flow_control'].value[r+shift][0] )
            seq.append(r)
            
            if row['states_0'] == False and row['states_1'] == False:
                actions.append('off')
            elif row['states_0'] == False and row['states_1'] == True:
                actions.append('left')
            elif row['states_0'] == True and row['states_1'] == False:
                actions.append('right')
        
        shape = self.odor_bag_data['phidgets_daq/digital_output'].value.shape
        new_odor_control = np.recarray(shape, dtype=[('header_seq', '<u4'), ('header_stamp_secs', '<i4'), ('header_stamp_nsecs', '<i4'), ('header_frame_id', 'S200'), ('header_stamp', '<f8'), ('action', 'S200'), ('flowrate', '<f4'), ('t_secs', '<i4'), ('t_nsecs', '<i4'), ('t', '<f8')])
        
        new_odor_control['action'] = actions
        new_odor_control['header_frame_id'] = ['' for i in range(shape[0])]
        new_odor_control['flowrate'] = flowrates
        new_odor_control['t'] = self.odor_bag_data['phidgets_daq/digital_output'].value['t']
        new_odor_control['t_secs'] = self.odor_bag_data['phidgets_daq/digital_output'].value['t_secs']
        new_odor_control['t_nsecs'] = self.odor_bag_data['phidgets_daq/digital_output'].value['t_nsecs']
        
        new_odor_control['header_seq'] = seq
        new_odor_control['header_stamp'] = self.odor_bag_data['phidgets_daq/digital_output'].value['t']
        
        new_odor_control['header_stamp_secs'] = self.odor_bag_data['phidgets_daq/digital_output'].value['t_secs']
        new_odor_control['header_stamp_nsecs'] = self.odor_bag_data['phidgets_daq/digital_output'].value['t_nsecs']
        
        self.odor_control = new_odor_control
            
    def process_alicat_data(self):
        alicat_bag_filename_contains = self.identifiercode + '_alicat.bag'
        alicat_bag_filename = get_filename(self.path, alicat_bag_filename_contains)
        self.alicat_data = load_bag_as_hdf5(alicat_bag_filename, skip_messages=self.skip_alicat_messages)
    
    def get_odor_status(self, epoch_time):
        if epoch_time < self.odor_control_time[0]:
            return self.odor_control_time[0]
        d = epoch_time - self.odor_control_time
        index = np.where(d>0)[0][-1]
        return self.odor_control_action[index]
        
    def draw(self, img, epoch_time):
        odor_status = self.get_odor_status(epoch_time)
        
        if odor_status == 'off':
            #color_pos = [255,255,0]
            #color_neg = [255,255,0]
            color_neg = [0,255,255]
            color_pos = [0,255,255]
        else:
            if odor_status == 'right':
                #color_pos = [0,0,255]
                color_pos = [255,0,0]
                #color_neg = [255,255,0]
                color_neg = [0,255,255]
            elif odor_status == 'left':
                #color_pos = [255,255,0]
                color_pos = [0,255,255]
                #color_neg = [0,0,255]
                color_neg = [255,0,0]
        
        cv2.circle(img, ( int(self.circular_region_pos['center'][0]), int(self.circular_region_pos['center'][1])), int(self.circular_region_pos['radius']), color_pos, 2)
        cv2.circle(img, ( int(self.circular_region_neg['center'][0]), int(self.circular_region_neg['center'][1])), int(self.circular_region_neg['radius']), color_neg, 2)
        
        localtime_hour = time.localtime(epoch_time).tm_hour
        localtime_min = time.localtime(epoch_time).tm_min
        
        text = str(localtime_hour) + ':' + str(localtime_min)
        cv2.putText(img, text, (20,20), cv2.FONT_HERSHEY_PLAIN, 2, [0,255,0])
    
    def draw_mpl(self, ax, epoch_time):
        odor_status = self.get_odor_status(epoch_time)
        
        if odor_status[1] == 0:
            color_pos = 'cyan'
            color_neg = 'cyan'	
        else:
            if odor_status[2] > 0:
                color_pos = 'red'
                color_neg = 'cyan'
            elif odor_status[2] < 0:
                color_pos = 'cyan'
                color_neg = 'red'
        
        c1 = patches.Circle(( int(self.circular_region_pos['center'][0]), int(self.circular_region_pos['center'][1])), int(self.circular_region_pos['radius']), edgecolor=color_pos, facecolor='none', linewidth=2)
        c2 = patches.Circle(( int(self.circular_region_neg['center'][0]), int(self.circular_region_neg['center'][1])), int(self.circular_region_neg['radius']), edgecolor=color_neg, facecolor='none', linewidth=2)
        
        ax.add_artist(c1)
        ax.add_artist(c2)
        
        
        localtime_hour = time.localtime(epoch_time).tm_hour
        localtime_min = time.localtime(epoch_time).tm_min
        
        text_local = str(localtime_hour) + ':' + str(localtime_min)
        text_epoch = str(epoch_time)
        
        ax.text(20,20, text_epoch, color='green')
        
    def preprocess_data(self, pandas_dataframe):
        print 'Preprocessing data - see config file for details!'
        pandas_dataframe = mta.read_hdf5_file_to_pandas.remove_rows_above_speed_threshold(pandas_dataframe, speed_threshold=10)
        return pandas_dataframe
        
    def postprocess_data(self, pandas_dataframe):
        localtimefloat_t0 = time.localtime(self.data.time_epoch.min()).tm_hour + time.localtime(self.data.time_epoch.min()).tm_min/60. + time.localtime(self.data.time_epoch.min()).tm_sec/3600. + self.data.time_epoch_nsecs.min()*1e-9/3600.
        pandas_dataframe['time_local'] = (pandas_dataframe['time_epoch'] - pandas_dataframe['time_epoch'].min())/3600. + localtimefloat_t0
        
        for region_name, region_parameters in self.circular_regions.items(): 
            pandas_dataframe = mta.data_slicing.calc_frames_with_object_in_circular_region(   pandas_dataframe, 
                                                                                              region_parameters['center'], 
                                                                                              region_parameters['radius'], 
                                                                                              region_name=region_name)
        return pandas_dataframe                                                                    
                                                                      
                                                                  
                                                                  
