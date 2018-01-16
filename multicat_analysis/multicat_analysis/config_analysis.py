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
            
def load_bag_as_hdf5(bag, skip_messages={}, topics=None):
    output_fname = bag.split('.')[0] + '.hdf5'
    print output_fname
    if not os.path.exists(output_fname):
        mta.bag2hdf5.bag2hdf5(   bag,
                                 output_fname,
                                 topics=topics,
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
        self.legacy = False
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


        self.actual_circular_region_pos = {'center': [585,240], 'radius': 40}
        self.actual_circular_region_neg = {'center': [80,242], 'radius': 40}
        self.actual_circular_regions = {'pos': self.circular_region_pos,
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
        
        if '_N1' in self.identifiercode:
            self.flip = False
        elif '_N2' in self.identifiercode:
            self.flip = True                                      

    def set_topics(self):
        self.topics = self.topics = ['/phidgets_interface_ssr', '/alicat_bb9']
        if '_N1' in self.identifiercode:
            rig = '/rig1'
        elif '_N2' in self.identifiercode:
            rig = '/rig2'  
        self.topics.append(rig)
        
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
        odor_control_bag_filename_contains = self.identifiercode + '_alicat_flow_control.bag'
        odor_control_bag_filename = get_filename(self.path, odor_control_bag_filename_contains)
        if 'skip_odor_control_messages' not in self.__dict__.keys():
            self.skip_odor_control_messages = {}
        if 'topics' not in self.__dict__.keys():
            self.set_topics()
            
        self.odor_bag_data = load_bag_as_hdf5(odor_control_bag_filename, skip_messages=self.skip_odor_control_messages, topics=self.topics)
        if '_N1' in self.identifiercode:
            rig = 'rig1'
        elif '_N2' in self.identifiercode:
            rig = 'rig2'  
        try:
            self.odor_control = self.odor_bag_data[rig]
        except:
            self.odor_control = self.odor_bag_data['rig1'] # for windy orchard two pad choice case
            
        # delta video player has problems with this ^ format for some reason. The following helps
        self.odor_control_time = np.array(self.odor_control['header_stamp'])
        self.odor_control_action = np.array(self.odor_control['action'])
        self.odor_control_flowrate = np.array(self.odor_control['flowrate'])

        self.sensory_stimulus_on = []
        self.sensory_stimulus_rgba = []
        for index, t in enumerate(self.odor_control_time):
            if index < len(self.odor_control_time)-1: # to deal with broken ends
                if self.odor_control_action[index] != 'off':
                    self.sensory_stimulus_on.append([t, self.odor_control_time[index+1]])
                    if self.odor_control_flowrate[index] != 0:
                        if not self.flip:
                            if self.odor_control_action[index] == 'left':
                                self.sensory_stimulus_rgba.append((255,0,0,150))
                            if self.odor_control_action[index] == 'right':
                                self.sensory_stimulus_rgba.append((0,0,255,150))
                        else:
                            if self.odor_control_action[index] == 'left':
                                self.sensory_stimulus_rgba.append((0,0,255,150))
                            if self.odor_control_action[index] == 'right':
                                self.sensory_stimulus_rgba.append((255,0,0,150))
                    else:
                        self.sensory_stimulus_rgba.append((125,125,125,150))
        
    def process_alicat_data(self):
        alicat_bag_filename_contains = self.identifiercode + '_alicat_flow_data.bag'
        alicat_bag_filename = get_filename(self.path, alicat_bag_filename_contains)
        if os.stat(alicat_bag_filename).st_size < 5000:
            raise ValueError('Check Alicat Data! Might be Erroneous!!')

        if 'skip_alicat_messages' not in self.__dict__.keys():
            self.skip_alicat_messages = {}

        self.alicat_data = load_bag_as_hdf5(alicat_bag_filename, skip_messages=self.skip_alicat_messages)
    
    def get_odor_status(self, epoch_time):
        if epoch_time < self.odor_control_time[0]:
            return self.odor_control_action[0]
        d = epoch_time - self.odor_control_time
        index = np.where(d>0)[0][-1]
        return self.odor_control_action[index]
        
    def draw(self, img, epoch_time, color_order='rgb'):
        odor_status = self.get_odor_status(epoch_time)
        
        if odor_status == 'off':
            #color_pos = [0,0,255]
            color_pos = [255,0,0]
            #color_neg = [0,0,255]
            color_neg = [255,0,0]
        else:
            if not self.flip:
                if odor_status == 'right':
                    color_pos = [0,0,255]#[255,0,0]
                    color_neg = [255,0,0]#[0,0,255]
                elif odor_status == 'left':
                    color_pos = [255,0,0]#[0,0,255]
                    color_neg = [0,0,255]#[255,0,0]
            else:
                if odor_status == 'right':
                    color_pos = [255,0,0]#[0,0,255]
                    color_neg = [0,0,255]#[255,0,0]
                elif odor_status == 'left':
                    color_pos = [0,0,255]#[255,0,0]
                    color_neg = [255,0,0]#[0,0,255]

        cv2.circle(img, ( int(self.circular_region_pos['center'][0]), int(self.circular_region_pos['center'][1])), int(self.circular_region_pos['radius']), color_pos, 2)
        cv2.circle(img, ( int(self.circular_region_neg['center'][0]), int(self.circular_region_neg['center'][1])), int(self.circular_region_neg['radius']), color_neg, 2)
        
        localtime_hour = time.localtime(epoch_time).tm_hour
        localtime_min = time.localtime(epoch_time).tm_min
        
        text = str(localtime_hour) + ':' + str(localtime_min)
        cv2.putText(img, text, (20,20), cv2.FONT_HERSHEY_PLAIN, 2, [0,255,0])
    
    def draw_mpl(self, ax, epoch_time, actual_region=True):
        print epoch_time
        odor_status = self.get_odor_status(epoch_time)

        print odor_status
        
        off_color = 'blue'

        if odor_status == 'off':
            color_pos = off_color
            color_neg = off_color
        else:
            if not self.flip:
                if odor_status == 'right':
                    color_pos = 'red'
                    color_neg = off_color
                elif odor_status == 'left':
                    color_pos = off_color
                    color_neg = 'red'
            else:
                if odor_status == 'right':
                    color_pos = off_color
                    color_neg = 'red'
                elif odor_status == 'left':
                    color_pos = 'red'
                    color_neg = off_color
        
        if actual_region is False:
            c1 = patches.Circle(( int(self.circular_region_pos['center'][0]), int(self.circular_region_pos['center'][1])), int(self.circular_region_pos['radius']), edgecolor=color_pos, facecolor='none', linewidth=2)
            c2 = patches.Circle(( int(self.circular_region_neg['center'][0]), int(self.circular_region_neg['center'][1])), int(self.circular_region_neg['radius']), edgecolor=color_neg, facecolor='none', linewidth=2)
        else:
            c1 = patches.Circle(( int(self.actual_circular_region_pos['center'][0]), int(self.actual_circular_region_pos['center'][1])), int(self.actual_circular_region_pos['radius']), edgecolor='none', facecolor=color_pos, linewidth=2, alpha=0.25)
            c2 = patches.Circle(( int(self.actual_circular_region_neg['center'][0]), int(self.actual_circular_region_neg['center'][1])), int(self.actual_circular_region_neg['radius']), edgecolor='none', facecolor=color_neg, linewidth=2, alpha=0.25)

        ax.add_artist(c1)
        ax.add_artist(c2)
        
        localtime_hour = time.localtime(epoch_time).tm_hour
        localtime_min = time.localtime(epoch_time).tm_min
        
        #text_local = str(localtime_hour) + ':' + str(localtime_min)
        #text_epoch = str(epoch_time)
        #ax.text(20,20, text_epoch, color='green')
        
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
                                                                      
                                                                  
    def get_local_time_range_around_nth_odor_pulse(self, nth=1):
        # nth=1 is first pulse
        flowrates = np.unique(self.odor_control.value['flowrate'])
        actions = np.unique(self.odor_control.value['action'])

        keyframe_times = []
        for action in actions:
            if action == 'off':
                continue
            for flowrate in flowrates:
                if flowrate == 0:
                    continue
                keyframe_times.extend(self.get_keyframe_times(action, flowrate))

        keyframe_times.sort()
        if self.odor_control.value['flowrate'][0]==0:
            t0 = self.odor_control.value['t'][0]
        else:
            t0 = self.alicat_data.values()[0].value['t'][0]
        t = keyframe_times[nth-1] - t0

        localtimefloat_t0 = time.localtime(t0).tm_hour + time.localtime(t0).tm_min/60. + time.localtime(t0).tm_sec/3600.
        localtimefloat_t1 = localtimefloat_t0 + t/3600. + 30/60.
        return [localtimefloat_t0, localtimefloat_t1]
