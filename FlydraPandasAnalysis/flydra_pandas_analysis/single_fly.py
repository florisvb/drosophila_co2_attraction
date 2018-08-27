import multi_tracker_analysis as mta
import numpy as np

import fly_plot_lib as fpl

def get_time_local(t):
    tl = time.localtime(t)
    lt = tl.tm_hour + tl.tm_min/60. + tl.tm_sec/3600.
    return lt

def calc_localtime(pd):
    lt = [get_time_local(t) for t in pd.time_epoch.values]
    pd['time_local'] = lt
    
def calc_relative_time(pd):
    lt = get_time_local( np.min(pd.time_epoch.values) )
    dt = 23 - lt
    relative_timepoint = np.min(pd.time_epoch.values) + dt
    pd['time_relative'] = pd.time_epoch - relative_timepoint
    

def get_pd_for_odor_on(pd):
    if 'time_relative' not in pd.keys():
        calc_relative_time(pd)
    
    pd_odor = pd.query('time_relative > 0 & time_relative < 14400')
    
    return pd_odor
    
def calc_pad_blob(pd):
    x_range = [-.055, .055]
    y_range = [-.08, .08]
    z_range = [-.02, .007]
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(pd, x_range, y_range, z_range, region_name='onpad')
    
    
    x_range = [.316, .404]
    y_range = [-.1,.1]
    z_range = [-.15,-.089]
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(pd, x_range, y_range, z_range, region_name='blob')
    
def plot_pad_dot(pd):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    pd_q = pd.query('onpad > 0')
    chunks, breaks = fpl.flymath.get_continuous_chunks(pd_q.time_relative.values, jump=1, return_index=False)
    for chunk in chunks:
        ax.fill_betweenx([0,1], chunk[0], chunk[-1], facecolor='red', edgecolor='none')
        
    
    pd_q = pd.query('blob > 0')
    chunks, breaks = fpl.flymath.get_continuous_chunks(pd_q.time_relative.values, jump=1, return_index=False)
    for chunk in chunks:
        ax.fill_betweenx([0,1], chunk[0], chunk[-1], facecolor='black', edgecolor='none')
        
    pd_q = pd.query('blob == 0 & onpad == 0')
    chunks, breaks = fpl.flymath.get_continuous_chunks(pd_q.time_relative.values, jump=1, return_index=False)
    for chunk in chunks:
        ax.fill_betweenx([0,1], chunk[0], chunk[-1], facecolor='cyan', edgecolor='none')
        
        
        
        
        
        
        
        
        
        
        
        
