import os, sys
import imp
import numpy as np
import matplotlib.pyplot as plt
import pandas
import multi_tracker_analysis as mta
import figurefirst
import scipy.signal

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def get_confidence_intervals_for_circadian(data, iterations=1000, avg_func='mean'):
    reps = []
    if type(data) is list:
        data = np.array(data)

    for row in range(len(data)):
        nans, x= nan_helper(data[row])
        data[row][nans] = np.interp(x(nans), x(~nans), data[row][~nans])
        
        
    if avg_func == 'mean':
        np.avg_func = np.nanmean
    elif avg_func == 'median':
        np.avg_func = np.median
    
    for i in range(iterations):
        indices = np.random.randint(0,len(data),len(data))
        randoms = data[indices]
        m = np.avg_func(randoms, axis=0)
        reps.append(m)
    
    reps = np.array(reps)
    reps = np.sort(reps, axis=0)

    conf_hi = reps[int(iterations*(1-.025))]
    conf_lo = reps[int(iterations*(1-.975))]
    
    return conf_lo, conf_hi
    
def get_mean_speed_for_day(path):
    pd_filename = mta.read_hdf5_file_to_pandas.get_filename(path, 'pd_circadian.pickle')
    pd = pandas.read_pickle(pd_filename)
    times = pd.circadian_time[0] - pd.circadian_time[0].min()
    speed = pd.circadian_speed[0]
    
    if len(times) < 1200:
        time_padding = np.arange(times[-1], 20, 1/60.)[1:]
        zero_padding = np.zeros(1200-len(times))*np.nan
        times = np.hstack((times, time_padding))
        speed = np.hstack((speed, zero_padding))
    
    return times[0:1200], speed[0:1200]
                
def plot_circadian(ax, directory):
    
    paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='day')
    
    times = []
    speeds = []
    
    for path in paths:
        t, s = get_mean_speed_for_day(path)
        speeds.append(s)
        times.append(t)
    
    mean = np.nanmean(speeds, axis=0)
    mean = np.nan_to_num(mean)
    #std = np.std(speeds, axis=0)
    conf_lo, conf_hi = get_confidence_intervals_for_circadian(speeds)
    
    wn = 0.7
    b, a = scipy.signal.butter(3, wn)
    conf_lo = scipy.signal.filtfilt(b,a,conf_lo)
    conf_hi = scipy.signal.filtfilt(b,a,conf_hi)
    mean = scipy.signal.filtfilt(b,a,mean)
    
    
    t = np.nanmean(times, axis=0)
    
    ax.fill_between(t, conf_lo, conf_hi, edgecolor='none', facecolor=(0.001,0.001,0.001), alpha=0.3)
    ax.plot(t, mean, linewidth=1, color='black')
    
    start = 0
    sunset = 11-3
    sunrise = (24-15+7)
    end = 20
        
    ax.fill_between([start, sunset], [20,20], [0,0], facecolor='yellow', edgecolor='none', alpha=0.1, zorder=-100)
    ax.fill_between([sunset, sunrise], [20,20], [0,0], facecolor='gray', edgecolor='none', alpha=0.3, zorder=-100)
    ax.fill_between([sunrise, end], [20,20], [0,0], facecolor='yellow', edgecolor='none', alpha=0.1, zorder=-100)
    
    config = mta.read_hdf5_file_to_pandas.load_config_from_path(paths[0])
    if config.odor_control['action'][0] != 'off':
        print 'First message is not OFF!!'
        starttime = config.odor_control['t'][0]/3600. - 10*60/3600.
    else:
        starttime = config.odor_control['t'][0]/3600.
        
    if not config.legacy: 
        v = config.odor_control.value
    else:
        v = config.odor_control
    for r, row in enumerate(v):
        if row[5] != 'off':
            color = 'green'
            ax.fill_between([row[4]/3600.-starttime, row[4]/3600.-starttime+600/3600.], [20, 20], [0,0], facecolor='white', edgecolor='none', zorder=-50)
            print row[4]/3600., starttime
            ax.fill_between([row[4]/3600.-starttime, row[4]/3600.-starttime+600/3600.], [20, 20], [0,0], facecolor=color, edgecolor='none', zorder=-50, alpha=0.3)
    
    ax.set_ylim(config.ylim_speed[0], config.ylim_speed[1])
    ax.set_xlim(start, end)
    
    xticks = [0, sunset, sunrise, end]
    yticks = []
    #figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], xticks=xticks, yticks=yticks)
    figurefirst.mpl_functions.adjust_spines(ax, [], xticks=xticks)

def plot_circadian_for_tmaze_exps(ax, directory, t_range=[0,14]):

    paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='day')
    
    times = []
    speeds = []
    
    for path in paths:
        t, s = get_mean_speed_for_day(path)
        speeds.append(s)
        times.append(t)
    
    mean = np.nanmean(speeds, axis=0)
    mean = np.nan_to_num(mean)
    #std = np.std(speeds, axis=0)
    conf_lo, conf_hi = get_confidence_intervals_for_circadian(speeds)
    
    wn = 0.7
    b, a = scipy.signal.butter(3, wn)
    conf_lo = scipy.signal.filtfilt(b,a,conf_lo)
    conf_hi = scipy.signal.filtfilt(b,a,conf_hi)
    mean = scipy.signal.filtfilt(b,a,mean)
    
    
    t = np.nanmean(times, axis=0)
    
    ax.fill_between(t, conf_lo, conf_hi, edgecolor='none', facecolor=(0.001,0.001,0.001), alpha=0.3)
    ax.plot(t, mean, linewidth=1, color='black')


    ax.fill_between([np.min(t), np.max(t)], [20,20], [0,0], facecolor='yellow', edgecolor='none', alpha=0.1, zorder=-100)



    config = mta.read_hdf5_file_to_pandas.load_config_from_path(paths[0])
    if config.odor_control['action'][0] != 'off':
        print 'First message is not OFF!!'
        starttime = config.odor_control['t'][0]/3600. - 10*60/3600.
    else:
        starttime = config.odor_control['t'][0]/3600.
        
    if not config.legacy: 
        v = config.odor_control.value
    else:
        v = config.odor_control
    for r, row in enumerate(v):
        if row[5] != 'off':
            color = 'green'
            ax.fill_between([row[4]/3600.-starttime, row[4]/3600.-starttime+600/3600.], [20, 20], [0,0], facecolor='white', edgecolor='none', zorder=-50)
            print row[4]/3600., starttime
            ax.fill_between([row[4]/3600.-starttime, row[4]/3600.-starttime+600/3600.], [20, 20], [0,0], facecolor=color, edgecolor='none', zorder=-50, alpha=0.3)
    
    ax.set_xlim(t_range[0]/60., t_range[1]/60.)
    ax.set_ylim(config.ylim_speed[0], config.ylim_speed[1])