import os, sys
import imp
import numpy as np
import matplotlib.pyplot as plt
import pandas
import figurefirst
import multi_tracker_analysis as mta
import matplotlib
import fly_plot_lib.plot as fpl
import flystat
import scipy.signal
import co2_paper_locations

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
                
def get_one_bootstrapped_mean(data):
    indices = np.random.randint(0, len(data), len(data))
    selection = [data[i] for i in indices]
    return np.mean(selection, axis=0)    

def get_95_confidence_intervals(data, iterations=100):
    bootstrapped_data = np.mean(data, axis=0)
    for iteration in range(iterations):
        bootstrapped_data = np.vstack( (bootstrapped_data, get_one_bootstrapped_mean(data)) )
    bootstrapped_data.sort(axis=0)
    
    index_hi = int(iterations*0.975)
    index_lo = int(iterations*0.025)
    
    return bootstrapped_data[index_lo, :], bootstrapped_data[index_hi, :]

def plot_mean_and_95_confidence(ax, x, y, color, iterations=100, normalize_to_one=False, show_confidence=True, color_mean=False, zorder=1):
    if type(ax) is list:
        layout, figure, axis = ax
        ax = layout.axes[(figure, axis)]
        fifidatafile = layout.fifidatafile
    else:
        fifidatafile = ''
    wn = 0.05
    
    if normalize_to_one:
        mean_y = np.mean(y, axis=0)
        b, a = scipy.signal.butter(3, wn)
        mean_y = scipy.signal.filtfilt(b,a,mean_y)
        third = int(0.25*len(mean_y))
        baseline_shift = np.mean( mean_y[0:third] )
        mean_y -= baseline_shift
        norm_factor = np.max(mean_y)
    else:
        norm_factor = 1
        baseline_shift = 0

    if len(y) > 1:
        if show_confidence:
            lo, hi = get_95_confidence_intervals(y, iterations=iterations)
            
            b, a = scipy.signal.butter(3, wn)
            lo_filtered = scipy.signal.filtfilt(b,a,lo)
            hi_filtered = scipy.signal.filtfilt(b,a,hi)
            
            #ax.fill_between(np.mean(x, axis=0), (lo_filtered-baseline_shift)/norm_factor, (hi_filtered-baseline_shift)/norm_factor, facecolor=color, alpha=0.2, edgecolor='none', zorder=zorder-1)
            y_baseline_shifted = [(yi-baseline_shift)/norm_factor for yi in y]
            #figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_line',
            #                     ax[0], ax[1], ax[2], fifidatafile, 'data confidence interval ' + color, 
            #                    ['Time',
            #                     'List of n flies'],
            #                     np.mean(x, axis=0), y_baseline_shifted,
            #                     color=color, shading='95conf', show_lines=False, show_mean=False, alpha=0.2)
            if ax.breadcrumb['layout_filename'] == co2_paper_locations.figure_template_locations.figure4_walking_arena_activity:
                title = 'Mean number of flies in '+color+' zone'
                description2 = 'List of mean number of flies for each trial'
            else:
                title = 'Preference Index for flies'
                description2 = 'List of PI for each trial'
            ax._custom([title, 'Time (secs)', description2], 'fly_plot_lib.plot.scatter_line', np.mean(x, axis=0), y_baseline_shifted, 
                        color=color, shading='95conf', show_lines=False, show_mean=False, alpha=0.2)

    mean = (np.mean(y, axis=0)-baseline_shift)/norm_factor
    
    print '********************: ', np.max(mean), norm_factor 
    b, a = scipy.signal.butter(3, wn)
    mean_filtered = scipy.signal.filtfilt(b,a,mean)
    
    if color_mean is False:
        #ax.plot(np.mean(x, axis=0), mean_filtered, color=color, zorder=zorder)
        #figurefirst.deprecated_regenerate.mpl('plot', ax[0], ax[1], ax[2], fifidatafile, 'data mean ' + color, 
        #                      ['Time',
        #                       '2.75 percent Confidence',
        #                       '97.5 percent Confidence'], 
        #                       np.mean(x, axis=0), mean_filtered, color=color, zorder=zorder)
        if ax.breadcrumb['layout_filename'] == co2_paper_locations.figure_template_locations.figure4_walking_arena_activity:
            title = 'Mean number of flies in '+color+' zone'
            description2 = 'Mean number of flies'
        else:
            title = 'Mean Preference Index for flies'
            description2 = 'Mean PI flies'
        ax._plot([title, 'Time (secs)', description2], np.mean(x, axis=0), mean_filtered, color=color, zorder=zorder)

    else:
        indices = np.arange(0,len(mean_filtered),10)
        x = np.mean(x, axis=0)[indices]
        y = mean_filtered[indices]
        c = color_mean['color'][indices]
        ax._custom(['Colored Preference Index for flies', 
                    'Time (secs)', 'Mean PI flies', 'Color (Mean PI flies)'], 'fly_plot_lib.plot.colorline', 
                    x, 
                    y, 
                    c,
                    cmap=color_mean['cmap'],
                    norm=color_mean['norm'],
                    zorder=zorder)

def plot_n_flies_in_control_and_odor_regions(ax, paths, localtimerange, flowrate, average_within_paths=False, normalize_to_one=False, traces_to_plot=['odor', 'control'], 
                                             show_confidence=True, color_mean=False, colors={'odor': 'red', 'control':'blue'}, fill_odor=True):
    layout, figure, axis = ax
    fifidatafile = layout.fifidatafile
    if type(paths) is not list:
        print 'LOADING PATHS'
        paths = mta.read_hdf5_file_to_pandas.get_filenames(paths, contains='day')
        
    control = np.array([])
    odor = np.array([])
    speed = np.array([])
    t = np.array([])
    
    for p, path in enumerate(paths):
        pd_pickle_filename = get_filename(path, 'pd_data.pickle', does_not_contain=['~', '.pyc'])
        pd = pandas.read_pickle(pd_pickle_filename)
        config = mta.read_hdf5_file_to_pandas.load_config_from_path(path)
        if type(localtimerange) is str:
            localtimerange = config.portion_of_day_to_local_time_range[localtimerange]
        
        query = "(action == 'left' or action == 'right') and (flowrate < " + str(flowrate+.005) + " and flowrate > " + str(flowrate-0.005) + ") and (localtime < " + str(localtimerange[1]) + " and localtime > " + str(localtimerange[0]) + ")"
        pd_tmp = pd.query(query)
                
        if average_within_paths:
            raise NotImplementedError
        
        control = np.hstack((control, pd_tmp.n_flies_control.values/float(config.nflies)))
        odor = np.hstack((odor, pd_tmp.n_flies_odor.values/float(config.nflies)))
        speed = np.hstack((speed, pd_tmp.speed.values))
        t = np.hstack((t, pd_tmp.t.values))
    
    config = mta.read_hdf5_file_to_pandas.load_config_from_path(paths[-1])
    
    if fill_odor:
        #ax.fill_between([0, 600], 0, 10, facecolor='green', edgecolor='none', alpha=0.2)
        figurefirst.deprecated_regenerate.mpl('fill_between', layout, figure, axis, fifidatafile, 'Odor stimulus', 
                              ['Time odor is on (secs)'], 
                              [0, 600], 0, 10, facecolor='green', edgecolor='none', alpha=0.2)
    
    if 'odor' in traces_to_plot:
        plot_mean_and_95_confidence(ax, t, odor, colors['odor'], iterations=100, normalize_to_one=normalize_to_one, show_confidence=show_confidence, color_mean=color_mean)
    if 'control' in traces_to_plot:
        plot_mean_and_95_confidence(ax, t, control, colors['control'], iterations=100, normalize_to_one=normalize_to_one, show_confidence=show_confidence)
    if 'speed' in traces_to_plot:
        plot_mean_and_95_confidence(ax, t, speed, (0.001, 0.001, 0.001), iterations=100, normalize_to_one=normalize_to_one, show_confidence=show_confidence)
    
    if 'speed' in traces_to_plot:
        #ax.set_ylim(0,7)
        figurefirst.deprecated_regenerate.mpl('set_ylim', ax[0], ax[1], ax[2], fifidatafile, 'ylim', [], 
                              0,7)
    else:
        #ax.set_ylim(0,1)
        figurefirst.deprecated_regenerate.mpl('set_ylim', ax[0], ax[1], ax[2], fifidatafile, 'ylim', [], 
                              0,1)

    #ax.set_xlim(config.xlim_nflies[0], config.xlim_nflies[1]) #-600, 1200)
    figurefirst.deprecated_regenerate.mpl('set_xlim', ax[0], ax[1], ax[2], fifidatafile, 'xlim', [], 
                              config.xlim_nflies[0], config.xlim_nflies[1]) #-600, 1200)
    
    try:
        #figurefirst.mpl_functions.adjust_spines(ax, [])
        figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', ax[0], ax[1], ax[2], fifidatafile, 'adjust spines', [], 
                              [])

    except:
        pass # probably already did this

def get_nflies_and_speed_at_time_for_directory(directory, localtimerange, flowrate, time_secs, smooth=False, combine_cohorts=False):
    '''
    time_secs - can be a number (seconds relative to odor stimulus), or 'min' or 'max' to get peak/min. Or, if a list of two numbers, function will return mean in that time range.
    
    '''
    
    paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='day')
    print paths
    
    n_flies_odor = []
    n_flies_control = []
    speed = []
    for path in paths:
        pd_pickle_filename = get_filename(path, 'pd_data.pickle', does_not_contain=['~', '.pyc'])
        pd = pandas.read_pickle(pd_pickle_filename)
        config = mta.read_hdf5_file_to_pandas.load_config_from_path(path)
        if type(localtimerange) is str:
            localtimerange_value = config.portion_of_day_to_local_time_range[localtimerange]
        else:
            localtimerange_value = localtimerange

        query = "(action == 'left' or action == 'right') and (flowrate < " + str(flowrate+.005) + " and flowrate > " + str(flowrate-0.005) + ") and (localtime < " + str(localtimerange_value[1]) + " and localtime > " + str(localtimerange_value[0]) + ")"
        pd_tmp = pd.query(query)
        
        print path
        print '        SHAPE: ', pd_tmp.shape
        
        n_flies_odor_for_exp = []
        n_flies_control_for_exp = []
        speed_for_exp = []
        for index in pd_tmp.index:
            
            
            if type(time_secs) is list:
                time_secs_loc_0 = np.abs(pd_tmp.loc[index].t - time_secs[0]).argmin()
                time_secs_loc_1 = np.abs(pd_tmp.loc[index].t - time_secs[-1]).argmin()
                #time_range_length = float(len(pd_tmp.loc[index].n_flies_odor[time_secs_loc_0:time_secs_loc_1]))
                
                _n_flies_odor = pd_tmp.loc[index].n_flies_odor
                _n_flies_control = pd_tmp.loc[index].n_flies_control
                _speed = pd_tmp.loc[index].speed
                if smooth:
                    wn = 0.01
                    b, a = scipy.signal.butter(3, wn)
                    _n_flies_odor = scipy.signal.filtfilt(b,a,_n_flies_odor)
                    _n_flies_control = scipy.signal.filtfilt(b,a,_n_flies_control)
                    _speed = scipy.signal.filtfilt(b,a,_speed)
                    
                    if 0:
                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.plot(pd_tmp.loc[index].n_flies_odor)
                        ax.plot(_n_flies_odor)
                
                n_flies_odor_for_exp.append( np.mean(_n_flies_odor[time_secs_loc_0:time_secs_loc_1])/float(config.nflies) )
                n_flies_control_for_exp.append( np.mean(_n_flies_control[time_secs_loc_0:time_secs_loc_1])/float(config.nflies) )
                speed_for_exp.append( np.mean(_speed[time_secs_loc_0:time_secs_loc_1]) )
                print speed_for_exp
                
                
            else:
                if time_secs == 'max':
                    time_secs_loc = pd_tmp.loc[index].n_flies_odor.argmax()
                elif time_secs == 'min':
                    time_secs_loc = pd_tmp.loc[index].n_flies_odor.argmin()
                else:
                    time_secs_loc = np.abs(pd_tmp.loc[index].t - time_secs).argmin()
                
                _n_flies_odor = pd_tmp.loc[index].n_flies_odor/float(config.nflies)
                _n_flies_control = pd_tmp.loc[index].n_flies_control/float(config.nflies)
                _speed = pd_tmp.loc[index].speed
                if smooth:
                    wn = 0.01
                    b, a = scipy.signal.butter(3, wn)
                    _n_flies_odor = scipy.signal.filtfilt(b,a,_n_flies_odor)
                    _n_flies_control = scipy.signal.filtfilt(b,a,_n_flies_control)
                    _speed = scipy.signal.filtfilt(b,a,_speed)
                    
                    if 0:
                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        ax.plot(pd_tmp.loc[index].n_flies_odor)
                        ax.plot(_n_flies_odor)
                
                n_flies_odor_for_exp.append( _n_flies_odor[time_secs_loc] )
                n_flies_control_for_exp.append( _n_flies_control[time_secs_loc] )
                speed_for_exp.append( _speed[time_secs_loc] )
            
        if combine_cohorts:
            n_flies_odor.extend(n_flies_odor_for_exp)
            n_flies_control.extend(n_flies_control_for_exp)
            speed.extend(speed_for_exp)
        else:
            n_flies_odor.append(n_flies_odor_for_exp)
            n_flies_control.append(n_flies_control_for_exp)
            speed.append(speed_for_exp)

    return {'n_flies_odor': n_flies_odor, 'n_flies_control': n_flies_control, 'speed': speed, 'config': config}
    
def plot_raw_scatter_with_cmap(ax, scatter_data, key, average_within_paths=False):
    if key == 'n_flies_odor':
        cmap = 'cool'
        confidence_color = 'red'
    if key == 'n_flies_control':
        cmap = 'cool'
        confidence_color = 'blue'
    if key == 'speed':
        cmap = 'summer'
        confidence_color = (0.0001, 0.0001, 0.0001)
        
    cmap = matplotlib.cm.get_cmap(cmap)
    startx = 1
    for d, y_data in enumerate(scatter_data[key]):
        color = cmap( d/float(len(scatter_data[key])) )
        print key
        print y_data
        fpl.scatter_box(ax, startx+d, y_data, xwidth=0.5, ywidth=0.1, color='black', marker_linewidth=0, flipxy=False, shading='none', markersize=1.6, random_scatter=False)
        fpl.scatter_box(ax, startx+d, y_data, xwidth=0.5, ywidth=0.1, color=color, marker_linewidth=0, flipxy=False, shading='none', markersize=1.4, random_scatter=False)
    
    if average_within_paths:
        all_data_points = np.array([np.nanmean(d) for d in scatter_data[key]])
    else:
        all_data_points = np.hstack(scatter_data[key])
    
    conf_interval = flystat.resampling.bootstrap_confidence_intervals_from_data(all_data_points, iterations=5000, use='mean')
    ax.hlines([np.nanmean(all_data_points)], startx-1.5, len(scatter_data[key])+1.5, colors=[confidence_color], linewidth=1)
    ax.fill_between([startx-1.5, len(scatter_data[key])+1.5], [conf_interval[0], conf_interval[0]], [conf_interval[1], conf_interval[1]], facecolor=confidence_color, edgecolor='none', alpha=0.2)
    
    

    if key == 'n_flies_odor' or key == 'n_flies_control':
        n = 1#scatter_data['config'].nflies
        ax.set_ylim(0, n)
    elif key == 'speed':
        ax.set_ylim(0, 7)
        
    figurefirst.mpl_functions.adjust_spines(ax, [])

def plot_raw_scatter(ax, scatter_data, key, average_within_paths=False, hide_markers=False):
    layout, figure, axis = ax
    fifidatafile = layout.fifidatafile
    if key == 'n_flies_odor':
        color = 'red'
        confidence_color = 'red'
    if key == 'n_flies_control':
        color = 'blue'
        confidence_color = 'blue'
    if key == 'speed':
        color = (0.0001, 0.0001, 0.0001)
        confidence_color = (0.0001, 0.0001, 0.0001)
        
    startx = 1
    for d, y_data in enumerate(scatter_data[key]):
        #color = cmap( d/float(len(scatter_data[key])) )
        print key
        print y_data
        #fpl.scatter_box(ax, startx+d, y_data, xwidth=0.5, ywidth=0.1, color='black', marker_linewidth=0, flipxy=False, shading='none', markersize=1, random_scatter=False, hide_markers=hide_markers)
        #fpl.scatter_box(ax, startx+d, y_data, xwidth=0.5, ywidth=0.1, color=color, marker_linewidth=0, shading='95conf', flipxy=False, markersize=0.8, random_scatter=False, hide_markers=hide_markers)
        #figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', ax[0], ax[1], ax[2], fifidatafile, '95 percent CI for ' + key,
        #                          ['index',
        #                           'List of trajectories that entered control and binary value for whether they entered test volume'],
        #                           startx+d, y_data, xwidth=0.5, ywidth=0.1, color=color, marker_linewidth=0, shading='95conf', flipxy=False, markersize=0.8, random_scatter=False, hide_markers=hide_markers)



    if average_within_paths:
        all_data_points = np.array([np.nanmean(d) for d in scatter_data[key]])
    else:
        all_data_points = np.hstack(scatter_data[key])
    
    conf_interval = flystat.resampling.bootstrap_confidence_intervals_from_data(all_data_points, iterations=5000, use='mean')
    

    if 0:
        #ax.hlines([np.nanmean(all_data_points)], startx-1.5, len(scatter_data[key])+1.5, colors=[confidence_color], linewidth=1)
        figurefirst.deprecated_regenerate.mpl('hlines', ax[0], ax[1], ax[2], fifidatafile, 'Mean '+key,
                                  ['Mean '+key],
                                   [np.nanmean(all_data_points)], startx-1.5, len(scatter_data[key])+1.5, colors=[confidence_color], linewidth=1)
        #ax.fill_between([startx-1.5, len(scatter_data[key])+1.5], [conf_interval[0], conf_interval[0]], [conf_interval[1], conf_interval[1]], facecolor=confidence_color, edgecolor='none', alpha=0.2)
        figurefirst.deprecated_regenerate.mpl('fill_between', ax[0], ax[1], ax[2], fifidatafile, '95 percent confidence '+key,
                                  ['95 perc CI '+key],
                                  [startx-1.5, len(scatter_data[key])+1.5], [conf_interval[0], conf_interval[0]], [conf_interval[1], conf_interval[1]], facecolor=confidence_color, edgecolor='none', alpha=0.2)
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_box',
                                 ax[0], ax[1], ax[2], fifidatafile, 'Mean '+key,
                                 ['index',
                                  'Values for '+key+' for each trial'],
                                 0, all_data_points, 
                                 xwidth=3, ywidth=0.1, color=confidence_color, edgecolor='none', flipxy=False, shading='95conf', linewidth=1, use='mean', hide_markers=True, alpha=0.2)
        

    if key == 'n_flies_odor' or key == 'n_flies_control':
        n = 1#scatter_data['config'].nflies
        #ax.set_ylim(0, n)
        figurefirst.deprecated_regenerate.mpl('set_ylim', ax[0], ax[1], ax[2], fifidatafile, 'set_ylim', [], 
                              0, n)
    elif key == 'speed':
        #ax.set_ylim(0, 7)
        figurefirst.deprecated_regenerate.mpl('set_ylim', ax[0], ax[1], ax[2], fifidatafile, 'set_ylim', [], 
                              0, 7)
        
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', ax[0], ax[1], ax[2], fifidatafile, 'adjust spines', [], 
                              [])
        
