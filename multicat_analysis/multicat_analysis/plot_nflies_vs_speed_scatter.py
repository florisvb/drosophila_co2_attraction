import plot_n_flies
import multi_tracker_analysis as mta

import os, sys
import imp
import numpy as np
import matplotlib.pyplot as plt
import pandas
import figurefirst
import matplotlib
import fly_plot_lib.plot as fpl
import fly_plot_lib.text as flytext
import flystat
import scipy.signal
from optparse import OptionParser



def get_directories(path, contains):
    return mta.read_hdf5_file_to_pandas.get_filenames(path, contains, does_not_contain=['~', '.pyc', 'CTRL_AIR'])

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
    
def get_attraction_index_and_speed(directory, flowrate, time_of_day):
    
    scatter_data_black = plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, time_of_day, flowrate, [-10*60, -2*60], smooth=False, combine_cohorts=True)
    #AI_black = np.array([np.array(s) for s in scatter_data_black['n_flies_odor']])
    AI_black = np.array([np.array(s) for s in scatter_data_black['n_flies_odor']]) - np.array([np.array(s) for s in scatter_data_black['n_flies_control']])
    speed_black = np.array([np.array(s) for s in scatter_data_black['speed']])
    
    scatter_data_magenta = plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, time_of_day, flowrate, [1*60, 3*60], smooth=True, combine_cohorts=True)
    AI_magenta = np.array([np.array(s) for s in scatter_data_magenta['n_flies_odor']]) - np.array([np.array(s) for s in scatter_data_magenta['n_flies_control']])
    
    scatter_data_orange = plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, time_of_day, flowrate, [8*60, 10*60], smooth=True, combine_cohorts=True)
    AI_orange = np.array([np.array(s) for s in scatter_data_orange['n_flies_odor']]) - np.array([np.array(s) for s in scatter_data_orange['n_flies_control']])
    
    
    #return speed_black, AI_black
    
    #odor = np.array([np.array(s) for s in scatter_data_magenta['n_flies_odor']]) - np.array([np.array(s) for s in scatter_data_black['n_flies_odor']])
    #control = np.array([np.array(s) for s in scatter_data_magenta['n_flies_control']]) - np.array([np.array(s) for s in scatter_data_black['n_flies_control']])
    AI = AI_magenta + AI_orange
    
    #AI = AI_magenta #- AI_black# + AI_orange - AI_black
    
    return speed_black, AI    

def get_new_attraction_index(directories, flowrates=1, localtimerange=[14,36], average_within_paths=False, use_first_half_of_odor_presentation=False):
    control = np.array([])
    odor = np.array([])
    speed = np.array([])
    t = np.array([])
    experiment_type = []

    for p, directory in enumerate(directories):
        paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='day')
        print 'Loading:  '
        print paths
        if len(paths) == 0:
            continue
        for path in paths:

            if type(flowrates) is list: # if doing different flow rates for each directory
                flowrate = flowrates[p]
            else:
                flowrate = flowrates

            pd_pickle_filename = get_filename(path, 'pd_data.pickle', does_not_contain=['~', '.pyc'])
            print pd_pickle_filename
            pd = pandas.read_pickle(pd_pickle_filename)
            config = mta.read_hdf5_file_to_pandas.load_config_from_path(path)
            if type(localtimerange) is str:
                use_localtimerange = config.portion_of_day_to_local_time_range[localtimerange]
            else:
                use_localtimerange = localtimerange
            if 'hot' in directory and localtimerange==[14,36]: # if hot directory, and localtime range is default
                use_localtimerange = [14, 30]
            
            query = "(action == 'left' or action == 'right') and (flowrate < " + str(flowrate+.005) + " and flowrate > " + str(flowrate-0.005) + ") and (localtime <= " + str(use_localtimerange[1]) + " and localtime >= " + str(use_localtimerange[0]) + ")"
            pd_tmp = pd.query(query)


            print
            print
            print path
            print localtimerange
            print query
            print pd_tmp
            if len(pd_tmp) == 0:
                print 'NO DATA HERE'
                continue

            experiment_types = {'hot': 2, 'noonstarved': 0, '24hrstarved': 1}
            N = pd_tmp.n_flies_control.values.shape[0]
            for et, val in experiment_types.items():
                if et in directory:
                    experiment_type.extend([val]*N)
                

            if average_within_paths:
                pd_tmp = pd_tmp[['t', 'speed', 'n_flies_control', 'n_flies_odor']].sum(axis=0)/float(pd_tmp.shape[0])
                tmp_n_flies_control = np.reshape(pd_tmp.n_flies_control, [1, len(pd_tmp.n_flies_control)])
                tmp_n_flies_odor = np.reshape(pd_tmp.n_flies_odor, [1, len(pd_tmp.n_flies_control)])
                tmp_speed = np.reshape(pd_tmp.speed, [1, len(pd_tmp.n_flies_control)])
                tmp_t = np.reshape(pd_tmp.t, [1, len(pd_tmp.n_flies_control)])

                print 'averaging within paths!!!!!!'
                print type(tmp_n_flies_control)
                print tmp_n_flies_control.shape

            else:
                pd_tmp = pd_tmp[['t', 'speed', 'n_flies_control', 'n_flies_odor']]
                tmp_n_flies_control = np.vstack(pd_tmp.n_flies_control.values)
                tmp_n_flies_odor = np.vstack(pd_tmp.n_flies_odor.values)
                tmp_speed = np.vstack(pd_tmp.speed.values)
                tmp_t = np.vstack(pd_tmp.t.values)

            try:
                if len(control) == 0:
                    control = np.vstack(tmp_n_flies_control)
                else:
                    control = np.vstack((control, tmp_n_flies_control))

                print control

                if len(odor) == 0:
                    odor = np.vstack(tmp_n_flies_odor)
                else:
                    odor = np.vstack((odor, tmp_n_flies_odor))

                if len(speed) == 0:
                    speed = np.vstack(tmp_speed)
                else:
                    speed = np.vstack((speed, tmp_speed))

                if len(t) == 0:
                    t = np.vstack(tmp_t)
                else:
                    t = np.vstack((t, tmp_t))

            except:
                pass

            print t.shape
    
    # calc mean speed


    third = int(0.33*t.shape[1])
    sixth = int(0.16*t.shape[1])
    print sixth, third
    mean_speed = np.mean(speed[:,sixth:third], axis=1)
    print 'shape: ', t.shape, speed.shape, mean_speed.shape
    print
    indices_to_sort = np.argsort(mean_speed)

    response = (odor-control)[indices_to_sort, :]

    baseline_response = np.mean(response[:,sixth:third], axis=1)
    baseline_response = np.reshape(baseline_response, [len(baseline_response), 1])
    baseline_response_repeated = np.repeat(baseline_response, response.shape[1], axis=1)

    #return t, speed, response

    response = response - baseline_response

    odor_on = int(t.shape[1]/3.)
    odor_middle = int(t.shape[1]/3.*1.5)
    odor_off = int(t.shape[1]/3.*2)

    sp = np.mean(speed[indices_to_sort][:,sixth:third], axis=1)

    if use_first_half_of_odor_presentation:
        re = np.mean(response[:,odor_on:odor_middle], axis=1)
    else:
        re = np.mean(response[:,odor_on:odor_off], axis=1)
    

    return response, speed, indices_to_sort, sp, re, experiment_type, t

def plot_odor_response_speed_and_time_matrix(directories, flowrates=1, localtimerange=[14,36], average_within_paths=False, layout=None, 
    flowrate_label='flowrate_a', use_speed_intercept=None, show_regression_results=False,
    minimum_speed=0, # 2017 08 23, tried min_speed of 0.25, makes little difference in CO2 low speed results
    odor='co2',
    speed_layout=None):

    print '********'
    print directories
    print flowrates
    print localtimerange
    print '********'
    response, speed, indices_to_sort, sp, re, experiment_type, t = get_new_attraction_index(directories, flowrates, localtimerange, average_within_paths)

    slope, y_intercept, rval, pval, stderr = scipy.stats.linregress(sp, re)
    x_intercept = -1*y_intercept / slope

    if show_regression_results:
        svgitem = 'regression_' + flowrate_label
        layout.svgitems[svgitem].style['font-size'] = 8
        s = "rsq=%(rsq)0.2f \n pval=%(pval)0.3f intcpt=%(intercept)0.3f" % {'rsq': rval**2, 'pval':pval, 'intercept':x_intercept}
        layout.svgitems[svgitem].text = s

    experiment_type = np.array(experiment_type)[indices_to_sort]
    experiment_type = np.reshape(experiment_type, [len(experiment_type),1])


    if flowrates==0:
        x_intercept_index = int(len(sp)/2.)
    else:
        x_intercept_index = np.argmin(np.abs(sp-x_intercept))
        x_intercept_index_minspeed = np.argmin(np.abs(sp-minimum_speed))

    if use_speed_intercept is not None:
        x_intercept = use_speed_intercept
        x_intercept_index = np.argmin(np.abs(sp-x_intercept))
        x_intercept_index_minspeed = np.argmin(np.abs(sp-minimum_speed))

    #return t, speed, response
    if layout is not None:
        if (flowrate_label, 'correlation') in layout.axes:
            ax = layout.axes[(flowrate_label, 'correlation')]
            sp_indices = np.array([i for i in range(len(sp))])
            if 'single' not in ''.join(directories):
                ax.scatter(re, sp_indices, marker='.', c=re, cmap='seismic', vmax=7, vmin=-7, linewidth=0.2)
            else:
                ax.scatter(re, sp_indices, marker='.', c=re, cmap='seismic', vmax=0.7, vmin=-0.7, linewidth=0.2)

            
            x = np.linspace(np.min(sp), np.max(sp), len(sp_indices))
            y = slope*x + y_intercept
            x_indices = []
            for _x_ in x:
                ix = np.argmin(np.abs(sp-_x_))
                x_indices.append(ix)
            ax.plot(y,x_indices,color='black')
            #ax.vlines(0, 0, len(sp_indices), linewidth=0.5)

            if np.min(re)<-7:
                lo = -10
            else:
                lo = -7

            ax.set_ylim(0,len(sp_indices))
            ax.hlines(x_intercept_index,lo,7,color='magenta')
            if 'single' not in ''.join(directories):
                ax.set_xlim(lo,7)
                figurefirst.mpl_functions.adjust_spines(ax, ['top'], xticks=[lo,0,7], spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5)
                ax.tick_params(length=2.5)
            else:
                ax.set_xlim(lo/10.,0.7)
                figurefirst.mpl_functions.adjust_spines(ax, ['top'], xticks=[lo/10.,0,0.7], spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5)
                ax.tick_params(length=2.5)
            ax.set_xticklabels([lo/10., 0, 0.7])
            flytext.set_fontsize(ax.figure, 6)

        try:
            ax = layout.axes[(flowrate_label, 'nflies')]
        except:
            ax = None   
        if ax is not None:
            if 'single' not in ''.join(directories):
                ax.imshow(response, extent=[-10,20,0,response.shape[0]], aspect='auto', origin='lower', interpolation='nearest', cmap='bwr', vmax=7, vmin=-7)
            else:
                ax.imshow(response, extent=[-10,20,0,response.shape[0]], aspect='auto', origin='lower', interpolation='nearest', cmap='bwr', vmax=0.7, vmin=-0.7)
            
            ax.vlines(0, 0, response.shape[0], color='lime')
            ax.vlines(10, 0, response.shape[0], color='lime')
            ax.hlines(x_intercept_index,-10,20,color='magenta') 
            if minimum_speed > 0:
                ax.hlines(x_intercept_index_minspeed,-10,20,color='cyan')
            ax.set_xlim(-10,20)
            #ax.set_xlabel('time, min')
            #ax.set_ylabel('speed, ordered (top = high)')
            #return t, speed, response
            figurefirst.mpl_functions.adjust_spines(ax, [])

        if (flowrate_label, 'experiment_type') in layout.axes:
            ax = layout.axes[(flowrate_label, 'experiment_type')]
            ax.imshow(experiment_type, extent=[0,1,0,response.shape[0]], aspect='auto', origin='lower', interpolation='nearest', cmap='jet', vmax=2, vmin=0)
            figurefirst.mpl_functions.adjust_spines(ax, [])

        try:
            ax = layout.axes[(flowrate_label, 'speed')]
        except:
            ax = None   
        if ax is not None:
            ax.imshow(speed[indices_to_sort], extent=[-10,20,0,response.shape[0]], aspect='auto', origin='lower', interpolation='nearest', cmap='hot', vmax=8, vmin=0)
            wn = 0.01
            b, a = scipy.signal.butter(3, wn)
            sp_array = np.repeat(sp.reshape(len(sp),1), speed.shape[1], axis=1)
            print '************************'
            print speed[indices_to_sort].shape, sp.shape
            speed_sorted = speed[indices_to_sort]
            indices_where_speed_is_not_too_high_or_low = np.where( (sp<12)*(sp>x_intercept) ) # only use speed in the ok range
            _speed = np.mean(speed_sorted[indices_where_speed_is_not_too_high_or_low]-sp_array[indices_where_speed_is_not_too_high_or_low], axis=0)
            _speed = scipy.signal.filtfilt(b,a,_speed)
            _speed /= np.max(_speed)
            _speed *= 0.3*response.shape[0]
            _speed += 0.15*response.shape[0]
            #ax.plot(t[0,:]/60., _speed, color='white', linewidth=0.5)
            ax.set_ylim(0,response.shape[0])
            ax.vlines(0, 0, response.shape[0], color='lime')
            ax.vlines(10, 0, response.shape[0], color='lime')
            ax.hlines(x_intercept_index,-10,20,color='magenta')
            if minimum_speed > 0:
                ax.hlines(x_intercept_index_minspeed,-10,20,color='cyan')
            #ax.set_ylabel('speed, ordered\n(top = high)')
            #ax.set_xlabel('time, min')
            ax.set_xlim(-10,20)
            figurefirst.mpl_functions.adjust_spines(ax, [])

        if speed_layout is not None:
            ax = speed_layout.axes[(flowrate_label, 'high_speed')]
            wn = 0.01
            b, a = scipy.signal.butter(3, wn)
            sp_array = np.repeat(sp.reshape(len(sp),1), speed.shape[1], axis=1)
            speed_sorted = speed[indices_to_sort]
            indices_where_speed_is_not_too_high_or_low = np.where( (sp<12)*(sp>x_intercept) ) # only use speed in the ok range

            _speed = speed_sorted[indices_where_speed_is_not_too_high_or_low]-sp_array[indices_where_speed_is_not_too_high_or_low]
            _speed_lo, _speed_hi = plot_n_flies.get_95_confidence_intervals(_speed, iterations=100)
            _speed_mean = np.mean(_speed, axis=0)

            _speed_mean = scipy.signal.filtfilt(b,a,_speed_mean)
            _speed_lo = scipy.signal.filtfilt(b,a,_speed_lo)
            _speed_hi = scipy.signal.filtfilt(b,a,_speed_hi)

            ax.plot(t[0,:], _speed_mean, color='black', linewidth=0.5, zorder=1)
            ax.fill_between(t[0,:], _speed_lo, _speed_hi, facecolor='gray', alpha=0.2, edgecolor='none', zorder=0)
            ax.fill_between([0, 600], -1, 2, facecolor='green', edgecolor='none', alpha=0.3, zorder=-1)

            ax.set_xlim(-10*60,20*60)
            ax.set_ylim(-0.5,1.3)
            #figurefirst.mpl_functions.adjust_spines(ax, ['left'], yticks=[-0.5, 0, 1])
            figurefirst.mpl_functions.adjust_spines(ax, [])
            ax = speed_layout.axes[(flowrate_label, 'low_speed')]
            figurefirst.mpl_functions.adjust_spines(ax, [])
            layout.append_figure_to_layer(speed_layout.figures[flowrate_label], flowrate_label, cleartarget=True)

        # high speed
        high_speed_response = np.mean(response[x_intercept_index:,:], axis=0)
        ax = layout.axes[(flowrate_label, 'high_speed')]
        #ax.plot(t[0,:], high_speed_response)
        if 'single' not in ''.join(directories):
            color_mean = {'color': high_speed_response,
                      'cmap': 'seismic',
                      'norm': [-7,7]}
        else:
            color_mean = {'color': high_speed_response,
                      'cmap': 'seismic',
                      'norm': [-0.7,0.7]}
        
        plot_n_flies.plot_mean_and_95_confidence(ax, t[x_intercept_index:,:], response[x_intercept_index:,:], 
            (0.001, 0.001, 0.001), color_mean=color_mean,
            zorder=10)
        ax.set_rasterization_zorder(12)
        ax.fill_between([0, 600], -7, 7, facecolor='green', edgecolor='none', alpha=0.3)
        if 'single' not in ''.join(directories):
            ax.set_ylim(-7,7)
        else:
            ax.set_ylim(-0.7,0.7)
        ax.set_xlim(-10*60,20*60)
        figurefirst.mpl_functions.adjust_spines(ax, [])

        # low speed
        low_speed_response = np.mean(response[x_intercept_index_minspeed:x_intercept_index,:], axis=0)
        ax = layout.axes[(flowrate_label, 'low_speed')]
        if odor == 'co2':
            #ax.plot(t[0,:], low_speed_response)
            if 'single' not in ''.join(directories):
                color_mean = {'color': low_speed_response,
                          'cmap': 'seismic',
                          'norm': [-7,7]}
            else:
                color_mean = {'color': low_speed_response,
                          'cmap': 'seismic',
                          'norm': [-0.7,0.7]}
            plot_n_flies.plot_mean_and_95_confidence(ax, t[x_intercept_index_minspeed:x_intercept_index,:], 
                response[x_intercept_index_minspeed:x_intercept_index,:], 
                (0.001, 0.001, 0.001), color_mean=color_mean,
                zorder=10)
            ax.set_rasterization_zorder(12)
            ax.fill_between([0, 600], -7, 7, facecolor='green', edgecolor='none', alpha=0.3)
            if 'single' not in ''.join(directories):
                ax.set_ylim(-7,7)
            else:
                ax.set_ylim(-0.7,0.7)
            ax.set_xlim(-10*60,20*60)
        figurefirst.mpl_functions.adjust_spines(ax, [])


        layout.append_figure_to_layer(layout.figures[flowrate_label], flowrate_label, cleartarget=True)
        #figure_output = os.path.join(directories[0], 'speed_time_nflies_automatic_figure_output.svg')
        #layout.write_svg(figure_output)

    print 'X INTERCEPT: ', x_intercept
    s = "rsq=%(rsq)0.2f \n pval=%(pval)0.3f intcpt=%(intercept)0.3f" % {'rsq': rval**2, 'pval':pval, 'intercept':x_intercept}
    print s

    #return speed[indices_to_sort], response

def get_paper_layout():
    svg = 'speed_time_nflies_scatter_automatic_summary.svg'
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout

def plot_odor_response_speed_and_time_matrix_for_all_flowrates(directories, localtimerange=[14,36], average_within_paths=False, use_speed_intercept=None):
    if type(directories) is not list:
        directories = [directories]

    layout = get_paper_layout()
    
    pd_data_sets = []

    for directory in directories:
        paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='day')
        for path in paths:
            pd_filename = mta.read_hdf5_file_to_pandas.get_filename(path, 'pd_data.pickle')
            if pd_filename is None:
                s = 'Cannot find pd_data.pickle file. Run analyze on path:\n' + path
                raise ValueError(s)
            pd_data = pandas.read_pickle(pd_filename)
            pd_data_sets.append(pd_data)
    pd_data = pandas.concat(pd_data_sets)
    
    flowrates = pd_data.flowrate.unique()
    flowrates.sort()
    flowrate_axes = ['flowrate_a', 'flowrate_b', 'flowrate_c', 'flowrate_d']

    if len(flowrates) > 4:
        test_flowrates = [0, 0.4, 1, 3.5]
        actual_flowrates = []
        for f in test_flowrates:
            if f in flowrates:
                actual_flowrates.append(f)
        flowrates = actual_flowrates

    for i, flowrate in enumerate(flowrates):
        plot_odor_response_speed_and_time_matrix(directories, flowrates=flowrate, localtimerange=localtimerange, average_within_paths=False, layout=layout, flowrate_label=flowrate_axes[i], use_speed_intercept=use_speed_intercept, show_regression_results=True)
        svgitem = 'text_' + flowrate_axes[i]
        layout.svgitems[svgitem].style['font-size'] = 8
        layout.svgitems[svgitem].text = "%0.1f" % flowrate
        layout.append_figure_to_layer(layout.figures[flowrate_axes[i]], flowrate_axes[i], cleartarget=True)
    

    figure_output = os.path.join(directories[0], 'speed_time_nflies_automatic_figure_output.svg')
    layout.apply_svg_attrs()
    layout.write_svg(figure_output)


def plot_cross_directory_figure(directory, localtimerange=[14,36], average_within_paths=False,use_speed_intercept=None):

    paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='', does_not_contain=['~', '.pyc', 'day','.svg','tmp'])
    if len(paths) > 0:
        plot_odor_response_speed_and_time_matrix_for_all_flowrates(paths, localtimerange, average_within_paths,use_speed_intercept)
    else:
        plot_odor_response_speed_and_time_matrix_for_all_flowrates(directory, localtimerange, average_within_paths,use_speed_intercept)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--path", type="str", dest="path", default='',
                        help="path to data")
    parser.add_option("--intercept", type="float", dest="intercept", default=None,
                        help="speed intercept to use")
    (options, args) = parser.parse_args()  
    
    # note: to combine multiple directories you must call function directly (e.g. from ipython)
    plot_cross_directory_figure(options.path, use_speed_intercept=options.intercept)