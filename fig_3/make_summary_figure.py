import figurefirst
import multi_tracker_analysis as mta
import matplotlib.pyplot as plt
import matplotlib
import fly_plot_lib.plot as fpl
import flystat
import fly_plot_lib.text as flytext
import os
import pickle
import pandas
import scipy.stats
import flystat.resampling
import numpy as np
from multicat_analysis import plot_nflies_vs_speed_scatter

import multicat_analysis as mcat

import co2_paper_locations
from co2_paper_locations import data_locations
ffname_to_directory = data_locations.walking_arena_activity

ffname_to_flowrate =      {'high_flow_co2_dusk': 5,
                           'high_flow_warm_dusk': 5,
                           'high_flow_arista_dusk': 5,
                           'high_flow_ethanol_dusk': 5,
                           'balanced_flow_0': 0,
                           'balanced_flow_1': 1,
                           }

ffname_to_time_of_day = {'example': 'dusk',
                         'balanced_flow_0': 'dusk',
                         'balanced_flow_1': 'dusk',
                         'singlepulse_2hrs_nflies': 'firstpulse',
                         'singlepulse_2min_nflies': [14,36]}

ffname_to_scatter_time = { 'high_flow_co2_dusk': None,
                           'low_flow_co2_starved_afternoon': None,
                           'low_flow_co2_afternoon': None,
                           'low_flow_co2_night': None, 
                           'high_flow_co2_dusk': None,
                           'singlepulse_2hrs_nflies': None,
                           'singlepulse_2min_nflies': None,
                            }

ffname_to_speed_scatter_time = { 'placeholder': None,
                                 #'singlepulse_2hrs_nflies': None,
                                 #'singlepulse_2min_nflies': None,
                            }

ffname_to_color =     {'high_flow_co2_dusk': 0,
                       'high_flow_ethanol_dusk': 1,
                       'low_flow_ethanol_dusk': 1,
                       'low_flow_co2_dusk': 0,
                       'high_flow_arista_dusk': 2,
                       'low_flow_arista_dusk': 2,
                       'high_flow_warm_dusk': 3,
                       'low_flow_warm_dusk': 3,
                       'afternoon': 0,
                       'dusk': 0,
                       'night': 0,
                       'morning': 0,
                       'low_flow_ethanol': 0,
                       '24hr_starved_afternoon': 0,
                       'ethanol_afternoon': 0,
                       'ethanol_dusk': 0,
                       }
                       
                       
AVERAGE_WITHIN_PATHS = False

def plot_cmaps():
  layout = get_paper_layout()

  ax_cool = layout.axes[('cmaps', 'cool_cmap')]
  fpl.colorbar(ax_cool, ticks=None, ticklabels=None, colormap='cool', aspect='auto', orientation='horizontal', show_spine=False)

  ax_summer = layout.axes[('cmaps', 'summer_cmap')]
  fpl.colorbar(ax_summer, ticks=None, ticklabels=None, colormap='summer', aspect='auto', orientation='horizontal', show_spine=False)

  layout.append_figure_to_layer(layout.figures['cmaps'], 'cmaps', cleartarget=True)
  layout.write_svg(co2_paper_locations.figure_template_locations.figure4_walking_arena_activity)


def check_paths():
  for ffname, directory in ffname_to_directory.items():
    try:
      paths = mta.read_hdf5_file_to_pandas.get_filenames(ffname_to_directory[ffname], contains='day') 
      print directory, ': ', len(paths)
    except:
      print 'FAILED: ', directory

def get_paper_layout():
    svg = co2_paper_locations.figure_template_locations.figure4_walking_arena_activity
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    fifidatafile  = os.path.join(os.path.dirname(svg), 'figure_data.pickle')
    layout.fifidatafile = fifidatafile
    return layout
    
def print_mean_speed_for_low_flow_dusk():
    scatter_data = mcat.plot_n_flies.get_nflies_and_speed_at_time_for_directory(ffname_to_directory['low_flow_co2'], 'dusk', 1, -5*60)
    print 'Speed: ', np.mean(np.hstack(scatter_data['speed']))
    
def plot_circadian():
    layout = get_paper_layout()
    fifidatafile = layout.fifidatafile
    
    mcat.plot_circadian.plot_circadian([layout, 'control_circadian', u'control_circadian'], ffname_to_directory['control_dusk'])
    mcat.plot_circadian.plot_circadian([layout, 'hcs_circadian', u'hcs_circadian'], ffname_to_directory['low_flow_co2_dusk'])
    mcat.plot_circadian.plot_circadian([layout, 'starved_circadian', u'starved_circadian'], ffname_to_directory['low_flow_co2_starved_dusk'])
    mcat.plot_circadian.plot_circadian([layout, 'ethanol_circadian', u'ethanol_circadian'], ffname_to_directory['ethanol_starved_dusk'])
    
    # color time ranges            
    paths = mta.read_hdf5_file_to_pandas.get_filenames(ffname_to_directory['low_flow_co2_dusk'], contains='day') 
    config = mta.read_hdf5_file_to_pandas.load_config_from_path(paths[0])
    axes = [[layout, 'control_circadian_colors', 'control_circadian_colors'],
            [layout, 'hcs_circadian_colors', 'hcs_circadian_colors'],
            [layout, 'starved_circadian_colors', 'starved_circadian_colors'],
            [layout, 'ethanol_circadian_colors', 'ethanol_circadian_colors'],
            ]
    portion_of_day_to_color = {'afternoon': 'white',
                               'dusk': 'black',
                               'night': 'white',
                               'morning': 'black',
                               }
    for ax in axes:
        for portion_of_day, time_range in config.portion_of_day_to_local_time_range.items():
            #ax.fill_between(time_range, 0, 1, facecolor=portion_of_day_to_color[portion_of_day], edgecolor='black')\
            layout, figure, axis = ax
            figurefirst.deprecated_regenerate.mpl('fill_between', layout, figure, axis, fifidatafile, 'color portion of day '+portion_of_day, 
                                      [], 
                                      time_range, 0, 1, facecolor=portion_of_day_to_color[portion_of_day], edgecolor='black')

        #ax.set_xlim(15,35)
        figurefirst.deprecated_regenerate.mpl('set_xlim', layout, figure, axis, fifidatafile, 'xlim', [], 
                              15,35)
        #ax.set_ylim(0,1)
        figurefirst.deprecated_regenerate.mpl('set_ylim', layout, figure, axis, fifidatafile, 'ylim', [], 
                              0,1)
        figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust ABV spines', [], 
                              [])

    
    layout.append_figure_to_layer(layout.figures['control_circadian'], 'control_circadian', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['control_circadian_colors'], 'control_circadian_colors', cleartarget=True)
    
    layout.append_figure_to_layer(layout.figures['hcs_circadian'], 'hcs_circadian', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['hcs_circadian_colors'], 'hcs_circadian_colors', cleartarget=True)
    
    layout.append_figure_to_layer(layout.figures['starved_circadian'], 'starved_circadian', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['starved_circadian_colors'], 'starved_circadian_colors', cleartarget=True)
    
    layout.append_figure_to_layer(layout.figures['ethanol_circadian'], 'ethanol_circadian', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['ethanol_circadian_colors'], 'ethanol_circadian_colors', cleartarget=True)
    layout.write_svg(co2_paper_locations.figure_template_locations.figure4_walking_arena_activity)
    
    
def plot_time_series_and_scatters(figurename, directory, flowrate, time_of_day, scatter_time_secs=75, speed_scatter_time=-300): # scatter_time_secs=75, speed_scatter_time=-300
    layout = get_paper_layout()
    fifidatafile = layout.fifidatafile
    print directory
    mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions([layout, figurename, 'nflies'], directory, time_of_day, flowrate, average_within_paths=AVERAGE_WITHIN_PATHS)
    
    print directory
    print flowrate, time_of_day
    if scatter_time_secs is not None:
        scatter_data = mcat.plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, time_of_day, flowrate, scatter_time_secs)
        mcat.plot_n_flies.plot_raw_scatter(layout.axes[(figurename, 'nflies_scatter_odor')],       scatter_data, 'n_flies_odor',       average_within_paths=False)
        mcat.plot_n_flies.plot_raw_scatter(layout.axes[(figurename, 'nflies_scatter_control')],    scatter_data, 'n_flies_control',    average_within_paths=False)
    
    if speed_scatter_time is not None:
        scatter_data = mcat.plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, time_of_day, flowrate, speed_scatter_time)
        mcat.plot_n_flies.plot_raw_scatter([layout, figurename, 'speed_scatter'],             scatter_data, 'speed',              average_within_paths=False, hide_markers=True)
        ax = layout.axes[(figurename, 'nflies')]
        #draw_triangle(ax, speed_scatter_time, 'black') # speed
    else:
        ax = layout.axes[(figurename, 'speed_scatter')]
        figurefirst.mpl_functions.adjust_spines(ax, [])
    
    if scatter_time_secs is not None:
        if scatter_time_secs > 150:
            #draw_triangle(ax, scatter_time_secs, 'orange') # speed
            draw_triangle(ax, scatter_time_secs, 'magenta') # speed
        else:
            draw_triangle(ax, scatter_time_secs, 'magenta') # speed

    if 1: # this dashed line
        #ax = layout.axes[(figurename, 'arrow')]
        figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figurename, 'arrow', fifidatafile, 'adjust spines', [], 
                              [])
        if speed_scatter_time is not None:
            #ax = layout.axes[(figurename, 'speed_scatter')]
            #speed_arrow = matplotlib.patches.FancyArrow(1, 3.25, -1, 0, width=1, length_includes_head=True, head_width=0.7, head_length=1, shape='full', overhang=0, head_starts_at_zero=False, facecolor='black', edgecolor='none')
            ax.plot([-0.5, 0.5], [3.3, 3.3], color='black', dashes=[1,1])
            figurefirst.deprecated_regenerate.mpl('plot', layout, figurename, 'speed_scatter', fifidatafile, 'speed dashed line',
                                      [],
                                      [-0.5, 0.5], [3.3, 3.3], color='black', dashes=[1,1])
            #ax.set_xlim(-0.5, 0.5)
            figurefirst.deprecated_regenerate.mpl('set_xlim', layout, figurename, 'speed_scatter', fifidatafile, 'set_xlim', [], 
                              -0.5, 0.5)
            #ax.set_ylim(0,7)
            figurefirst.deprecated_regenerate.mpl('set_ylim', layout, figurename, 'speed_scatter', fifidatafile, 'set_ylim', [], 
                              0,7)
            #figurefirst.mpl_functions.adjust_spines(ax, [])
        else:
            pass
            #ax = layout.axes[(figurename, 'speed_scatter')]
            #figurefirst.mpl_functions.adjust_spines(ax, [])
    
    layout.append_figure_to_layer(layout.figures[figurename], figurename, cleartarget=True)
    layout.write_svg(co2_paper_locations.figure_template_locations.figure4_walking_arena_activity )
    

def update_panels_for_full_day(panel_basename):
    panels = [  panel_basename+'_afternoon', 
                panel_basename+'_dusk', 
                panel_basename+'_night', 
                panel_basename+'_morning', 
                ]
    update_panels(panels)
    
def update_panels(panels):
    
    for panel in panels:
        directory = ffname_to_directory[panel]
        try:
            flowrate = ffname_to_flowrate[panel]
        except:
            flowrate = 1
            
        try:
            time_of_day = ffname_to_time_of_day[panel]
        except:
            time_of_day = panel.split('_')[-1]
        
        try:
            scatter_time = ffname_to_scatter_time[panel]
        except:
            scatter_time = None # 75

        try:
            speed_scatter_time = ffname_to_speed_scatter_time[panel]
        except:
            speed_scatter_time = -300

        print panel
        print directory
        print flowrate
        print time_of_day
        
        plot_time_series_and_scatters(panel, directory, flowrate, time_of_day, scatter_time_secs=scatter_time, speed_scatter_time=speed_scatter_time)
    
def plot_attraction_index_vs_speed():
    layout = get_paper_layout()
    fifidatafile = layout.fifidatafile
    ax = [layout, 'activity_correlation_relative', 'activity_correlation_relative']
    layout, figure, axis = ax
    
    AI = np.array([])
    speed = np.array([])
    
    AI_means = []
    speed_means = []
    
    
    for time_of_day in ['afternoon', 'dusk', 'night', 'morning']:
        AI_tmp, speed_tmp = plot_attraction_index_vs_speed_for_time_of_day(ax, 'high_flow_co2_dusk', 5, time_of_day)
        AI = np.hstack((AI,AI_tmp))
        AI_means.append(np.mean(AI_tmp))
        speed = np.hstack((speed,speed_tmp))
        speed_means.append(np.mean(speed_tmp))
        AI_tmp, speed_tmp = plot_attraction_index_vs_speed_for_time_of_day(ax, 'high_flow_arista_dusk', 5, time_of_day)
        AI = np.hstack((AI,AI_tmp))
        AI_means.append(np.mean(AI_tmp))
        speed = np.hstack((speed,speed_tmp))
        speed_means.append(np.mean(speed_tmp))
        AI_tmp, speed_tmp = plot_attraction_index_vs_speed_for_time_of_day(ax, 'low_flow_co2_dusk', 1, time_of_day)
        AI = np.hstack((AI,AI_tmp))
        AI_means.append(np.mean(AI_tmp))
        speed = np.hstack((speed,speed_tmp))
        speed_means.append(np.mean(speed_tmp))
        AI_tmp, speed_tmp = plot_attraction_index_vs_speed_for_time_of_day(ax, 'low_flow_co2_starved_dusk', 1, time_of_day)
        AI = np.hstack((AI,AI_tmp))
        AI_means.append(np.mean(AI_tmp))
        speed = np.hstack((speed,speed_tmp))
        speed_means.append(np.mean(speed_tmp))
    for time_of_day in ['afternoon', 'dusk']:
        AI_tmp, speed_tmp = plot_attraction_index_vs_speed_for_time_of_day(ax, 'high_flow_warm_dusk', 5, time_of_day)
        AI = np.hstack((AI,AI_tmp))
        AI_means.append(np.mean(AI_tmp))
        speed = np.hstack((speed,speed_tmp))
        speed_means.append(np.mean(speed_tmp))
        AI_tmp, speed_tmp = plot_attraction_index_vs_speed_for_time_of_day(ax, 'low_flow_warm_dusk', 1, time_of_day)
        AI = np.hstack((AI,AI_tmp))
        AI_means.append(np.mean(AI_tmp))
        speed = np.hstack((speed,speed_tmp))
        speed_means.append(np.mean(speed_tmp))
    
    # regression
    if 0:
        slope, intercept, r, p, stderr = scipy.stats.linregress(speed, AI)
        if p < 0.05:
            x = np.linspace( np.min(speed), np.max(speed), 50 )
            y = x*slope + intercept
            ax.plot(x, y, color='red', linewidth=1)
        print 'All data: '
        print slope, intercept, r, p, stderr
    
    # regression 2
    slope, intercept, rval, pval, stderr = scipy.stats.linregress(np.array(speed_means), np.array(AI_means))
    if pval < 0.05:
        x = np.linspace( np.min(speed), np.max(speed), 100 )
        y = x*slope + intercept
        #ax.plot(x, y, color='black', linewidth=1)
        figurefirst.deprecated_regenerate.mpl('plot', layout, figure, axis, fifidatafile, 'plot correlation line', 
                                  ['speed', 'attraction index'], 
                                  x, y, color='black', linewidth=1)
    print 'means: '
    print slope, intercept, rval, pval, stderr

    #ax.set_xlim(0,8)
    figurefirst.deprecated_regenerate.mpl('set_xlim', layout, figure, axis, fifidatafile, 'xlim', [], 
                              0,8)
    #ax.set_ylim(-4,4)
    figurefirst.deprecated_regenerate.mpl('set_ylim', layout, figure, axis, fifidatafile, 'ylim', [], 
                              -4,4)

    #ax.hlines(0,0,8, linestyle='--', color='black', linewidth=0.5)
    figurefirst.deprecated_regenerate.mpl('hlines', layout, figure, axis, fifidatafile, 'hlines', [], 
                              0,0,8, linestyle='--', color='black', linewidth=0.5)
    regression_intersection = np.argmin( np.abs(y-0) )
    speed_intersection = int(x[regression_intersection]*100)/100.
    #ax.vlines(speed_intersection, -4, 4, linestyle='--', color='black', linewidth=0.5)
    figurefirst.deprecated_regenerate.mpl('vlines', layout, figure, axis, fifidatafile, 'vlines', [], 
                              speed_intersection, -4, 4, linestyle='--', color='black', linewidth=0.5)
    #figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], spine_locations={'left': 5, 'bottom': 5}, xticks=[0,speed_intersection, 8], yticks=[-4, 0, 4], linewidth=0.5)
    #ax.tick_params(length=2.5)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust x spines', [], 
                              ['left', 'bottom'], spine_locations={'left': 5, 'bottom': 5}, xticks=[0,speed_intersection, 8], yticks=[-4, 0, 4], linewidth=0.5, tick_length=2.5)

    figurefirst.deprecated_regenerate.mpl('set_yticklabels', layout, figure, axis, fifidatafile, 'y tick labels', [], 
                              [])

    xticklabels = ['0', str(speed_intersection), '8']
    #ax.set_xticklabels(xticklabels)
    figurefirst.deprecated_regenerate.mpl('set_xticklabels', layout, figure, axis, fifidatafile, 'x tick labels', [], 
                              [])
    #flytext.set_fontsize(ax.figure, 6)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, figure, axis, fifidatafile, 'set fontsize', [], 
                                  6)
    
    layout.svgitems['co2_activity_regression'].style['font-size'] = 6
    if pval < 0.001:
        pval = 0.001
    layout.svgitems['co2_activity_regression'].text = 'p<{0:.3f}'.format(pval) + ', r2={0:.2f}'.format(rval**2)
            
    layout.append_figure_to_layer(layout.figures['activity_correlation_relative'], 'activity_correlation_relative', cleartarget=True)
    layout.apply_svg_attrs()
    layout.write_svg(co2_paper_locations.figure_template_locations.figure4_walking_arena_activity )
    

def plot_attraction_index_vs_speed_ethanol():
    layout = get_paper_layout()
    ax = layout.axes['activity_correlation_relative_ethanol', 'activity_correlation_relative_ethanol']
    
    AI = np.array([])
    speed = np.array([])
    
    AI_means = []
    speed_means = []
    
    
    for time_of_day in ['afternoon', 'dusk', 'night', 'morning']:
        AI_tmp, speed_tmp = plot_attraction_index_vs_speed_for_time_of_day(ax, 'low_flow_ethanol', 1, time_of_day)
        AI = np.hstack((AI,AI_tmp))
        AI_means.append(np.mean(AI_tmp))
        speed = np.hstack((speed,speed_tmp))
        speed_means.append(np.mean(speed_tmp))
    
    # regression
    if 1:
        slope, intercept, r, p, stderr = scipy.stats.linregress(speed, AI)
        if p < 0.05:
            x = np.linspace( np.min(speed), np.max(speed), 50 )
            y = x*slope + intercept
            ax.plot(x, y, color='red')
        print 'All data: '
        print slope, intercept, r, p, stderr
    
    # regression 2
    slope, intercept, r, p, stderr = scipy.stats.linregress(np.array(speed_means), np.array(AI_means))
    if p < 0.05:
        x = np.linspace( np.min(speed), np.max(speed), 50 )
        y = x*slope + intercept
        ax.plot(x, y, color='black')
    print 'means: '
    print slope, intercept, r, p, stderr

    ax.set_xlim(0,12)
    ax.set_ylim(-5,10)
    figurefirst.mpl_functions.adjust_spines(ax, [])
    

    layout.append_figure_to_layer(layout.figures['activity_correlation_relative_ethanol'], 'activity_correlation_relative_ethanol', cleartarget=True)
    layout.write_svg(co2_paper_locations.figure_template_locations.figure4_walking_arena_activity )
    
def get_color_from_ffname(ffname):
#    if 'low_flow' in ffname:
#        cmap = matplotlib.cm.get_cmap('winter')
#        color = cmap( ffname_to_color[ffname]/3. )
    if ffname == 'high_flow_co2_dusk':
        color = 'maroon'
    elif ffname == 'high_flow_arista_dusk':
        color = 'purple'
    elif ffname == 'high_flow_warm_dusk':
        color = 'yellow'
    
    if ffname == 'ethanol_starved_dusk':
        color = (0.001, 0.001, 0.001)
        
    if ffname == 'low_flow_co2_starved_dusk':
        color = (0.001, 0.001, 0.001)
    if ffname == 'low_flow_co2_dusk':
        color = 'cyan'
    if ffname == 'low_flow_warm_dusk':
        color = 'teal'
        
    return color
    
def plot_attraction_index_vs_speed_for_time_of_day(ax, ffname, flowrate, time_of_day):
    layout, figure, axis = ax
    fifidatafile = layout.fifidatafile

    directory = ffname_to_directory[ffname]
    AI_and_speed = get_attraction_index_and_speed(directory, flowrate, time_of_day)
    
    time_of_day_to_marker = {'afternoon': "^",
                                'dusk': 's',
                                'night': 'x',
                                'morning': 'o',
                                }
    time_of_day_to_size = {  'afternoon': 30,
                                'dusk': 30,
                                'night': 30,
                                'morning': 40,
                                }
                                
    color = get_color_from_ffname(ffname)
    
    
    if time_of_day is not 'night':
        #ax.scatter(AI_and_speed['speed']['mean'], AI_and_speed['AI']['mean'], facecolor=color, edgecolor='black', linewidth=0.25, marker=time_of_day_to_marker[time_of_day], s=time_of_day_to_size[time_of_day])
        figurefirst.deprecated_regenerate.mpl('scatter', layout, figure, axis, fifidatafile, 'speed vs attraction index '+ ffname.replace('_dusk','') + ' ' + time_of_day,
                                  ['Speed',
                                   'Mean attraction index'],
                                   AI_and_speed['speed']['mean'], AI_and_speed['AI']['mean'], facecolor=color, edgecolor='black', linewidth=0.25, marker=time_of_day_to_marker[time_of_day], s=time_of_day_to_size[time_of_day])
    else:
        #ax.scatter(AI_and_speed['speed']['mean'], AI_and_speed['AI']['mean'], facecolor=color, edgecolor=color, linewidth=1, marker=time_of_day_to_marker[time_of_day], s=time_of_day_to_size[time_of_day])
        figurefirst.deprecated_regenerate.mpl('scatter', layout, figure, axis, fifidatafile, 'speed vs attraction index '+ ffname.replace('_dusk','') + ' ' + time_of_day,
                                  ['Speed',
                                   'Mean attraction index'],
                                   AI_and_speed['speed']['mean'], AI_and_speed['AI']['mean'], facecolor=color, edgecolor=color, linewidth=1, marker=time_of_day_to_marker[time_of_day], s=time_of_day_to_size[time_of_day])
    if 0:    
        el = matplotlib.patches.Ellipse( (AI_and_speed['speed']['mean'], AI_and_speed['AI']['mean']), 
                                        AI_and_speed['speed']['conf_hi']-AI_and_speed['speed']['conf_lo'], 
                                        AI_and_speed['AI']['conf_hi']-AI_and_speed['AI']['conf_lo'],
                                        angle=0.0, facecolor=color, edgecolor='none', alpha=0.3)
        ax.add_artist(el)
    else:
        figurefirst.deprecated_regenerate.mpl_patch('Ellipse', layout, figure, axis, fifidatafile, 'speed vs attraction index 95 perc confidence '+ ffname.replace('_dusk','') + ' ' + time_of_day, 
                                        [], 
                                        (AI_and_speed['speed']['mean'], AI_and_speed['AI']['mean']), 
                                        AI_and_speed['speed']['conf_hi']-AI_and_speed['speed']['conf_lo'], 
                                        AI_and_speed['AI']['conf_hi']-AI_and_speed['AI']['conf_lo'],
                                        angle=0.0, facecolor=color, edgecolor='none', alpha=0.3)

    return AI_and_speed['AI']['raw'], AI_and_speed['speed']['raw']
        
    
def get_attraction_index_and_speed(directory, flowrate, time_of_day):
    
    scatter_data_black = mcat.plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, time_of_day, flowrate, -5*60)
    AI_black = np.array([np.array(s) for s in scatter_data_black['n_flies_odor']])
    speed_black = np.array([np.array(s) for s in scatter_data_black['speed']])
    
    scatter_data_magenta = mcat.plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, time_of_day, flowrate, [0,600])
    AI_magenta = np.array([np.array(s) for s in scatter_data_magenta['n_flies_odor']])
    
    AI = AI_magenta - AI_black
    
    if AVERAGE_WITHIN_PATHS:
        raise NotImplementedError
    
    AI_mean = np.mean(np.hstack(AI))
    AI_conf_lo, AI_conf_hi = flystat.resampling.bootstrap_confidence_intervals_from_data(np.hstack(AI), use='mean')
    
    speed_mean = np.mean(np.hstack(speed_black))
    speed_conf_lo, speed_conf_hi = flystat.resampling.bootstrap_confidence_intervals_from_data(np.hstack(speed_black), use='mean')
    
    return {'AI': {'mean': AI_mean*10, 'conf_lo': AI_conf_lo*10, 'conf_hi': AI_conf_hi*10, 'raw': np.hstack(AI)*10}, # the *10 is because the attraction indices are scaled to 1, but I want to compare to the N=10 flies
            'speed': {'mean': speed_mean, 'conf_lo': speed_conf_lo, 'conf_hi': speed_conf_hi, 'raw': np.hstack(speed_black)},
            }
    
def draw_triangle(ax, location, color):#, nflies_rawindex, index_odor, index_ctrl):
    print 'DRAWING TRIANGLE: ', location, color
    print ax.get_ylim()
    y = ax.get_ylim()[-1]
    nflies_triangle = matplotlib.patches.FancyArrow(location, y, 0, -0.1*y, width=1, length_includes_head=True, head_width=200, head_length=.1*y, shape='full', overhang=0, head_starts_at_zero=False, facecolor=color, edgecolor='none')
    
    #nflies_triangle = matplotlib.patches.RegularPolygon((location,10), 3, radius=5, orientation=np.pi, facecolor=color, edgecolor='none')
    ax.add_artist(nflies_triangle)


def plot_odor_response_speed_and_time_matrix():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    directories = []
    flowrates = []

    ffnames = [ 'high_flow_co2_dusk',
                'high_flow_arista_dusk',
                'low_flow_co2_dusk',
                'low_flow_co2_starved_dusk',
                'high_flow_warm_dusk',
                'low_flow_warm_dusk',
                ]

    for ffname in ffnames:
      directories.append(ffname_to_directory[ffname])
      if ffname in ffname_to_flowrate.keys():
        flowrates.append(ffname_to_flowrate[ffname])
      else:
        flowrates.append(1)

    plot_nflies_vs_speed_scatter.plot_odor_response_speed_and_time_matrix(directories, flowrates=flowrates, localtimerange=[14,36], average_within_paths=False, ax=ax)
    
    
def make_example():
  update_panels(['example'])

  directory = ffname_to_directory['example']
  layout = get_paper_layout()
  mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('example_nflies_odor', 'timetrace')], directory, 'dusk', 1, traces_to_plot=['odor'])
  mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('example_nflies_control', 'timetrace')], directory, 'dusk', 1, traces_to_plot=['control'])
  mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('example_speed', 'timetrace')], directory, 'dusk', 1, traces_to_plot=['speed'])

  speed_scatter_data = mcat.plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, 'dusk', 1, -300)
  mcat.plot_n_flies.plot_raw_scatter(layout.axes[('example_speed', 'scatter')], speed_scatter_data, 'speed',  average_within_paths=False, hide_markers=True)
  draw_triangle(layout.axes[('example_speed', 'timetrace')], -300, 'black')

  odor_scatter_data = mcat.plot_n_flies.get_nflies_and_speed_at_time_for_directory(directory, 'dusk', 1, 75)
  mcat.plot_n_flies.plot_raw_scatter(layout.axes[('example_nflies_odor', 'scatter')], odor_scatter_data, 'n_flies_odor',  average_within_paths=False)
  mcat.plot_n_flies.plot_raw_scatter(layout.axes[('example_nflies_control', 'scatter')], odor_scatter_data, 'n_flies_control',  average_within_paths=False)
  draw_triangle(layout.axes[('example_nflies_odor', 'timetrace')], 75, 'magenta')
  draw_triangle(layout.axes[('example_nflies_control', 'timetrace')], 75, 'magenta')


  for figurename in ['example_nflies_odor', 'example_nflies_control', 'example_speed']:
    layout.append_figure_to_layer(layout.figures[figurename], figurename, cleartarget=True)
  layout.write_svg(co2_paper_locations.figure_template_locations.figure4_walking_arena_activity )

  #plot_cmaps()


def make_single_pulse_plots():

  update_panels(['singlepulse_2hrs_nflies', 'singlepulse_2min_nflies'])


  layout = get_paper_layout()

  mcat.plot_circadian.plot_circadian_for_tmaze_exps([layout, 'singlepulse_2hrs', u'speed'], ffname_to_directory['singlepulse_2hrs'], t_range=[0,2*60+30])
  mcat.plot_circadian.plot_circadian_for_tmaze_exps([layout, 'singlepulse_2min', u'speed'], ffname_to_directory['singlepulse_2min'], t_range=[0,30])
  
  #figurefirst.mpl_functions.adjust_spines(layout.axes[(u'singlepulse_2hrs', u'speed')], [])
  figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, 'singlepulse_2hrs', u'speed', layout.fifidatafile, 'adjust spines', [], 
                              [])
  #figurefirst.mpl_functions.adjust_spines(layout.axes[(u'singlepulse_2min', u'speed')], [])
  figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, 'singlepulse_2min', u'speed', layout.fifidatafile, 'adjust spines', [], 
                              [])

  layout.append_figure_to_layer(layout.figures['singlepulse_2hrs'], 'singlepulse_2hrs', cleartarget=True)
  layout.append_figure_to_layer(layout.figures['singlepulse_2min'], 'singlepulse_2min', cleartarget=True)

  layout.write_svg(co2_paper_locations.figure_template_locations.figure4_walking_arena_activity )

if __name__ == '__main__':

    if 1:
      make_single_pulse_plots()
      make_example()

      plot_circadian()
      #plot_attraction_index_vs_speed()

      panel_basenames = ['control', 'low_flow_co2', 'low_flow_co2_starved', 'ethanol_starved']
      for panel_basename in panel_basenames:
          update_panels_for_full_day(panel_basename)
      
      high_flow_panels = ['high_flow_ethanol_dusk', 'high_flow_warm_dusk', 'high_flow_arista_dusk', 'high_flow_co2_dusk']
      update_panels(high_flow_panels)
    
    #plot_circadian()
    #make_example()