from optparse import OptionParser
import sys, os
import imp

import pickle

import time
import numpy as np
import scipy.stats

from multi_tracker_analysis import read_hdf5_file_to_pandas
from multi_tracker_analysis import data_slicing

import matplotlib.pyplot as plt

import multi_tracker_analysis as mta
import copy

import fly_plot_lib.plot as fpl
import flystat
from fly_plot_lib.colormaps import viridis as viridis

import pandas

import data_fit

import orchard.annotated_time_on_pad_analysis as ann
import orchard

import figurefirst

import multicat_analysis as mcat

import co2_paper_locations

AVERAGE_WITHIN_PATHS = False

def get_paper_layout():
    svg = co2_paper_locations.figure_template_locations.figure4_walking_arena_activity
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout

def update_all_and_save():
    layout = get_paper_layout()
    try:
        layout.clear_fflayer('pad')
    except:
        pass
        

    plot_lengths('co2_60sccm', layout.axes_groups['wind_co2']['scatter'])
    plot_lengths('eth_60sccm', layout.axes_groups['wind_eth']['scatter'])

    # co2
    co2_directory = co2_paper_locations.data_locations.walking_arena_activity['low_flow_co2_dusk'] 
    portion_of_day = 'dusk'
    mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('wind_walking_co2', 'nflies')], co2_directory, 'dusk', 1, average_within_paths=AVERAGE_WITHIN_PATHS, normalize_to_one=True, traces_to_plot=['odor'], show_confidence=False)
    plot_bootstrapped_nflies_windtunnel('co2_60sccm', layout.axes_groups['wind_walking_co2']['nflies'])
    figurefirst.mpl_functions.adjust_spines(layout.axes_groups['wind_walking_co2']['scatter'], [])
        
    # ethanol
    eth_directory = co2_paper_locations.data_locations.walking_arena_activity['ethanol_dusk']
    portion_of_day = 'dusk'
    mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('wind_walking_eth', 'nflies')], eth_directory, 'dusk', 1, average_within_paths=AVERAGE_WITHIN_PATHS, normalize_to_one=True, traces_to_plot=['odor'], show_confidence=False)
    plot_bootstrapped_nflies_windtunnel('eth_60sccm', layout.axes_groups['wind_walking_eth']['nflies'])
    figurefirst.mpl_functions.adjust_spines(layout.axes_groups['wind_walking_eth']['scatter'], [])

    # walking comparison
    mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('walking_comparison', 'nflies')], co2_directory, 'dusk', 1, average_within_paths=AVERAGE_WITHIN_PATHS, normalize_to_one=True, traces_to_plot=['odor'], show_confidence=False, colors={'odor': 'magenta', 'control':'blue'})
    mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('walking_comparison', 'nflies')], eth_directory, 'dusk', 1, average_within_paths=AVERAGE_WITHIN_PATHS, normalize_to_one=True, traces_to_plot=['odor'], show_confidence=False, colors={'odor': 'black', 'control':'blue'}, fill_odor=False)
    figurefirst.mpl_functions.adjust_spines(layout.axes_groups['walking_comparison']['scatter'], [])

    # wind comparison
    mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('wind_eth', 'nflies')], co2_directory, 'dusk', 1, average_within_paths=AVERAGE_WITHIN_PATHS, normalize_to_one=True, traces_to_plot=['odor'], show_confidence=False, colors={'odor': 'none', 'control':'blue'})
    mcat.plot_n_flies.plot_n_flies_in_control_and_odor_regions(layout.axes[('wind_co2', 'nflies')], co2_directory, 'dusk', 1, average_within_paths=AVERAGE_WITHIN_PATHS, normalize_to_one=True, traces_to_plot=['odor'], show_confidence=False, colors={'odor': 'none', 'control':'blue'})
    plot_bootstrapped_nflies_windtunnel('eth_60sccm', layout.axes_groups['wind_eth']['nflies'])
    plot_bootstrapped_nflies_windtunnel('co2_60sccm', layout.axes_groups['wind_co2']['nflies'])
    #figurefirst.mpl_functions.adjust_spines(layout.axes_groups['wind_comparison']['scatter'], [])

    layout.append_figure_to_layer(layout.figures['wind_walking_eth'], 'wind_walking_eth', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['wind_walking_co2'], 'wind_walking_co2', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['walking_comparison'], 'walking_comparison', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['wind_co2'], 'wind_co2', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['wind_eth'], 'wind_eth', cleartarget=True)
    layout.write_svg(co2_paper_locations.figure_template_locations.figure4_walking_arena_activity)
    
def plot_lengths(label, ax):
    data_dict_lengths = ann.load_lengths()
    color_dict = {'orco_co260': 'lime',
                  'co2_5sccm': 'magenta',
                'co2_15sccm': 'magenta',
                'co2_60sccm': 'magenta',
                'co2_hungry': 'magenta',
             'eth_60sccm': (0.001, 0.001, 0.001),
             'co2_200sccm': 'magenta',
             'eth_200sccm': (0.001, 0.001, 0.001),
             'eth_15sccm': (0.001, 0.001, 0.001),
             'eth_60_co2_15': (0.001,0.001,0.001),
             'vinegar_15sccm': 'brown',
             'vinegar_60sccm': 'brown',
             'vinegar_200sccm': 'brown',
             'vinegar_60sccm_15co2': (0.001,0.001,0.001),
             'h2o_15sccm': 'blue',
             'h2o_60sccm': 'blue',
             'h2o_200sccm': 'blue',
             'ethacet_15_50': 'green',
             'ethyl_acetate_15sccm': 'green',
             'ethyl_acetate_60sccm': 'green'}
     
    print 'Plotting: ', label
    length = np.array(data_dict_lengths[label]) / 30. / 60. + 75/60. # start at peak time
    length = length.tolist()
    fpl.scatter_box(ax, 1, np.array(length), xwidth=1, color=color_dict[label], shading='95conf', markersize=0.3, edgecolor=color_dict[label], flipxy=True)

    ax.set_ylim(0.5,1.5)
    ax.set_xlim(-10,20)
    figurefirst.mpl_functions.adjust_spines(ax, [])

def plot_maxn(directory, experiment_routines, group, layout, name, portion_of_day, average_across_trials, fill_ranges, fill_colors, ylim):
    orchard.plot_max_nflies.plot_max_in_regions_from_directory_experiment_routine_analysis(
                layout.axes_groups[group][name]['axis'], 
                directory, 
                portion_of_day=portion_of_day, 
                experiment_routines=experiment_routines, 
                include='day', 
                analysis='mean', 
                average_across_trials=average_across_trials, 
                fill_ranges=fill_ranges, 
                fill_colors=fill_colors,
                xlim=[-10,20],
                ylim=ylim,
                show_confidence_intervals=True,
                odor_color='red',
                odor_only=True,
                normalize_to_one=True)
    figurefirst.mpl_functions.adjust_spines(layout.axes_groups[group][name]['axis'], [], xticks=[-10,0,10,20], yticks=[0,1])

def plot_bootstrapped_nflies_windtunnel(label, ax, peak_time=75):
    data_dict_lengths = ann.load_lengths()
    color_dict = {'orco_co260': (0.001, 0.001, 0.001),
                  'co2_5sccm': 'magenta',
                'co2_15sccm': 'magenta',
                'co2_60sccm': 'magenta',
                'co2_hungry': 'magenta',
             'eth_60sccm': (0.001, 0.001, 0.001),
             'co2_200sccm': 'magenta',
             'eth_200sccm': (0.001, 0.001, 0.001),
             'eth_15sccm': (0.001, 0.001, 0.001),
             'eth_60_co2_15': (0.001,0.001,0.001),
             'vinegar_15sccm': 'brown',
             'vinegar_60sccm': 'brown',
             'vinegar_200sccm': 'brown',
             'vinegar_60sccm_15co2': (0.001,0.001,0.001),
             'h2o_15sccm': 'blue',
             'h2o_60sccm': 'blue',
             'h2o_200sccm': 'blue',
             'ethacet_15_50': 'green',
             'ethyl_acetate_15sccm': 'green',
             'ethyl_acetate_60sccm': 'green'}
     
    print 'Plotting: ', label

    # note: co2 tag is a general placeholder, actual data determined by label
    def get_bootstrapped_data(label, data_dict):
        co2_data = np.array(data_dict[label])
        co2_data.sort()
        # bootstrap confidence intervals
        co2_n_bootstrapped = []
        co2_t_bootstrapped = []
        for i in range(1000):
            indices = np.random.randint(0,len(co2_data), len(co2_data))
            indices.sort()
            co2_n = [len(indices)]
            co2_t = [0]
            for trajec in co2_data[indices]:
                co2_n.append(co2_n[-1]-1)
                co2_t.append(trajec/30.)
            co2_n = np.array(co2_n)
            co2_n = co2_n / float(np.max(co2_n))
            co2_time_peak = peak_time # comes from peak time at 75 seconds
            co2_t = np.array(co2_t) + co2_time_peak
            co2_n_bootstrapped.append(co2_n)
            co2_t_bootstrapped.append(co2_t)
        co2_n_bootstrapped_interpolated = []
        co2_t_bootstrapped_interpolated = np.arange(0, int(co2_data[-1]/30.)+1,1)
        for i, c in enumerate(co2_n_bootstrapped):
            c_interpolated = np.interp(co2_t_bootstrapped_interpolated, co2_t_bootstrapped[i], co2_n_bootstrapped[i])
            co2_n_bootstrapped_interpolated.append(c_interpolated)
        co2_n_bootstrapped_interpolated = np.array(co2_n_bootstrapped_interpolated)
        co2_n_bootstrapped_interpolated = np.sort(co2_n_bootstrapped_interpolated, axis=0)
        print int(co2_n_bootstrapped_interpolated.shape[0]*0.025), int(co2_n_bootstrapped_interpolated.shape[0]*0.975)
        co2_conf_low = co2_n_bootstrapped_interpolated[ int(co2_n_bootstrapped_interpolated.shape[0]*0.025) ]
        co2_conf_high = co2_n_bootstrapped_interpolated[ int(co2_n_bootstrapped_interpolated.shape[0]*0.975) ]
        return co2_t_bootstrapped_interpolated, co2_n_bootstrapped_interpolated, co2_conf_low, co2_conf_high
        
    co2_t_bootstrapped_interpolated, co2_n_bootstrapped_interpolated, co2_conf_low, co2_conf_high = get_bootstrapped_data(label, data_dict_lengths)

    ax.plot(co2_t_bootstrapped_interpolated, np.mean(co2_n_bootstrapped_interpolated, axis=0), color=color_dict[label], linewidth=1)
    ax.fill_between(co2_t_bootstrapped_interpolated, co2_conf_low, co2_conf_high, edgecolor='none', facecolor=color_dict[label], alpha=0.3, zorder=-100)





