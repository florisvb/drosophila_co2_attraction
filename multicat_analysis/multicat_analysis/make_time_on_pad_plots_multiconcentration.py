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
import fly_plot_lib.text as flytext
import flystat
from fly_plot_lib.colormaps import viridis as viridis

import pandas

import data_fit

import orchard.annotated_time_on_pad_analysis as ann


import figurefirst

import co2_paper_locations
import pairs2groups

def get_paper_layout(version='old'):
    if version == 'old':
        svg = co2_paper_locations.figure_template_locations.figure3_windtunnel_walking
    elif version == 'supplemental':
        'print supplemental figure'
        svg = co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_supp
    else:
        svg = co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_simple
    layout = figurefirst.svg_to_axes.FigureLayout(str(svg))
    layout.make_mplfigures()
    return layout

def plot_all_lengths(version='old'):

    data_dict_lengths = ann.load_lengths()
    data_dict_entries = ann.load_lengths(load_data='odor_entries')
    data_dict_distances = ann.load_lengths(load_data='distances')
    data_dict_speeds = ann.load_lengths(load_data='speeds')
    
    #odors = ['h2o', 'co2', 'eth', 'eth-co2', 'vinegar', 'ethylacetate', 'choice-eth-only', 'choice-eth-and-co2']
    #odors = ['h2o', 'co2', 'eth', 'eth-co2', 'vinegar', 'ethylacetate']
    if version == 'old':
        plot_lengths([data_dict_lengths, data_dict_entries, data_dict_distances, data_dict_speeds])
    else:
        plot_lengths_simple([data_dict_lengths, data_dict_entries, data_dict_distances, data_dict_speeds])

def plot_lengths_simple(data=None):
    fifidatafile = os.path.join(os.path.dirname(co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_simple), 'figure_data.pickle')

    if data is None:
        data_dict_lengths = ann.load_lengths()
        data_dict_entries = ann.load_lengths(load_data='odor_entries')
        data_dict_distances = ann.load_lengths(load_data='distances')
        data_dict_speeds = ann.load_lengths(load_data='speeds')
    else:
        data_dict_lengths, data_dict_entries, data_dict_distances, data_dict_speeds = data
    color_dict = {  'co2_60sccm': 'magenta',
                    'eth_60sccm': (0.001,0.001,0.001),
                    'eth-co2_60_15': 'purple',
                    'vinegar_60sccm': 'brown',
                    'h2o_60sccm': 'blue'}

    labels = ['h2o_60sccm', 'co2_60sccm', 'eth_60sccm', 'eth-co2_60_15', 'vinegar_60sccm']
    #labels = ['h2o_60sccm', 'co2_60sccm', 'eth_60sccm', 'eth-co2_60_15', 'vinegar_60sccm', 'ethylacetate_15sccm']
    labels.reverse()
    
    #'vinegar_60sccm_15co2',  'h2o_60sccm', 'orco_co260', 

    layout = get_paper_layout('new')
    #layout_figure = 'pad_' + odor
    #try:
    #    layout.clear_fflayer(layout_figure)
    #except:
    #    pass

    yticks = []
    
    nlabels = 0

    data_for_stats = { 'time': {},
                       'distance': {},
                       'approaches': {},
                       'speed': {}}
    for x, label in enumerate(labels):
    
        #if label == 'none':
        #    continue
        #if odor + '_' not in label:
        #    continue
        #if label.startswith(odor) is False:
        #    continue
        
        if label[-2:] == '_2':
            ff_label = copy.copy(label)
            label = label[0:-2]
        else:
            ff_label = label
    
        MARKERSIZE = 0.4

        print '**************'
        print label, len(np.array(data_dict_lengths[label]) )
        print '**************'

        def set_lims_spines(layout, figure, axis, label, ymin, ymax, yticks):
            figurefirst.deprecated_regenerate.mpl('set_xlim', layout, figure, axis, fifidatafile, 'xlim', [], 
                                      -0.5, 0.5)
            figurefirst.deprecated_regenerate.mpl('set_ylim', layout, figure, axis, fifidatafile, 'ylim', [], 
                                      ymin, ymax)

            if label == 'h2o_60sccm':
                #figurefirst.mpl_functions.adjust_spines(layout.axes[('time', ff_label)], ['left'], spine_locations={'bottom': 5, 'bottom': 5}, yticks=[0,5,10], linewidth=0.5)
                #layout.axes[('time', ff_label)].tick_params(length=2.5)
                figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust spines', [], 
                              ['left'], yticks=yticks, spine_locations={'left': 5, 'bottom': 5}, linewidth=0.5, tick_length=2.5)
            else:
                #figurefirst.mpl_functions.adjust_spines(layout.axes[('time', ff_label)], [])
                figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust spines', [], 
                              [])
        
        print 'Plotting: ', label, ff_label
        length = np.array(data_dict_lengths[label]) / 30. / 60.
        figure = 'time'
        #fpl.scatter_box(layout.axes[('time', ff_label)], 0, np.array(length), color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, xwidth=0.4)
        figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', layout, figure, ff_label, fifidatafile, 'stats for ' + ff_label + ' ' + figure,
                                      ['index',
                                       'List of ' + figure + ' for ' + ff_label],
                                       0, length, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, xwidth=0.4, use='median')
        set_lims_spines(layout, figure, ff_label, label, 0, 10, [0,5,10])
        data_for_stats['time'][label.split('_')[0]] = length

        distances = np.array(data_dict_distances[label])*30./46
        figure = 'distances'
        #fpl.scatter_box(layout.axes[('distances', ff_label)], 0, distances, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, use='median', xwidth=0.4)
        figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', layout, figure, ff_label, fifidatafile, 'stats for ' + ff_label + ' ' + figure,
                                      ['index',
                                       'List of ' + figure + ' for ' + ff_label],
                                       0, distances, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, xwidth=0.4, use='median')
        set_lims_spines(layout, figure, ff_label, label, 0, 100, [0,50,100])
        data_for_stats['distance'][label.split('_')[0]] = distances

        speeds = np.array(data_dict_speeds[label])*30./46*10 # raw speed is in pixels per frame
        figure = 'speeds'
        #fpl.scatter_box(layout.axes[('speeds', ff_label)], 0, speeds, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, use='median', xwidth=0.4)
        figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', layout, figure, ff_label, fifidatafile, 'stats for ' + ff_label + ' ' + figure,
                                      ['index',
                                       'List of ' + figure + ' for ' + ff_label],
                                       0, speeds, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, xwidth=0.4, use='median')
        set_lims_spines(layout, figure, ff_label, label, 0, 4, [0,2,4])
        data_for_stats['speed'][label.split('_')[0]] = speeds
        
        entries = data_dict_entries[label]
        r = scipy.stats.uniform(-0.25,0.25)
        entries = np.array(entries) + r.rvs(len(entries))
        figure = 'approaches'
        #fpl.scatter_box(layout.axes[('approaches', ff_label)], 0, entries, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, use='mean', xwidth=0.4)
        figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', layout, figure, ff_label, fifidatafile, 'stats for ' + ff_label + ' ' + figure,
                                      ['index',
                                       'List of ' + figure + ' for ' + ff_label],
                                       0, entries, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, xwidth=0.4, use='mean')
        set_lims_spines(layout, figure, ff_label, label, 0, 10, [0,5,10])
        data_for_stats['approaches'][label.split('_')[0]] = entries

        print 'done with: ', label
        
        nlabels += 1
        
    #for figurename in ['approaches', 'distances', 'time', 'speeds']:
    #    #flytext.set_fontsize(layout.figures[figurename], 6)
    #    axis = layout.figures[figurename].keys()[0]
    #    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, figurename, axis, fifidatafile, 'set_fontsize ' + axis, [], 
    #                          6)
    for ax in layout.axes.values():
        ax._set_fontsize([], 6)

    #########################################################
    # STATISTICS
    #########################################################
    odor_labels = ['h2o', 'co2', 'eth', 'eth-co2', 'vinegar'] 
    for analysis_string in ['time', 'distance', 'approaches', 'speed']:
        populations = [np.array(data_for_stats[analysis_string][k]) for k in odor_labels]
        group_info = pairs2groups.label_homogeneous_groups(populations, significance_level=0.05)
        group_strings = group_info['group_strings']

        for i, odor_label in enumerate(odor_labels):
            svgtext = layout.svgitems[analysis_string][odor_label]
            svgtext.text = group_strings[i]
            svgtext.style['font-size'] = 6

        print group_info
    ########################################################

    layout.apply_svg_attrs()

    layout.append_figure_to_layer(layout.figures['approaches'], 'approaches', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['distances'], 'distances', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['time'], 'time', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['speeds'], 'speeds', cleartarget=True)
    layout.write_svg(co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_simple )

def plot_lengths_supplemental(odor, data=None):
    fifidatafile = os.path.join(os.path.dirname(co2_paper_locations.figure_template_locations.figure3_co2_calibration), 'supplemental_figure_data.pickle')
    if data is None:
        data_dict_lengths = ann.load_lengths()
        data_dict_entries = ann.load_lengths(load_data='odor_entries')
        data_dict_distances = ann.load_lengths(load_data='distances')
        data_dict_speeds = ann.load_lengths(load_data='speeds')
    else:
        data_dict_lengths, data_dict_entries, data_dict_distances, data_dict_speeds = data
    color_dict = {  'co2': 'magenta',
                    'eth': (0.001,0.001,0.001),
                    'vinegar': 'brown',
                    'h2o': 'blue'}

    labels = ['5sccm', '15sccm', '60sccm', '200sccm']

    layout = get_paper_layout('supplemental')

    yticks = []
    
    nlabels = 0

    n_trajectories = {}
    for x, label in enumerate(labels):

        #if label == 'none':
        #    continue
        #if odor + '_' not in label:
        #    continue
        #if label.startswith(odor) is False:
        #    continue
        
        ff_label = copy.copy(label)
        label = odor+'_'+ff_label
        

        if odor != 'co2' and ff_label == '5sccm':
            for m in ['speeds', 'distances', 'approaches', 'time']:
                #figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+m, ff_label)], [])
                figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, odor+'_'+m, ff_label, fifidatafile, 'adjust spines ' + odor+'_'+m + ff_label, [], 
                              [])

            continue
    
    
        MARKERSIZE = 1

        n_trajectories[ff_label] = len(np.array(data_dict_lengths[label]) )
    
        print 'Plotting: ', label, ff_label
        '''
        length = np.array(data_dict_lengths[label]) / 30. / 60.
        length = length.tolist()
        fpl.scatter_box(layout.axes[(odor+'_'+'time', ff_label)], 0, np.array(length), color=color_dict[odor], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=False, flipxy=False, xwidth=0.4)
        layout.axes[(odor+'_'+'time', ff_label)].set_xlim(-0.5,0.5)
        layout.axes[(odor+'_'+'time', ff_label)].set_ylim(0, 10)
        if nlabels == 0:
            figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+'time', ff_label)], ['left'], spine_locations={'bottom': 5, 'bottom': 5}, yticks=[0,5,10], linewidth=0.5)
            layout.axes[(odor+'_'+'time', ff_label)].tick_params(length=2.5)
        else:
            figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+'time', ff_label)], [])
        
        distances = np.array(data_dict_distances[label])*30./46
        fpl.scatter_box(layout.axes[(odor+'_'+'distances', ff_label)], 0, distances, color=color_dict[odor], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=False, flipxy=False, use='median', xwidth=0.4)
        layout.axes[(odor+'_'+'distances', ff_label)].set_xlim(-0.5,0.5)
        layout.axes[(odor+'_'+'distances', ff_label)].set_ylim(0, 100)
        if nlabels == 0:
            figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+'distances', ff_label)], ['left'], yticks=[0,50,100], spine_locations={'left': 5, 'bottom': 5}, linewidth=0.5)
            layout.axes[(odor+'_'+'distances', ff_label)].tick_params(length=2.5)
        else:
            figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+'distances', ff_label)], [])

        speeds = np.array(data_dict_speeds[label])*30./46*10 # raw speed is in pixels per frame
        fpl.scatter_box(layout.axes[(odor+'_'+'speeds', ff_label)], 0, speeds, color=color_dict[odor], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=False, flipxy=False, use='median', xwidth=0.4)
        layout.axes[(odor+'_'+'speeds', ff_label)].set_xlim(-0.5,0.5)
        layout.axes[(odor+'_'+'speeds', ff_label)].set_ylim(0, 4)
        if nlabels == 0:
            figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+'speeds', ff_label)], ['left'], spine_locations={'bottom': 5, 'bottom': 5}, yticks=[0,2,4], linewidth=0.5)
            layout.axes[(odor+'_'+'speeds', ff_label)].tick_params(length=2.5)
        else:
            figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+'speeds', ff_label)], [])
        
        entries = data_dict_entries[label]
        r = scipy.stats.uniform(-0.25,0.25)
        entries = np.array(entries) + r.rvs(len(entries))
        fpl.scatter_box(layout.axes[(odor+'_'+'approaches', ff_label)], 0, entries, color=color_dict[odor], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=False, flipxy=False, use='mean', xwidth=0.4)
        layout.axes[(odor+'_'+'approaches', ff_label)].set_xlim(-0.5,0.5)
        layout.axes[(odor+'_'+'approaches', ff_label)].set_ylim(0, 20)
        if nlabels == 0:
            figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+'approaches', ff_label)], ['left'], spine_locations={'bottom': 5, 'bottom': 5}, yticks=[0,10,20], linewidth=0.5)
            layout.axes[(odor+'_'+'approaches', ff_label)].tick_params(length=2.5)
        else:
            figurefirst.mpl_functions.adjust_spines(layout.axes[(odor+'_'+'approaches', ff_label)], [])
        '''
        def set_lims_spines(layout, figure, axis, label, ymin, ymax, yticks):
            figurefirst.deprecated_regenerate.mpl('set_xlim', layout, figure, axis, fifidatafile, 'xlim', [], 
                                      -0.5, 0.5)
            figurefirst.deprecated_regenerate.mpl('set_ylim', layout, figure, axis, fifidatafile, 'ylim', [], 
                                      ymin, ymax)

            if label == 'h2o_15sccm' or label == 'co2_5sccm' or label == 'eth_15sccm' or label == 'vinegar_15sccm':
                #figurefirst.mpl_functions.adjust_spines(layout.axes[('time', ff_label)], ['left'], spine_locations={'bottom': 5, 'bottom': 5}, yticks=[0,5,10], linewidth=0.5)
                #layout.axes[('time', ff_label)].tick_params(length=2.5)
                figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust spines', [], 
                              ['left'], yticks=yticks, spine_locations={'left': 5, 'bottom': 5}, linewidth=0.5, tick_length=2.5)
            else:
                #figurefirst.mpl_functions.adjust_spines(layout.axes[('time', ff_label)], [])
                figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust spines', [], 
                              [])

            figurefirst.deprecated_regenerate.mpl('set_rasterization_zorder', layout, figure, axis, fifidatafile, 'set_rasterization_zorder', [], 
                              1000)
            figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, figure, axis, fifidatafile, 'set_fontsize ' + axis, [], 
                              6)

        length = np.array(data_dict_lengths[label]) / 30. / 60.
        figure = label.split('_')[0] + '_time'
        axis = label.split('_')[-1]
        print figure, axis
        figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', layout, figure, axis, fifidatafile, 'stats for ' + ff_label + ' ' + figure,
                                      ['index',
                                       'List of ' + figure + ' for ' + ff_label],
                                       0, length, color=color_dict[label.split('_')[0]], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=False, flipxy=False, xwidth=0.4, use='median')
        set_lims_spines(layout, figure, ff_label, label, 0, 10, [0,5,10])

        distances = np.array(data_dict_distances[label])*30./46
        figure = label.split('_')[0] + '_distances'
        axis = label.split('_')[-1]
        #fpl.scatter_box(layout.axes[('distances', ff_label)], 0, distances, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, use='median', xwidth=0.4)
        figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', layout, figure, axis, fifidatafile, 'stats for ' + ff_label + ' ' + figure,
                                      ['index',
                                       'List of ' + figure + ' for ' + ff_label],
                                       0, distances, color=color_dict[label.split('_')[0]], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=False, flipxy=False, xwidth=0.4, use='median')
        set_lims_spines(layout, figure, ff_label, label, 0, 100, [0,50,100])

        speeds = np.array(data_dict_speeds[label])*30./46*10 # raw speed is in pixels per frame
        figure = label.split('_')[0] + '_speeds'
        axis = label.split('_')[-1]
        #fpl.scatter_box(layout.axes[('speeds', ff_label)], 0, speeds, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, use='median', xwidth=0.4)
        figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', layout, figure, axis, fifidatafile, 'stats for ' + ff_label + ' ' + figure,
                                      ['index',
                                       'List of ' + figure + ' for ' + ff_label],
                                       0, speeds, color=color_dict[label.split('_')[0]], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=False, flipxy=False, xwidth=0.4, use='median')
        set_lims_spines(layout, figure, ff_label, label, 0, 4, [0,2,4])
        
        entries = data_dict_entries[label]
        r = scipy.stats.uniform(-0.25,0.25)
        entries = np.array(entries) + r.rvs(len(entries))
        figure = label.split('_')[0] + '_approaches'
        axis = label.split('_')[-1]
        #fpl.scatter_box(layout.axes[('approaches', ff_label)], 0, entries, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=True, flipxy=False, use='mean', xwidth=0.4)
        figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_box', layout, figure, axis, fifidatafile, 'stats for ' + ff_label + ' ' + figure,
                                      ['index',
                                       'List of ' + figure + ' for ' + ff_label],
                                       0, entries, color=color_dict[label.split('_')[0]], shading='95conf', markersize=MARKERSIZE, edgecolor='none', hide_markers=False, flipxy=False, xwidth=0.4, use='mean')
        set_lims_spines(layout, figure, ff_label, label, 0, 20, [0,10,20])
        print 'done with: ', label
        
        nlabels += 1
        
    #for figurename in [odor+'_'+'approaches', odor+'_'+'distances', odor+'_'+'time', odor+'_'+'speeds']:
    #    axis = layout.figures[figurename].keys()[0]
    #    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, figurename, axis, fifidatafile, 'set_fontsize ' + axis, [], 
    #                          6)

    layout.append_figure_to_layer(layout.figures[odor+'_'+'approaches'], odor+'_'+'approaches', cleartarget=True)
    layout.append_figure_to_layer(layout.figures[odor+'_'+'distances'], odor+'_'+'distances', cleartarget=True)
    layout.append_figure_to_layer(layout.figures[odor+'_'+'time'], odor+'_'+'time', cleartarget=True)
    layout.append_figure_to_layer(layout.figures[odor+'_'+'speeds'], odor+'_'+'speeds', cleartarget=True)
    layout.write_svg(co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_supp )

    print
    print
    print '****** ', odor, ' ******'
    print n_trajectories
    print '************************'

def plot_lengths(data=None):
    if data is None:
        data_dict_lengths = ann.load_lengths()
        data_dict_entries = ann.load_lengths(load_data='odor_entries')
        data_dict_distances = ann.load_lengths(load_data='distances')
        data_dict_speeds = ann.load_lengths(load_data='speeds')
    else:
        data_dict_lengths, data_dict_entries, data_dict_distances, data_dict_speeds = data
    color_dict = {  'orco_co260': 'orange',
                    'co2_5sccm': 'magenta',
                    'co2_15sccm': 'magenta',
                    'co2_60sccm': 'magenta',
                    'co2_hungry': 'magenta',
                    'eth_60sccm': (0.001,0.001,0.001),
                    'co2_200sccm': 'magenta',
                    'eth_200sccm': (0.001,0.001,0.001),
                    'eth_15sccm': (0.001,0.001,0.001),
                    'eth-co2_60_15': 'purple',
                    'vinegar_15sccm': 'brown',
                    'vinegar_60sccm': 'brown',
                    'vinegar_200sccm': 'brown',
                    'vinegar_60sccm_15co2': 'purple',
                    'h2o_15sccm': 'blue',
                    'h2o_60sccm': 'blue',
                    'h2o_200sccm': 'blue',
                    'ethylacetate_15_50': 'green',
                    'ethylacetate_15sccm': 'green',
                    'ethylacetate_60sccm': 'green',
                    'choice-eth-only_60sccm': (0.001, 0.001, 0.001),
                    'choice-eth-and-co2_60sccm_15sccm': 'purple'}
     
    labels = ['h2o_15sccm', 'h2o_60sccm', 'h2o_200sccm', 'co2_5sccm', 'co2_15sccm', 'co2_60sccm', 'co2_200sccm', 'eth_15sccm', 'eth_60sccm', 'eth_200sccm', 'eth-co2_60_15', 'vinegar_15sccm', 'vinegar_60sccm', 'vinegar_200sccm', 'ethylacetate_15_50', 'ethylacetate_15sccm', 'ethylacetate_60sccm', 'co2_15sccm_2', 'eth_60sccm_2']
    #labels = ['co2_5sccm']
    #labels = ['h2o_60sccm', 'co2_60sccm', 'eth_60sccm', 'eth-co2_60_15', 'vinegar_60sccm', 'ethylacetate_15sccm', 'choice-eth-only_60sccm', 'choice-eth-and-co2_60sccm_15sccm']
    #labels = ['h2o_60sccm', 'co2_60sccm', 'eth_60sccm', 'eth-co2_60_15', 'vinegar_60sccm', 'ethylacetate_15sccm']
    labels.reverse()
    
    #'vinegar_60sccm_15co2',  'h2o_60sccm', 'orco_co260', 

    layout = get_paper_layout()
    #layout_figure = 'pad_' + odor
    #try:
    #    layout.clear_fflayer(layout_figure)
    #except:
    #    pass

    yticks = []
    
    nlabels = 0
    for x, label in enumerate(labels):
    
        #if label == 'none':
        #    continue
        #if odor + '_' not in label:
        #    continue
        #if label.startswith(odor) is False:
        #    continue
        
        if label[-2:] == '_2':
            ff_label = copy.copy(label)
            label = label[0:-2]
        else:
            ff_label = label
    
        MARKERSIZE = 2
    
        print 'Plotting: ', label, ff_label
        length = np.array(data_dict_lengths[label]) / 30. / 60.
        length = length.tolist()
        fpl.scatter_box(layout.axes[('time_on_pad', ff_label)], 0, np.array(length), color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', flipxy=True, xwidth=0.4)
        layout.axes[('time_on_pad', ff_label)].set_ylim(-0.5,0.5)
        layout.axes[('time_on_pad', ff_label)].set_xlim(0, 10)
        if label == 'ethylacetate_60sccm':
            figurefirst.mpl_functions.adjust_spines(layout.axes[('time_on_pad', ff_label)], ['bottom'], spine_locations={'bottom': 5, 'bottom': 5}, xticks=[0,10])
        else:
            figurefirst.mpl_functions.adjust_spines(layout.axes[('time_on_pad', ff_label)], [])
        
        distances = np.array(data_dict_distances[label])*30./46
        fpl.scatter_box(layout.axes[('distances', ff_label)], 0, distances, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', flipxy=True, use='median', xwidth=0.4)
        layout.axes[('distances', ff_label)].set_ylim(-0.5,0.5)
        layout.axes[('distances', ff_label)].set_xlim(0, 100)
        if label == 'ethylacetate_60sccm':
            figurefirst.mpl_functions.adjust_spines(layout.axes[('distances', ff_label)], ['bottom'], spine_locations={'bottom': 5, 'bottom': 5}, xticks=[0,100])
        else:
            figurefirst.mpl_functions.adjust_spines(layout.axes[('distances', ff_label)], [])
        # single pad: 47.2 pixels / cm
        # left and right: 45.1 pixels / cm

        speeds = np.array(data_dict_speeds[label])*30./46*10 # raw speed is in pixels per frame
        fpl.scatter_box(layout.axes[('speeds', ff_label)], 0, speeds, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', flipxy=True, use='median', xwidth=0.4)
        layout.axes[('speeds', ff_label)].set_ylim(-0.5,0.5)
        layout.axes[('speeds', ff_label)].set_xlim(0, 4)
        if label == 'ethylacetate_60sccm':
            figurefirst.mpl_functions.adjust_spines(layout.axes[('speeds', ff_label)], ['bottom'], spine_locations={'bottom': 5, 'bottom': 5}, xticks=[0,4])
        else:
            figurefirst.mpl_functions.adjust_spines(layout.axes[('speeds', ff_label)], [])
        # single pad: 47.2 pixels / cm
        # left and right: 45.1 pixels / cm
        
        entries = data_dict_entries[label]
        r = scipy.stats.uniform(-0.25,0.25)
        entries = np.array(entries) + r.rvs(len(entries))
        fpl.scatter_box(layout.axes[('approaches', ff_label)], 0, entries, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', flipxy=True, use='mean', xwidth=0.4)
        layout.axes[('approaches', ff_label)].set_ylim(-0.5,0.5)
        layout.axes[('approaches', ff_label)].set_xlim(0, 20)
        if label == 'ethylacetate_60sccm':
            figurefirst.mpl_functions.adjust_spines(layout.axes[('approaches', ff_label)], ['bottom'], spine_locations={'bottom': 5, 'bottom': 5}, xticks=[0,20])
        else:
            figurefirst.mpl_functions.adjust_spines(layout.axes[('approaches', ff_label)], [])

        if 1:
            directory = co2_paper_locations.data_locations.windtunnel_walking[label]
            fraction_of_time_spent_walking, mean_walking_speed = get_fraction_of_time_spent_walking_and_mean_walking_speed(directory)
            r = scipy.stats.uniform(-0.25,0.25)
            fraction_of_time_spent_walking = np.array(fraction_of_time_spent_walking) + r.rvs(len(fraction_of_time_spent_walking))
            fpl.scatter_box(layout.axes[('fraction_of_time_spent_walking', ff_label)], 0, fraction_of_time_spent_walking, color=color_dict[label], shading='95conf', markersize=MARKERSIZE, edgecolor='none', flipxy=True, use='mean', xwidth=0.4)
            layout.axes[('fraction_of_time_spent_walking', ff_label)].set_ylim(-0.5,0.5)
            layout.axes[('fraction_of_time_spent_walking', ff_label)].set_xlim(0, 1)
            if label == 'ethylacetate_60sccm':
                figurefirst.mpl_functions.adjust_spines(layout.axes[('fraction_of_time_spent_walking', ff_label)], ['bottom'], spine_locations={'bottom': 5, 'bottom': 5}, xticks=[0,1])
            else:
                figurefirst.mpl_functions.adjust_spines(layout.axes[('fraction_of_time_spent_walking', ff_label)], [])
        
        print 'done with: ', label
        
        nlabels += 1
        
    for figurename in ['approaches', 'distances', 'time_on_pad', 'speeds', 'fraction_of_time_spent_walking']:
        flytext.set_fontsize(layout.figures[figurename], 6)

    layout.append_figure_to_layer(layout.figures['approaches'], 'approaches', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['distances'], 'distances', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['time_on_pad'], 'time_on_pad', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['speeds'], 'speeds', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['fraction_of_time_spent_walking'], 'fraction_of_time_spent_walking', cleartarget=True)
    layout.write_svg(co2_paper_locations.figure_template_locations.figure3_windtunnel_walking )
    

def get_fraction_of_time_spent_walking_and_mean_walking_speed(directory):
    paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='day')
    fraction_of_time_spent_walking = []
    mean_walking_speed = []

    for path in paths:
        dataset = mta.read_hdf5_file_to_pandas.load_dataset_from_path(path, load_saved=True)
        dataset.calculate_function_for_all_trajecs(mta.trajectory_analysis.calculate_fraction_of_time_spent_walking)
        dataset.calculate_function_for_all_trajecs(mta.trajectory_analysis.calculate_mean_walking_speed_while_walking)
        
        f = [dataset.trajecs[key].fraction_of_time_spent_walking for key in dataset.keys]
        fraction_of_time_spent_walking.extend(f)

        s = [dataset.trajecs[key].mean_walking_speed for key in dataset.keys]
        mean_walking_speed.extend(s)

    return fraction_of_time_spent_walking, mean_walking_speed

    
if __name__ == '__main__':
    plot_all_lengths()
