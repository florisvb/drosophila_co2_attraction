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

import plot_nflies_vs_speed_scatter
import co2_paper_locations

FLOWRATES = 1

def get_mean_peak_attraction_and_aversion(directories, flowrates=FLOWRATES, localtimerange=[14,36], average_within_paths=False, use_speed_intercept=2.3):
    minimum_speed = 0

    response, speed, indices_to_sort, sp, re, experiment_type, t = plot_nflies_vs_speed_scatter.get_new_attraction_index(directories, flowrates, localtimerange, average_within_paths)

    x_intercept = use_speed_intercept
    x_intercept_index = np.argmin(np.abs(sp-x_intercept))
    x_intercept_index_minspeed = np.argmin(np.abs(sp-minimum_speed))

    high_speed_response = np.mean(response[x_intercept_index:,:], axis=0)
    low_speed_response = np.mean(response[x_intercept_index_minspeed:x_intercept_index,:], axis=0)

    responses = np.hstack((high_speed_response, low_speed_response))

    attraction = np.max(high_speed_response[2400:4800])
    aversion = np.min(low_speed_response[2400:4800])

    if 'single' not in ''.join(directories):
        return attraction/10., aversion/10.
    else:
        return attraction, aversion

def get_distribution_peak_attraction_and_aversion(directories, flowrates=FLOWRATES, localtimerange=[14,36], average_within_paths=False, use_speed_intercept=2.3):
    minimum_speed = 0

    response, speed, indices_to_sort, sp, re, experiment_type, t = plot_nflies_vs_speed_scatter.get_new_attraction_index(directories, flowrates, localtimerange, average_within_paths)

    x_intercept = use_speed_intercept
    x_intercept_index = np.argmin(np.abs(sp-x_intercept))
    x_intercept_index_minspeed = np.argmin(np.abs(sp-minimum_speed))

    high_speed_response = response[x_intercept_index:,:]
    low_speed_response = response[x_intercept_index_minspeed:x_intercept_index,:]

    attraction = np.max(high_speed_response[:,2400:4800], axis=1)
    aversion = np.min(low_speed_response[:,2400:4800], axis=1)


    # means
    mean_high_speed_response = np.mean(response[x_intercept_index:,:], axis=0)
    mean_low_speed_response = np.mean(response[x_intercept_index_minspeed:x_intercept_index,:], axis=0)
    mean_attraction = np.max(mean_high_speed_response[2400:4800])
    mean_aversion = np.min(mean_low_speed_response[2400:4800])

    if 'single' not in ''.join(directories):
        return attraction/10., aversion/10., mean_attraction/10., mean_aversion/10.
    else:
        return attraction, aversion, mean_attraction, mean_aversion

def compare_labels(label1, label2):
    def get_data(label):
        directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics, label)
        if label in ['ir25a2_BAC', 'ir40a', 'ir21a', 'ir93a']:
            DNC = ['hot']
        else:
            DNC = []
        directories = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='', does_not_contain=DNC)
        attraction, aversion, mean_attraction, mean_aversion = get_distribution_peak_attraction_and_aversion(directories)
        return attraction, aversion, mean_attraction, mean_aversion

    attraction_1, aversion_1, mean_attraction_1, mean_aversion_1 = get_data(label1)
    attraction_2, aversion_2, mean_attraction_2, mean_aversion_2 = get_data(label2)

    ks_attraction_stat, ks_attraction_pval = scipy.stats.ks_2samp(attraction_1, attraction_2)
    ks_aversion_stat, ks_aversion_pval = scipy.stats.ks_2samp(aversion_1, aversion_2)

    return ks_attraction_stat, ks_attraction_pval, ks_aversion_stat, ks_aversion_pval, mean_attraction_2, mean_aversion_2, attraction_2, aversion_2

def save_statistical_data(layout):
    labels = ['hcs', 'ir64a', 'gr63a', 'double_gr63_ir64', 'anosmic',  
              'antennaless', 'orco', 'ir8aM120', 'M37ir25a2', 
              'ir25a2_BAC', 'M106', 'ir40a']

    pretty_label_name = {'hcs': 'hcs', 
                         'ir64a': 'ir64a', 
                         'gr63a': 'gr63a',
                         'double_gr63_ir64': 'gr63a_ir64a', 
                         'anosmic': 'anosmic',  
                         'antennaless': 'antennaless',
                         'orco': 'orco',
                         'ir8aM120': 'ir8a',
                         'M37ir25a2': 'ir25a', 
                         'ir25a2_BAC': 'ir25a_BAC',
                         'M106': 'orco_ir8a_gr63a', 
                         'ir40a': 'ir40a'}

    attraction_values = {}
    aversion_values = {}
    attraction_stat = {}
    aversion_stat = {}
    attraction_pval = {}
    aversion_pval = {}
    raw_attraction_values = {}
    raw_aversion_values = {}

    for label in labels:
        ks_attraction_stat, ks_attraction_pval, ks_aversion_stat, ks_aversion_pval, mean_attraction_2, mean_aversion_2, attraction_2, aversion_2 = compare_labels('hcs', label)
        attraction_values[label] = mean_attraction_2
        aversion_values[label] = mean_aversion_2
        attraction_stat[label] = ks_attraction_stat
        aversion_stat[label] = ks_aversion_stat
        attraction_pval[label] = ks_attraction_pval
        aversion_pval[label] = ks_aversion_pval
        raw_attraction_values[label] = attraction_2
        raw_aversion_values[label] = aversion_2

    layout.write_fifidata(['Summary Statistics for:\n###### hcs, ir64a, gr63a, gr63_ir64, anosmic, antennaless, orco, ir8a, ir25a, ir25a_BAC, orco_ir8a_gr63a, ir40a', 
                           'N trials Attraction',
                           'N trials Aversion',
                           'Mean Max Attraction Index',
                           'Mean Min Aversion Index',
                           'K-S test statistic compared to HCS for Attraction Index',
                           'K-S p-value compared to HCS for Attraction Index',
                           'K-S test statistic compared to HCS for Aversion Index',
                           'K-S p-value compared to HCS for Aversion Index'],
                           [len(raw_attraction_values[label]) for label in labels],
                           [len(raw_aversion_values[label]) for label in labels],
                           [attraction_values[label] for label in labels],
                           [aversion_values[label] for label in labels],
                           [attraction_stat[label] for label in labels],
                           [attraction_pval[label] for label in labels],
                           [aversion_stat[label] for label in labels],
                           [aversion_pval[label] for label in labels],
                           [aversion_values[label] for label in labels],
                           )

    
    for label in labels:
        layout.write_fifidata(['Complete data for: ' + pretty_label_name[label], 
                               'Attraction Index',
                               'Aversion Index'],
                               raw_attraction_values[label],
                               raw_aversion_values[label],
                               )

def write_in_stats(odor='co2', bonferoni_N=11):

    labels = ['hcs', 'ir64a', 'gr63a', 'double_gr63_ir64', 'anosmic',  
              'antennaless', 'orco', 'ir8aM120', 'M37ir25a2', 
              'ir25a2_BAC', 'M106', 'ir40a']#, 'ir21a', 'ir93aM98']

    if odor == 'co2':
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2
    else:
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2_control
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()


    svgitems = []
    meanvals = []
    p_level_1 = 0.05/float(bonferoni_N)
    p_level_2 = 0.01/float(bonferoni_N)
    p_level_3 = 0.005/float(bonferoni_N)

    for label in labels:
        if label == 'hcs':
            continue
        ks_attraction_stat, ks_attraction_pval, ks_aversion_stat, ks_aversion_pval, mean_attraction_2, mean_aversion_2, attraction_2, aversion_2 = compare_labels('hcs', label)

        def write_text(ks_attraction_pval, mean_attraction_label, suffix, svgitems, meanvals):
            if suffix == '_attr':
                svg_group_name = 'attraction_pval'
            else:
                svg_group_name = 'aversion_pval'
            p_level = 0
            if ks_attraction_pval < p_level_1:
                p_level = 1
            if ks_attraction_pval < p_level_2:
                p_level = 2
            if ks_attraction_pval < p_level_3:
                p_level = 3
            stars = '*'*p_level
            if ks_attraction_pval < 0.001:
                text = '<.001'+stars
            else:
                text = ''+"%.3f"%ks_attraction_pval+''+stars
            layout.svgitems[svg_group_name][label+'_attr'].text = text
            layout.svgitems[svg_group_name][label+'_attr'].style['font-size'] = 768
            layout.svgitems[svg_group_name][label+'_attr'].style['font-weight'] = 'bold'
            layout.svgitems[svg_group_name][label+'_attr'].style['font-family'] = 'Arial'
            layout.svgitems[svg_group_name][label+'_attr'].style['-inkscape-font-specification'] = 'Arial, Bold'
            
            #if ks_attraction_pval > 0.005 or np.abs(mean_attraction_label) > 0.2:
            #    layout.svgitems[('attraction_pval', label+suffix)].style['fill'] = '#ffffff'
            #else:
            layout.svgitems[svg_group_name][label+'_attr'].style['fill'] = '#000000'

            #svgitems.append(label+suffix)
            svgitems.append(svg_group_name)
            meanvals.append(np.abs(mean_attraction_label))
            return svgitems, meanvals
    
        svgitems, meanvals = write_text(ks_attraction_pval, mean_attraction_2, '_attr', svgitems, meanvals)
        svgitems, meanvals = write_text(ks_aversion_pval, mean_aversion_2, '_aver', svgitems, meanvals)
    
    print
    print svgitems
    print meanvals

    layout.apply_svg_attrs(svg_items_to_update=svgitems)
    layout.write_svg(svg)

def write_in_values(odor='co2'):

    labels = ['hcs', 'ir64a', 'gr63a', 'double_gr63_ir64', 'anosmic',  
              'antennaless', 'orco', 'ir8aM120', 'M37ir25a2', 
              'ir25a2_BAC', 'M106', 'ir40a']#, 'ir21a', 'ir93aM98']

    if odor == 'co2':
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2
    else:
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2_control
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()


    svgitems = []
    meanvals = []

    for label in labels:
        ks_attraction_stat, ks_attraction_pval, ks_aversion_stat, ks_aversion_pval, mean_attraction_2, mean_aversion_2, attraction_2, aversion_2 = compare_labels('hcs', label)

        def write_text(ks_attraction_pval, mean_attraction_label, suffix, svgitems, meanvals):
            if suffix == '_attr':
                svg_group_name = 'attraction'
            else:
                svg_group_name = 'aversion'
            text = ''+"%.3f"%mean_attraction_label+''
            layout.svgitems[svg_group_name][label+'_attr'].text = text
            layout.svgitems[svg_group_name][label+'_attr'].style['font-size'] = 768
            layout.svgitems[svg_group_name][label+'_attr'].style['font-weight'] = 'bold'
            layout.svgitems[svg_group_name][label+'_attr'].style['font-family'] = 'Arial'
            layout.svgitems[svg_group_name][label+'_attr'].style['-inkscape-font-specification'] = 'Arial, Bold'
            if np.abs(mean_attraction_label) > 0.1:
                layout.svgitems[svg_group_name][label+'_attr'].style['fill'] = '#ffffff'
            else:
                layout.svgitems[svg_group_name][label+'_attr'].style['fill'] = '#000000'

            #svgitems.append(label+suffix)
            svgitems.append(svg_group_name)
            meanvals.append(np.abs(mean_attraction_label))
            return svgitems, meanvals
    
        svgitems, meanvals = write_text(ks_attraction_pval, mean_attraction_2, '_attr', svgitems, meanvals)
        svgitems, meanvals = write_text(ks_aversion_pval, mean_aversion_2, '_aver', svgitems, meanvals)
    
    print
    print svgitems
    print meanvals

    layout.apply_svg_attrs(svg_items_to_update=svgitems)
    layout.write_svg(svg)

def get_attraction_aversion_array():

    labels = ['hcs', 'ir64a', 'gr63a', 'double_gr63_ir64', 'anosmic',  
              'antennaless', 'orco', 'ir8aM120', 'M37ir25a2', 
              'ir25a2_BAC', 'M106', 'ir40a']#, 'ir21a', 'ir93aM98']
    attraction_aversion_array = np.zeros([2, len(labels)])

    for i, label in enumerate(labels):
        if label is None:
            continue

        directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics, label)
        if label in ['ir25a2_BAC', 'ir40a', 'ir21a']:
            DNC = ['hot']
        else:
            DNC = []

        directories = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='', does_not_contain=DNC)

        attraction, aversion = get_mean_peak_attraction_and_aversion(directories)

        attraction_aversion_array[0, i] = attraction
        attraction_aversion_array[1, i] = aversion

    return attraction_aversion_array

def plot_attraction_aversion_array_on_layout(odor='co2', attraction_aversion_array=None):
    if attraction_aversion_array is None:
        attraction_aversion_array = get_attraction_aversion_array()

    if odor == 'co2':
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2
    else:
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2_control
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()

    # add white rows
    attr_aver_array_with_pvals = np.zeros([4,attraction_aversion_array.shape[1]])
    attr_aver_array_with_pvals[0,:] = attraction_aversion_array[0,:]
    attr_aver_array_with_pvals[2,:] = attraction_aversion_array[1,:]

    layout.axes['summary_array','summary_array']._imshow(['Summary statistics', 'Mean Attraction Indices and P-value'], attr_aver_array_with_pvals, cmap='seismic', vmin=-0.7, vmax=0.7, interpolation='nearest', aspect='auto')
    
    ax = layout.axes['summary_array','summary_array']
    ax.record = True
    ax.set_xticks(np.linspace(-0.5,11.5, 13))
    ax.set_xticklabels([])
    ax.set_yticks([-0.5, 0.5, 1.5, 2.5, 3.5])
    ax.set_yticklabels([])
    ax.grid(color='black', linewidth=0.5, linestyle='-')
    ax.record = False

    def set_spine_widths(ax):
        for loc, spine in ax.spines.items():
            spine.set_linewidth(0.5)
    ax._custom([], set_spine_widths)

    layout.append_figure_to_layer(layout.figures['summary_array'], 'summary_array', cleartarget=True)
    layout.write_svg(svg)

    #write_in_stats(bonferoni_N=13)

    return attraction_aversion_array

if __name__ == '__main__':
    odor = 'co2'
    plot_attraction_aversion_array_on_layout(odor=odor)
    write_in_stats(odor=odor, bonferoni_N=11)
    write_in_values(odor=odor)
