from multicat_analysis import plot_n_flies
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
import copy

from multicat_analysis import plot_nflies_vs_speed_scatter

import co2_paper_locations

USE_SPEED_INTERCEPT = 2.3
USE_SPEED_INTERCEPT_CONCENTRATION = 2.5
USE_ETH_SPEED_INTERCEPT = 1

DOES_NOT_CONTAIN=['~', '.pyc', 'day','.svg','tmp']

def make_cross_directory_figure(directory, label, flowrate=1, use_speed_intercept=None, figure='co2'):
    speed_svg = None

    if 'control' in figure:
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_control
        speed_svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_speed_control
    elif figure == 'co2' or figure == 'ethanol' or figure == 'vinegar':
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only
        speed_svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_speed

    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()

    if speed_svg is not None:
        speed_layout = figurefirst.svg_to_axes.FigureLayout(speed_svg)
        speed_layout.make_mplfigures()
    else:
        speed_layout = None

    odor = 'co2'
    if 'ethanol' in figure or 'eth' in figure:
        odor = 'ethanol'
    if 'vinegar' in figure:
        odor = 'vinegar'

    DNC = copy.copy(DOES_NOT_CONTAIN)
    if '_nothot' in label:
        DNC.append('hot')
    if '_hot' in label:
        DNC.append('24hrs')
        DNC.append('noon')
    if label == 'ir25a2_BAC':
        DNC.append('hot')
    directories = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='', does_not_contain=DNC)
    plot_nflies_vs_speed_scatter.plot_odor_response_speed_and_time_matrix(directories, flowrates=flowrate, localtimerange=[14,36], 
        average_within_paths=False, layout=layout, flowrate_label=label, use_speed_intercept=use_speed_intercept, odor=odor, speed_layout=None)

    #layout.apply_svg_attrs()
    layout.append_figure_to_layer(layout.figures[label], label, cleartarget=True)
    layout.write_svg(svg)

    if 0: #speed_layout is not None:
        speed_layout.append_figure_to_layer(speed_layout.figures[label], label, cleartarget=True)
        speed_layout.write_svg(speed_svg)

def make_summary_of_genetic_experiments(figure='co2'):

    if 'control' in figure:
        flowrate = 0
    if figure == 'co2':
        flowrate = 1
    elif figure == 'ethanol' or figure == 'vinegar':
        flowrate = 3

    if 'co2' in figure:
        labels = ['hcs', 'anosmic', 'orco', 'M37ir25a2', 'antennaless', 'M106', 'ir25a2_BAC', 'ir8aM120']
        labels_supp = ['hcs_00', 'hcs_04', 'hcs_12', 'hcs_35', 'gr63a', 'ir64a', 'double_gr63_ir64', 'single_flies_h_c_s']
        labels_ir25a = ['ir25a_nothot', 'ir25a_hot', 'M37ir25a2_nothot', 'M37ir25a2_hot', 'ir25a2_BAC_nothot', 'ir25a2_BAC_hot']
        labels_ir21a_ir93a = ['ir21a', 'ir93aM98']
        labels.extend(labels_supp)
        labels.extend(labels_ir25a)
        labels.extend(labels_ir21a_ir93a)



    elif 'ethanol' in figure:
        labels = ['eth_hcs', 'eth_anosmic', 'eth_orco', 'eth_M37ir25a2', 'eth_antennaless', 'eth_M106', 'eth_ir25a2_BAC', 'eth_ir21a', 'eth_ir8aM120']#, 'eth_ir93aM98']
    elif 'vinegar' in figure:
        labels = ['vinegar_hcs', 'vinegar_anosmic', 'vinegar_orco', 'vinegar_M37ir25a2', 'vinegar_M106', 'vinegar_ir8aM120']#, 'eth_antennaless', 'eth_M106', 'eth_ir25a2_BAC']

    for label in labels:
        if 'ethanol' in figure or 'eth' in figure:
            directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics_ethanol, label)
            make_cross_directory_figure(directory, 
                label, flowrate=flowrate, use_speed_intercept=USE_ETH_SPEED_INTERCEPT, figure=figure)
        elif 'vinegar' in figure:
            directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics_vinegar, label)
            make_cross_directory_figure(directory, 
                label, flowrate=flowrate, use_speed_intercept=USE_ETH_SPEED_INTERCEPT, figure=figure)
        elif 'hcs_' in label:
            directory = co2_paper_locations.data_locations.walking_arena_hcs_concentration
            if 'control' not in figure:
                f = float(label.split('_')[-1])/10.
            else:
                f = 0
            make_cross_directory_figure(directory, 
                label, flowrate=f, use_speed_intercept=USE_SPEED_INTERCEPT, figure=figure)
        else:
            if '_hot' in label:
                l = label.split('_hot')[0]
                directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics, l)  
            elif '_nothot' in label:
                l = label.split('_nothot')[0]
                directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics, l)    
            else:
                directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics, label)
            make_cross_directory_figure(directory, 
                label, flowrate=flowrate, use_speed_intercept=USE_SPEED_INTERCEPT, figure=figure)

def make_supplemental_summary_of_genetic_experiments(figure='supplemental'):
    labels = ['hcs', 'gr21a_x_tnt', 'parental_control_tnt']

    if figure == 'experiment':
        flowrate = 1
    elif figure == 'control':
        flowrate = 0
    elif figure == 'supplemental':
        flowrate = 1

    for label in labels:
        directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics, label)
        make_cross_directory_figure(directory, 
            label, flowrate=flowrate, use_speed_intercept=USE_SPEED_INTERCEPT, figure=figure)

def make_colormaps(svg):
    
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()

    #fpl.colorbar(ax=layout.axes[('cmap', 'speed_cmap')], colormap='hot')
    #fpl.colorbar(ax=layout.axes[('cmap', 'pref_ind_cmap')], colormap='bwr')

    #if 'only' in svg:
    #    fpl.colorbar(ax=layout.axes[('cmap_h', 'speed_cmap_h')], colormap='hot', orientation='horizontal')
    #    fpl.colorbar(ax=layout.axes[('cmap_h', 'pref_ind_cmap_h')], colormap='seismic', orientation='vertical')
    #else:

    fpl.colorbar(ax=layout.axes[('cmap_h', 'speed_cmap_h')], colormap='hot', orientation='horizontal')
    fpl.colorbar(ax=layout.axes[('cmap_h', 'pref_ind_cmap_h')], colormap='bwr', orientation='horizontal')

    #figurefirst.mpl_functions.adjust_spines(layout.axes[('cmap_h', 'speed_cmap_h')], [])
    #figurefirst.mpl_functions.adjust_spines(layout.axes[('cmap_h', 'pref_ind_cmap_h')], [])
    #layout.append_figure_to_layer(layout.figures['cmap'], 'cmap', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['cmap_h'], 'cmap_h', cleartarget=True)
    layout.write_svg(svg)

def make_hcs_full(concentration_dataset=False, odor='co2', figure='co2'):
    if 'control' in figure:
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_control
    elif 'co2' in figure:
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only
    elif 'ethanol' in figure:
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_hcs_full
    elif figure == 'gr63a':
        pass#svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_gr63a
    elif figure == 'concentration':
        pass#svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_concentration
    elif 'ir21a' in figure:
        pass#svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_ir21a

    print svg

    if odor == 'co2':
        if concentration_dataset:
            paths = mta.read_hdf5_file_to_pandas.get_filenames(co2_paper_locations.data_locations.walking_arena_hcs_concentration, 
                contains='', does_not_contain=DOES_NOT_CONTAIN)
            flowrate = 1.2
        else:
            paths = mta.read_hdf5_file_to_pandas.get_filenames(co2_paper_locations.data_locations.walking_arena_tripleparadigm_hcs, 
                contains='', does_not_contain=DOES_NOT_CONTAIN)
            flowrate = 1
    elif odor == 'ethanol':
        paths = mta.read_hdf5_file_to_pandas.get_filenames(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics_ethanol_hcs, 
                contains='', does_not_contain=DOES_NOT_CONTAIN)
        flowrate = 3
    elif odor == 'vinegar':
        paths = mta.read_hdf5_file_to_pandas.get_filenames(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics_vinegar_hcs, 
                contains='', does_not_contain=DOES_NOT_CONTAIN)
        flowrate = 3


    use_speed_intercept = None
    if 'control' in figure:
        flowrate = 0
        if odor == 'co2':
            use_speed_intercept = USE_SPEED_INTERCEPT

    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()

    if odor == 'co2':
        plot_nflies_vs_speed_scatter.plot_odor_response_speed_and_time_matrix(paths, flowrates=flowrate, flowrate_label='hcs_full', 
                                                                                use_speed_intercept=use_speed_intercept, layout=layout)
        layout.append_figure_to_layer(layout.figures['hcs_full'], 'hcs_full', cleartarget=True)
    elif odor == 'ethanol':
        plot_nflies_vs_speed_scatter.plot_odor_response_speed_and_time_matrix(paths, flowrates=flowrate, flowrate_label='eth_hcs_full', 
                                                                                use_speed_intercept=1, layout=layout)
        layout.append_figure_to_layer(layout.figures['eth_hcs_full'], 'eth_hcs_full', cleartarget=True)
    elif odor == 'vinegar':
        plot_nflies_vs_speed_scatter.plot_odor_response_speed_and_time_matrix(paths, flowrates=flowrate, flowrate_label='vinegar_hcs_full', 
                                                                                use_speed_intercept=1, layout=layout)
        layout.append_figure_to_layer(layout.figures['vinegar_hcs_full'], 'vinegar_hcs_full', cleartarget=True)
    
    layout.write_svg(svg)

    make_colormaps(svg)

def make_figure():
    make_hcs_full()
    make_summary_of_concentration_experiments()
    make_summary_of_genetic_experiments()

    #make_colormaps()

def make_control_figure():
    make_summary_of_genetic_experiments(figure='control')

def make_hcs_full_for_specific_directory(directory, flowrate, name, use_speed_intercept):
    svg_template = co2_paper_locations.figure_template_locations.figure5_full_template
    svg_output = os.path.join(os.path.dirname(svg_template), name)

    paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, 
        contains='', does_not_contain=DOES_NOT_CONTAIN)


    print paths

    layout = figurefirst.svg_to_axes.FigureLayout(svg_template)
    layout.make_mplfigures()

    plot_nflies_vs_speed_scatter.plot_odor_response_speed_and_time_matrix(paths, flowrates=flowrate, flowrate_label='hcs_full', 
                                                                            use_speed_intercept=use_speed_intercept, layout=layout)
    layout.append_figure_to_layer(layout.figures['hcs_full'], 'hcs_full', cleartarget=True)
    
    layout.write_svg(svg_output)

    