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

import figurefirst
import data_fit
import matplotlib.patches

import orchard.annotated_time_on_pad_analysis as ann

import co2_paper_locations


def get_paper_layout(version='old'):
    if version == 'old':
        svg = co2_paper_locations.figure_template_locations.figure3_windtunnel_walking
    else:
        svg = co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_simple
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout


def get_filename(path, contains):
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
            return os.path.join(path, filename)

def plot_trajectories(ax, paths, keys, key_lengths=None, downsample=10):
    # dead flies in 200 sccm co2 day 4: 5932, 4394; use lengths: {5932: 9000, 4394: 3000}
    # ethanol trajectory: day2, trajectory 13437 '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_ethanol/day2'
    # h20 trajectory: day2, trajec 299  '/media/Orchard2/singlepad_windtunnel_blacktop/h20_60sccm/day2'
    if type(paths) is not list:
        paths = [paths]
        keys = [keys]
        key_lengths = [key_lengths]
        
    path = paths[0]
    bgimg_filename = get_filename(path, 'bgimg_N1.png')
    bgimg = plt.imread(bgimg_filename)

    # equalize bgimgs to match
    r = 0.34/np.mean(bgimg)
    bgimg *= r
    

    ax.set_aspect('equal')
    ax.set_xlim(150,630)

    #ax.imshow(bgimg, cmap='gray', zorder=1)
    mids = 1- np.abs(bgimg - 0.5)
    darks = 1- np.abs(0.2-bgimg) 
    #ax.imshow(bgimg + mids*2 - darks*2, cmap='gray', zorder=1)

    for p, path in enumerate(paths):
        dataset_filename = os.path.join(path, 'dataset.pickle')
        if not os.path.exists(dataset_filename):
            ann.save_annotated_dataset(path)
        f = open(dataset_filename)
        dataset = pickle.load(f)
        f.close()
        
        for key in keys[p]: 
            trajec = dataset.trajec(key)
            
            if key_lengths is not None:
                l = key_lengths[p][key]
            else:
                l = len(trajec.speed)
            
            interpolated_indices = np.where(trajec.interpolated==1)[0]
            print len(interpolated_indices)
            r = np.arange(0, len(interpolated_indices), 5)
            interpolated_indices = interpolated_indices[r]
            print len(interpolated_indices)
            
            trajec.position_x[interpolated_indices] = np.nan
            trajec.position_y[interpolated_indices] = np.nan
            
            x = trajec.position_x[0:l:downsample]
            y = trajec.position_y[0:l:downsample]
            t = trajec.time_epoch[0:l:downsample]-trajec.time_epoch[0:l:downsample][0]
            
            if 0:
                fpl.colorline(ax, x, y, t, linewidth=0.5, norm=None, zorder=10, alpha=1, linestyle='solid', cmap=viridis, axis_size_inches=(1.5,1.5), hack_round_caps=True)
                circle_inner = matplotlib.patches.Circle((384, 240), 75, facecolor='red', edgecolor='none', linestyle='solid', linewidth=0.5, zorder=20, alpha=0.25)
                circle_outer = matplotlib.patches.Circle((384, 240), 220, facecolor='gray', edgecolor='black', linestyle='solid', linewidth=0.5, zorder=20, alpha=0.25) 
                ax.add_artist(circle_outer)
                ax.add_artist(circle_inner)
                ax.set_xlim(144,624)
                ax.set_ylim(0,480)
            elif 0:
                fpl.colorline(ax, 480-y, x, t, linewidth=0.5, norm=None, zorder=10, alpha=1, linestyle='solid', cmap=viridis, axis_size_inches=(1.5,1.5), hack_round_caps=True)
                circle_inner = matplotlib.patches.Circle((240, 384), 75, facecolor='red', edgecolor='none', linestyle='solid', linewidth=0.5, zorder=20, alpha=0.25)
                circle_outer = matplotlib.patches.Circle((240, 384), 220, facecolor='gray', edgecolor='black', linestyle='solid', linewidth=0.5, zorder=20, alpha=0.25)
                ax.add_artist(circle_outer)
                ax.add_artist(circle_inner)
                ax.set_ylim(144,624)
                ax.set_xlim(0,480)
            else:
                fpl.colorline(ax, y, x, t, linewidth=0.5, norm=None, zorder=10, alpha=1, linestyle='solid', cmap=viridis, axis_size_inches=(1.5,1.5), hack_round_caps=True)
                circle_inner = matplotlib.patches.Circle((240, 384), 75, facecolor='red', edgecolor='none', linestyle='solid', linewidth=0.5, zorder=20, alpha=0.25)
                circle_outer = matplotlib.patches.Circle((240, 384), 220, facecolor='gray', edgecolor='black', linestyle='solid', linewidth=0.5, zorder=20, alpha=0.25)
                ax.add_artist(circle_outer)
                ax.add_artist(circle_inner)
                ax.set_ylim(144,624)
                ax.set_xlim(0,480)

            

            

    fpl.adjust_spines(ax, [])
    ax.set_frame_on(False)
    
    ax.set_rasterization_zorder(1000)
    #ax.set_ylim(25,455)
    
    
def plot_h2o_trajectory(layout):
    ax = layout.axes_groups['trajectories']['h2o']
    paths = [os.path.join(co2_paper_locations.data_locations.windtunnel_walking['h2o_60sccm'], 'day2')]
    keys = [[299]]
    plot_trajectories(ax, paths, keys, key_lengths=None)
    # dead flies in 200 sccm co2 day 4: 5932, 4394; use lengths: {5932: 9000, 4394: 3000}
    # ethanol trajectory: day2, trajectory 13437 '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_ethanol/day2'
    # h20 trajectory: day2, trajec 299  '/media/Orchard2/singlepad_windtunnel_blacktop/h20_60sccm/day2'
    
def plot_eth_trajectory(layout):
    ax = layout.axes_groups['trajectories']['ethanol']
    paths = [os.path.join(co2_paper_locations.data_locations.windtunnel_walking['eth_60sccm'], 'day3')]
    keys = [[12908]]
    plot_trajectories(ax, paths, keys, key_lengths=None)
    # dead flies in 200 sccm co2 day 4: 5932, 4394; use lengths: {5932: 9000, 4394: 3000}
    # ethanol trajectory: day2, trajectory 13437 '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_ethanol/day2'
    # h20 trajectory: day2, trajec 299  '/media/Orchard2/singlepad_windtunnel_blacktop/h20_60sccm/day2'
    
def plot_co2_trajectory(layout):
    ax = layout.axes_groups['trajectories']['co2']
    paths = [os.path.join(co2_paper_locations.data_locations.windtunnel_walking['co2_60sccm'], 'day5')]
    keys = [[3042]]
    plot_trajectories(ax, paths, keys, key_lengths=None)
    # dead flies in 200 sccm co2 day 4: 5932, 4394; use lengths: {5932: 9000, 4394: 3000}
    # ethanol trajectory: day2, trajectory 13437 '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_ethanol/day2'
    # h20 trajectory: day2, trajec 299  '/media/Orchard2/singlepad_windtunnel_blacktop/h20_60sccm/day2'
    
def plot_vinegar_trajectory(layout):
    ax = layout.axes_groups['trajectories']['vinegar']
    paths = [os.path.join(co2_paper_locations.data_locations.windtunnel_walking['vinegar_60sccm'], 'day1')]
    keys = [[2674]]
    plot_trajectories(ax, paths, keys, key_lengths=None)

def plot_ethacet_trajectory(layout):
    ax = layout.axes_groups['trajectories']['ethyl-acetate']
    paths = [os.path.join(co2_paper_locations.data_locations.windtunnel_walking['ethylacetate_15sccm'], 'day1')]
    keys = [[1448]]
    plot_trajectories(ax, paths, keys, key_lengths=None)

def plot_all_trajectories_on_figure(version='old'):
    layout = get_paper_layout(version=version)
    plot_h2o_trajectory(layout)
    plot_co2_trajectory(layout)
    plot_eth_trajectory(layout)
    #plot_vinegar_trajectory(layout)
    #plot_ethacet_trajectory(layout)
    layout.append_figure_to_layer(layout.figures['trajectories'], 'trajectories', cleartarget=True)
    if version == 'old':
        layout.write_svg(co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_trajecs )
    else:
        layout.write_svg(co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_simple )
    
    
if __name__ == '__main__':
    
    
    plot_all_trajectories_on_figure()
    
    
    
    
    
