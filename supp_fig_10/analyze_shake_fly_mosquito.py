from optparse import OptionParser
import sys, os
import h5py

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np

import pickle
import csv 

import data_fit

import matplotlib.pyplot as plt

import figurefirst
import fly_plot_lib.plot as fpl
import fly_plot_lib.text as flytext

from multi_tracker_analysis import bag2hdf5
import multi_tracker_analysis
get_filenames = multi_tracker_analysis.read_hdf5_file_to_pandas.get_filenames

import co2_paper_locations


    
def get_calibration_factor():
    data = h5py.File(os.path.join(co2_paper_locations.data_locations.shake_fly_mosquito_data,'shake_calibration_co2_2016_09_01_15_34_22_400ppmGAS_327ppmREADING.hdf5'))
    co2 = np.mean(data['phidgets_daq/co2'].value['value'])
    factor = 327/co2
    return factor
    
def get_aligned_data(filename, calibration_factor):
    data = h5py.File(filename)
    print filename
    N_mosquito = float(filename.split('_')[-1].split('.hd')[0].split('N')[-1])
    maxtarg = np.argmax(data['phidgets_daq/co2'].value['value'])
    t = data['phidgets_daq/co2'].value['t'] - data['phidgets_daq/co2'].value['t'][maxtarg]
    co2 = data['phidgets_daq/co2'].value['value']*calibration_factor / N_mosquito
    
    
    t_sync = np.linspace(-60, 3*60, 300)
    co2_sync = np.interp(t_sync, t, co2)
    
    
    return t_sync, co2_sync*100/1000. # 100 mL / min; units in micro liter

def plot_all_mosquitoes():
    fig = plt.figure(figsize=(4,4))

    layout = figurefirst.FigureLayout('supp_fig_10_shakefly.svg')
    layout.make_mplfigures()

    ax = layout.axes[('shake', 'shake')]
    
    calibration_factor = get_calibration_factor()
    
    mosquito_files = get_filenames( co2_paper_locations.data_locations.shake_fly_mosquito_data, 'shake_mosquito', does_not_contain=['~', '.pyc', '.bag', 'calibration'])
    for i, filename in enumerate(mosquito_files):
        t, co2 = get_aligned_data(filename, calibration_factor)
        ax._plot(['Mosquito CO2 release trial: '+str(i+1),
                  'Time',
                  'CO2 concentration'], t, co2, color='blue')
        
    mosquito_files = get_filenames( co2_paper_locations.data_locations.shake_fly_mosquito_data, 'shake_fly', does_not_contain=['~', '.pyc', '.bag', 'calibration'])
    for filename in mosquito_files:
        t, co2 = get_aligned_data(filename, calibration_factor)
        ax._plot(['Fly (drosophila) CO2 release trial: '+str(i+1),
                  'Time',
                  'CO2 concentration'], t, co2, color='black')
        
        
    ax._fill_between(['Time of shaking', 'Time (secs)'], [-30, 0], 0, 1, facecolor='red', edgecolor='none', alpha=0.3)  
    
    
    ax.record = True

    ax.adjust_spines(['left', 'bottom'], spine_locations={'left': 5, 'bottom': 5}, xticks=[-60, -30, 0, 30, 60, 90, 120, 150], yticks=[0, 0.2, 0.4, 0.6, 0.8, 1], linewidth=0.5, tick_length=2.5)
    
    ax.set_ylabel('micro L CO2 / min / animal')
    ax.set_xlabel('Time, sec')
    
    ax.set_fontsize(6)
    
    layout.append_figure_to_layer(layout.figures['shake'], 'shake', cleartarget=True)
    layout.write_svg('supp_fig_10_shakefly.svg')
        
## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
        
    ## Read data #############################################################
    parser = OptionParser()
    parser.add_option('--filename', type=str, help="the .bag file")
    parser.add_option('--max_strlen', type=int, default=255,
                        help="maximum length of encoded strings")
    parser.add_option('--out', type=str, default=None,
                        help="name of output file")
    parser.add_option('--topic', type=str, nargs='*',
                        help="topic name to convert. defaults to all. "
                        "multiple may be specified.")
    (options, args) = parser.parse_args()
           
    fname = os.path.splitext(options.filename)[0]
            
    if options.out is not None:
        output_fname = options.out
    else:
        output_fname = fname + '.hdf5'
            
    print 'Output name: ', output_fname
            
    if not os.path.exists(output_fname):
        if not os.path.exists(options.filename):
            print >> sys.stderr, 'No file %s' % options.filename
            sys.exit(1)
        bag2hdf5.bag2hdf5(options.filename,
                 output_fname,
                 max_strlen=options.max_strlen,
                 topics=['/phidgets_daq/co2', '/signal'])            
    
    f = h5py.File(output_fname, 'r')  
            
    t = f['phidgets_daq']['co2']['t']
    force = f['phidgets_daq']['co2']['value']
    signal = f['signal']
    
    print 'loaded time and force'


    calibration_filename = fname + '_metadata.pickle'
    if os.path.exists(calibration_filename):
        cal = open(calibration_filename)
        regions = pickle.load(cal)
        cal.close()
    else:
        # guess time ranges:
        region_length = 0.5
        
        region1_start = t[0]
        region1_end = region1_start + region_length
        
        region2_start = t[int(len(t)/2.)]
        region2_end = region2_start + region_length
        
        region3_start = t[-1]
        region3_end = region3_start - region_length

        #############################################################################
        
        regions = {'1':  {  'values': [region1_start, region1_end],
                            'indices': [0, 0],
                            'mean': 0,
                            'std': 0},
                   '2':  {  'values': [region2_start, region2_end],
                            'indices': [0, 0],
                            'mean': 0,
                            'std': 0}}
                            
    #############################################################################
    
    try:
        motor_start_stop_filename = fname + '_motor_start_stop.csv'
        motor_start_stop = read_motor_start_stop_data(motor_start_stop_filename)
    except:
        motor_start_stop = None
    calibrationwidget = CalibrationWidget(regions, t, force, calibration_filename, motor_start_stop)
    calibrationwidget.run()
