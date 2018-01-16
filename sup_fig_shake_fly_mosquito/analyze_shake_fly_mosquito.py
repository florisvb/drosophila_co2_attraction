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


'''

class CalibrationWidget(object):
    def __init__(self, regions, t, force, calibration_filename, motor_start_stop):
        self.motor_start_stop = motor_start_stop
        self.regions = regions
        self.t = t
        self.force = force
        self.calibration_filename = calibration_filename
        
        ## create GUI
        self.app = QtGui.QApplication([])
        self.w = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.w.setLayout(self.layout)
        
        btn = QtGui.QPushButton('save calibration')
        btn.pressed.connect(self.save_calibration)
        self.layout.addWidget(btn, 0, 0)
        
        self.p1 = pg.PlotWidget(title="Basic array plotting", x=t, y=force)
        self.p1.enableAutoRange('xy', False)
        self.layout.addWidget(self.p1, 0, 1)
        
        for r, region in self.regions.items():
            lr = pg.LinearRegionItem(values=region['values'])
            f = 'update_linear_region_' + str(r)
            lr.sigRegionChanged.connect(self.__getattribute__(f))
            self.p1.addItem(lr)
        
        if motor_start_stop is not None:
            for l in self.motor_start_stop:
                print 'making infinite line: ', l, self.t[0], self.t[-1]
                il = pg.InfiniteLine(pos=l, angle=90)
                self.p1.addItem(il)
            
    def update_linear_region(self, r, linear_region):
        self.regions[r]['values'] = linear_region.getRegion()
        indices = [0, 0]
        indices[0] = np.argmin( np.abs( self.t - self.regions[r]['values'][0] ) )
        indices[1] = np.argmin( np.abs( self.t - self.regions[r]['values'][1] ) )
        self.regions[r]['indices'] = indices
        self.regions[r]['mean'] = np.mean(self.force[self.regions[r]['indices'][0]: self.regions[r]['indices'][-1]])
        self.regions[r]['std'] = np.std(self.force[self.regions[r]['indices'][0]: self.regions[r]['indices'][-1]])
        print 'region: ', r, ' :', self.regions[r]['mean']
        
    def update_linear_region_1(self, linear_region):
        self.update_linear_region('1', linear_region)

    def update_linear_region_2(self, linear_region):
        self.update_linear_region('2', linear_region)    

    def run(self):

        ## Display the widget as a new window
        self.w.show()

        ## Start the Qt event loop
        self.app.exec_()
        
    def save_calibration(self):
        f = open(self.calibration_filename, 'w')
        pickle.dump(self.regions, f)
        f.close()
        print 'calibration saved'

def read_motor_start_stop_data(filename):
    values = []
    with open(filename, 'rb') as csvfile:
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in datareader:
            values.append(float(row[1]))
    return values
    

def load_bottle_data_and_align_peaks(filename):
    fname = os.path.splitext(filename)[0]
    output_fname = fname + '.hdf5'
            
    print 'Output name: ', output_fname
            
    if not os.path.exists(output_fname):
        if not os.path.exists(filename):
            print >> sys.stderr, 'No file %s' % filename
            sys.exit(1)
        bag2hdf5.bag2hdf5(filename,
                 output_fname,
                 max_strlen=100,
                 topics=['/phidgets_daq/co2'])            
    
    f = h5py.File(output_fname, 'r')  
            
    t = f['phidgets_daq']['co2']['t']
    force = f['phidgets_daq']['co2']['value']
    
    # calibrate
    f = open(co2_paper_locations.data_locations.shake_fly_mosquito_data+'/400ppm_co2_calibration_2016_07_11_16_45_48_metadata.pickle')
    calibration_400ppm = pickle.load(f)['1']['mean']
    calibration_zeroppm = 0
    
    if 'bottle' in os.path.basename(filename):
        force *= 2
    
    force -= calibration_zeroppm
    force /= calibration_400ppm
    force *= 400
    
    
    return t, force

def plot_calibration_values(ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    co2_1percent_files = [co2_paper_locations.data_locations.shake_fly_mosquito_data+'/1percent_co2_2016_07_11_15_38_55.bag' , co2_paper_locations.data_locations.shake_fly_mosquito_data+'/1percent_co2_2016_07_11_15_51_11.bag' , co2_paper_locations.data_locations.shake_fly_mosquito_data+'/1percent_co2_2016_07_11_17_40_43.bag']
    co2_2000ppm_files = [co2_paper_locations.data_locations.shake_fly_mosquito_data+'/2000ppm_2_co2_2016_07_11_17_15_44.bag' , co2_paper_locations.data_locations.shake_fly_mosquito_data+'/2000ppm_co2_2016_07_11_16_07_01.bag' , co2_paper_locations.data_locations.shake_fly_mosquito_data+'/2000ppm_co2_2016_07_11_16_15_26.bag' , co2_paper_locations.data_locations.shake_fly_mosquito_data+'/2000ppm_co2_2016_07_11_16_25_22.bag' ]
    co2_400ppm_files = [co2_paper_locations.data_locations.shake_fly_mosquito_data+'/400ppm_co2_2016_07_11_16_50_23.bag' , co2_paper_locations.data_locations.shake_fly_mosquito_data+'/400ppm_co2_2016_07_11_16_56_37.bag' , co2_paper_locations.data_locations.shake_fly_mosquito_data+'/400ppm_co2_2016_07_11_17_03_21.bag']
    
    co2_1percent_peaks = []
    co2_2000ppm_peaks = []
    co2_400ppm_peaks = []
    
    for filename in co2_1percent_files:
        t, f = load_bottle_data_and_align_peaks(filename)
        co2_1percent_peaks.append(np.max(f))
    
    for filename in co2_2000ppm_files:
        t, f = load_bottle_data_and_align_peaks(filename)
        co2_2000ppm_peaks.append(np.max(f))
    
    for filename in co2_400ppm_files:
        t, f = load_bottle_data_and_align_peaks(filename)
        co2_400ppm_peaks.append(np.max(f))
        
    
    co2_1percent_actual = [10000 for i in range(len(co2_1percent_peaks))]
    co2_2000ppm_actual = [2000 for i in range(len(co2_2000ppm_peaks))]
    co2_400ppm_actual = [400 for i in range(len(co2_400ppm_peaks))]
    
    ax.plot(co2_400ppm_peaks, co2_400ppm_actual, 'o', color='black')
    ax.plot(co2_2000ppm_peaks, co2_2000ppm_actual, 'o', color='black')
    ax.plot(co2_1percent_peaks, co2_1percent_actual, 'o', color='black')
    
    all_measured_values = []
    all_measured_values.extend(co2_400ppm_peaks)
    all_measured_values.extend(co2_2000ppm_peaks)
    all_measured_values.extend(co2_1percent_peaks)
    all_measured_values = np.array(all_measured_values)
    
    all_actual_values = []
    all_actual_values.extend(co2_400ppm_actual)
    all_actual_values.extend(co2_2000ppm_actual)
    all_actual_values.extend(co2_1percent_actual)
    all_actual_values = np.array(all_actual_values)
    
    model = data_fit.models.LinearModel()
    model.fit(all_actual_values, inputs=all_measured_values)
    
    x = np.arange(0, 3000)
    y = model.get_val(x)
    ax.plot(x,y, color='black', linewidth=2)
    
    
    # bottles
    q200mg_bottle_files = ['/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_2dayold_25dincubator/bottle4_co2_2016_07_16_13_38_08.bag',
                     '/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_2dayold_25dincubator/bottle5_co2_2016_07_16_13_41_31.bag',
                    '/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_2dayold_25dincubator/bottle6_co2_2016_07_16_13_45_47.bag',
                    '/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_2dayold_my_unvented_incubator/bottle_novent_200mg_1_co2_2016_07_20_10_29_39.bag' , 
                    '/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_2dayold_my_unvented_incubator/bottle_novent_200mg_2_co2_2016_07_20_10_32_25.bag' , 
                    '/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_2dayold_my_unvented_incubator/bottle_novent_200mg_3_co2_2016_07_20_10_35_12.bag',
                    ]
    
    q400mg_bottle_files = [ '/media/Orchard3/fly_bottle_data/bottles/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_1_co2_2016_07_18_14_09_38.bag', 
                           '/media/Orchard3/fly_bottle_data/bottles/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_2_co2_2016_07_18_14_11_45.bag', 
                           '/media/Orchard3/fly_bottle_data/bottles/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_3_co2_2016_07_18_14_14_43.bag',
                    ]
    
    q200mg_bottle_flies_files = [  '/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_1_co2_2016_08_24_16_20_16.bag',
                                   '/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_2_co2_2016_08_24_16_22_32.bag',
                                   '/media/Orchard3/fly_bottle_data/bottles/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_3_co2_2016_08_24_16_25_19.bag', 
                    ]
                    
                    
    bottle_values_200mg = []
    for bottle in q200mg_bottle_files:
        t, f = load_bottle_data_and_align_peaks(bottle)
        ax.vlines( np.max(f), 0, model.get_val(np.max(f)), color='blue' )
        ax.hlines( model.get_val(np.max(f)), 0,  np.max(f), color='blue')
        bottle_values_200mg.append(model.get_val(np.max(f)))
        
    bottle_values_400mg = []
    for bottle in q400mg_bottle_files:
        t, f = load_bottle_data_and_align_peaks(bottle)
        ax.vlines( np.max(f), 0, model.get_val(np.max(f)), color='green' )
        ax.hlines( model.get_val(np.max(f)), 0,  np.max(f), color='green')
        bottle_values_400mg.append(model.get_val(np.max(f)))
        
    bottle_values_200mg_flies = []
    for bottle in q200mg_bottle_flies_files:
        t, f = load_bottle_data_and_align_peaks(bottle)
        ax.vlines( np.max(f), 0, model.get_val(np.max(f)), color='red' )
        ax.hlines( model.get_val(np.max(f)), 0,  np.max(f), color='red')
        bottle_values_200mg_flies.append(model.get_val(np.max(f)))
    
    yticks = [400, 5000, 10000]
    xticks = [0, 1000]
    figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], yticks=yticks, xticks=xticks)

    return np.array(bottle_values_200mg), np.array(bottle_values_400mg), np.array(bottle_values_200mg_flies)

def get_paper_layout():
    svg = '/media/Orchard2/fermentation/figure_output.svg'
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout

def make_plot():
    layout = get_paper_layout()
    
    bottle_values_200mg, bottle_values_400mg, bottle_values_200mg_flies = plot_calibration_values() # add ax if you want to direct the plot
    
    fpl.scatter_box(layout.axes_groups['bottles']['fly_bottles'], 1, bottle_values_200mg, xwidth=0.4, color='blue', use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)
    fpl.scatter_box(layout.axes_groups['bottles']['fly_bottles'], 2, bottle_values_400mg, xwidth=0.4, color='green', use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)
    fpl.scatter_box(layout.axes_groups['bottles']['fly_bottles'], 3, bottle_values_200mg_flies, xwidth=0.4, color='red', use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)
    
    yticks = [400, 5000, 10000]
    layout.axes_groups['bottles']['fly_bottles'].set_ylim(0,15000)
    #figurefirst.mpl_functions.adjust_spines(layout.axes_groups['bottles']['fly_bottles'], ['left'], yticks=yticks)
    figurefirst.mpl_functions.adjust_spines(layout.axes_groups['bottles']['fly_bottles'], [])
    
    #ytick_labels = [str(ytick/float(1e6)*100) + '%' for ytick in yticks]
    #layout.axes_groups['bottles']['fly_bottles'].set_yticklabels(ytick_labels)
    layout.axes_groups['bottles']['fly_bottles'].set_yticklabels([])
    
    layout.append_figure_to_layer(layout.figures['bottles'], 'bottles', cleartarget=True)
    layout.write_svg( '/media/Orchard2/fermentation/figure_output.svg')
'''
     
    
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
    ax = fig.add_axes([0.15,0.15, 0.8, 0.8])
    
    calibration_factor = get_calibration_factor()
    
    mosquito_files = get_filenames( co2_paper_locations.data_locations.shake_fly_mosquito_data, 'shake_mosquito', does_not_contain=['~', '.pyc', '.bag', 'calibration'])
    for filename in mosquito_files:
        t, co2 = get_aligned_data(filename, calibration_factor)
        ax.plot(t, co2, color='blue')
        
    mosquito_files = get_filenames( co2_paper_locations.data_locations.shake_fly_mosquito_data, 'shake_fly', does_not_contain=['~', '.pyc', '.bag', 'calibration'])
    for filename in mosquito_files:
        t, co2 = get_aligned_data(filename, calibration_factor)
        ax.plot(t, co2, color='black')
        
        
    ax.fill_between([-30, 0], 0, 1, facecolor='red', edgecolor='none', alpha=0.3)  
    
    
    
    fpl.adjust_spines(ax, ['left', 'bottom'], spine_locations={'left': 5, 'bottom': 5}, xticks=[-60, -30, 0, 30, 60, 90, 120, 150], yticks=[0, 0.2, 0.4, 0.6, 0.8, 1])
    
    ax.set_ylabel('micro L CO2 / min / animal')
    ax.set_xlabel('Time, sec')
    
    flytext.set_fontsize(fig, 8)
    
    fig.savefig('CO2_mosquito_fly.pdf', format='pdf')
        
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
