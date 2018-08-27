from optparse import OptionParser
import sys, os
import bag2hdf5, h5py

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

import figurefirst

import co2_paper_locations
data_locations = co2_paper_locations.data_locations
figure_template_locations = co2_paper_locations.figure_template_locations

fifidatafile = os.path.join(os.path.dirname(figure_template_locations.sup_figure_co2_calibration), 'figure_data.pickle')

def get_supplemental_layout():
    svg = figure_template_locations.sup_figure_co2_calibration
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout

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
    

def plot_calibration_supplemental_figure():
    layout = get_supplemental_layout()
    plot_calibration_concentrations(layout, 'co2_calibration', 'calibration_timecourse')
    plot_calibration_values(layout, 'co2_calibration', 'calibration_curve')
    #ax = layout.axes[('co2_calibration', 'calibration_curve')].figure
    #flytext.set_fontsize(ax.figure, 6)
    layout.append_figure_to_layer(layout.figures['co2_calibration'], 'co2_calibration', cleartarget=True)
    layout.write_svg(figure_template_locations.sup_figure_co2_calibration )

def plot_calibration_concentrations(layout, figure, axis):
    # plot_calibration_concentrations(layout, 'co2_calibration', 'calibration_timecourse')

    co2_1percent_files = [data_locations.fly_bottle_calibration+'/1percent_co2_2016_07_11_15_38_55.bag' , data_locations.fly_bottle_calibration+'/1percent_co2_2016_07_11_15_51_11.bag' , data_locations.fly_bottle_calibration+'/1percent_co2_2016_07_11_17_40_43.bag']
    co2_2000ppm_files = [data_locations.fly_bottle_calibration+'/2000ppm_2_co2_2016_07_11_17_15_44.bag' , data_locations.fly_bottle_calibration+'/2000ppm_co2_2016_07_11_16_07_01.bag' , data_locations.fly_bottle_calibration+'/2000ppm_co2_2016_07_11_16_15_26.bag' , data_locations.fly_bottle_calibration+'/2000ppm_co2_2016_07_11_16_25_22.bag' ]
    co2_400ppm_files = [data_locations.fly_bottle_calibration+'/400ppm_co2_2016_07_11_16_50_23.bag' , data_locations.fly_bottle_calibration+'/400ppm_co2_2016_07_11_16_56_37.bag' , data_locations.fly_bottle_calibration+'/400ppm_co2_2016_07_11_17_03_21.bag']
    
    for f, filename in enumerate(co2_1percent_files):
        t, f = load_bottle_data_and_align_peaks(filename)
        t0 = t[np.argmax(f)]
        #ax.plot(t-t0, f, color='black')
        figurefirst.deprecated_regenerate.mpl('plot', layout, figure, axis, fifidatafile, 'co2 1 percent calibration sample: '+str(f),
                              ['Time',
                               'Measured PPM for 1 percent calibration gas'],
                               t-t0, f, color='black')
    
    for f, filename in enumerate(co2_2000ppm_files):
        t, f = load_bottle_data_and_align_peaks(filename)
        t0 = t[np.argmax(f)]
        #ax.plot(t-t0, f, color='black')
        figurefirst.deprecated_regenerate.mpl('plot', layout, figure, axis, fifidatafile, 'co2 2000 ppm sample: '+str(f),
                              ['Time',
                               'Measured PPM for 2000 PPM calibration gas'],
                               t-t0, f, color='black')

    for f, filename in enumerate(co2_400ppm_files):
        t, f = load_bottle_data_and_align_peaks(filename)
        t0 = t[np.argmax(f)]
        #ax.plot(t-t0, f, color='black')
        figurefirst.deprecated_regenerate.mpl('plot', layout, figure, axis, fifidatafile, 'co2 400 ppm sample: '+str(f),
                              ['Time',
                               'Measured PPM for 400 PPM calibration gas'],
                               t-t0, f, color='black')

    #figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], xticks=[0, 300, 600], yticks=[0,1000,2000], linewidth=0.5, spine_locations={'left': 5, 'bottom': 5})
    #ax.tick_params(length=2.5)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust spines timecourse', [], 
                                  ['left', 'bottom'], xticks=[0, 300, 600], yticks=[0,1000,2000], linewidth=0.5, spine_locations={'left': 5, 'bottom': 5}, tick_length=2.5)

    #flytext.set_fontsize(ax.figure, 6)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, figure, axis, fifidatafile, 'set fontsize', [], 
                                  6)

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
    f = open(data_locations.fly_bottle_calibration + '/400ppm_co2_calibration_2016_07_11_16_45_48_metadata.pickle')
    calibration_400ppm = pickle.load(f)['1']['mean']
    calibration_zeroppm = 0
    
    if 'bottle' in os.path.basename(filename):
        force *= 2
    
    force -= calibration_zeroppm
    force /= calibration_400ppm
    force *= 400
    
    
    return t, force

def plot_calibration_values(layout, figure, axis):
        
    co2_1percent_files = [data_locations.fly_bottle_calibration + '/1percent_co2_2016_07_11_15_38_55.bag' , data_locations.fly_bottle_calibration + '/1percent_co2_2016_07_11_15_51_11.bag' , data_locations.fly_bottle_calibration + '/1percent_co2_2016_07_11_17_40_43.bag']
    co2_2000ppm_files = [data_locations.fly_bottle_calibration + '/2000ppm_2_co2_2016_07_11_17_15_44.bag' , data_locations.fly_bottle_calibration + '/2000ppm_co2_2016_07_11_16_07_01.bag' , data_locations.fly_bottle_calibration + '/2000ppm_co2_2016_07_11_16_15_26.bag' , data_locations.fly_bottle_calibration + '/2000ppm_co2_2016_07_11_16_25_22.bag' ]
    co2_400ppm_files = [data_locations.fly_bottle_calibration + '/400ppm_co2_2016_07_11_16_50_23.bag' , data_locations.fly_bottle_calibration + '/400ppm_co2_2016_07_11_16_56_37.bag' , data_locations.fly_bottle_calibration + '/400ppm_co2_2016_07_11_17_03_21.bag']
    
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
    
    #ax.plot(co2_400ppm_peaks, co2_400ppm_actual, 'o', color='black', markersize=2)
    figurefirst.deprecated_regenerate.mpl('plot', layout, figure, axis, fifidatafile, 'co2 400 ppm',
                              ['Measured PPM for 400 PPM calibration gas',
                               'Actual PPM'],
                               co2_400ppm_peaks, co2_400ppm_actual, 'o', color='black', markersize=2)

    #ax.plot(co2_2000ppm_peaks, co2_2000ppm_actual, 'o', color='black', markersize=2)
    figurefirst.deprecated_regenerate.mpl('plot', layout, figure, axis, fifidatafile, 'co2 2000 ppm',
                              ['Measured PPM for 2000 PPM calibration gas',
                               'Actual PPM'],
                               co2_2000ppm_peaks, co2_2000ppm_actual, 'o', color='black', markersize=2)

    #ax.plot(co2_1percent_peaks, co2_1percent_actual, 'o', color='black', markersize=2)
    figurefirst.deprecated_regenerate.mpl('plot', layout, figure, axis, fifidatafile, 'co2 1 percent',
                              ['Measured PPM for 1 percent calibration gas',
                               'Actual PPM'],
                               co2_1percent_peaks, co2_1percent_actual, 'o', color='black', markersize=2)
    
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
    #ax.plot(x,y, color='black', linewidth=1, linestyle=':')
    figurefirst.deprecated_regenerate.mpl('plot', layout, figure, axis, fifidatafile, 'linear model',
                              ['linear model (measured value)',
                               'linear model (actual value)'],
                               x,y, color='black', linewidth=1, linestyle=':')
    
    
    # bottles
    q200mg_bottle_files = [data_locations.fly_bottle_data + '/200mg_bakersyeast_2dayold_25dincubator/bottle4_co2_2016_07_16_13_38_08.bag',
                     data_locations.fly_bottle_data + '/200mg_bakersyeast_2dayold_25dincubator/bottle5_co2_2016_07_16_13_41_31.bag',
                    data_locations.fly_bottle_data + '/200mg_bakersyeast_2dayold_25dincubator/bottle6_co2_2016_07_16_13_45_47.bag',
                    data_locations.fly_bottle_data + '/200mg_bakersyeast_2dayold_my_unvented_incubator/bottle_novent_200mg_1_co2_2016_07_20_10_29_39.bag' , 
                    data_locations.fly_bottle_data + '/200mg_bakersyeast_2dayold_my_unvented_incubator/bottle_novent_200mg_2_co2_2016_07_20_10_32_25.bag' , 
                    data_locations.fly_bottle_data + '/200mg_bakersyeast_2dayold_my_unvented_incubator/bottle_novent_200mg_3_co2_2016_07_20_10_35_12.bag',
                    ]
    
    q400mg_bottle_files = [ data_locations.fly_bottle_data + '/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_1_co2_2016_07_18_14_09_38.bag', 
                           data_locations.fly_bottle_data + '/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_2_co2_2016_07_18_14_11_45.bag', 
                           data_locations.fly_bottle_data + '/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_3_co2_2016_07_18_14_14_43.bag',
                           data_locations.fly_bottle_data + '/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_4_co2_2016_08_10_16_44_32.bag', 
                           data_locations.fly_bottle_data + '/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_5_co2_2016_08_10_16_47_30.bag', 
                           data_locations.fly_bottle_data + '/400mg_bakersyeast_2dayold_25degincubator/bottle_400mg_6_co2_2016_08_10_16_50_48.bag', 
                    ]
    
    q200mg_bottle_flies_files = [  data_locations.fly_bottle_data + '/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_1_co2_2016_08_24_16_20_16.bag',
                                   data_locations.fly_bottle_data + '/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_2_co2_2016_08_24_16_22_32.bag',
                                   data_locations.fly_bottle_data + '/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_3_co2_2016_08_24_16_25_19.bag', 
                                   data_locations.fly_bottle_data + '/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_4_co2_2017_02_01_09_26_22.bag',
                                   data_locations.fly_bottle_data + '/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_5_co2_2017_02_01_09_28_28.bag',
                                   data_locations.fly_bottle_data + '/200mg_bakersyeast_12days_adult_flies_normalstudding/bottle_flies_200mg_6_co2_2017_02_01_09_30_46.bag',
                    ]
    
    qnoyeast_bottle_files = [  data_locations.fly_bottle_data + '/noyeast/noyeast_1_2016_11_05_15_37_01.bag' , 
                               data_locations.fly_bottle_data + '/noyeast/noyeast_2_2016_11_05_15_38_46.bag' , 
                                data_locations.fly_bottle_data + '/noyeast/noyeast_3_2016_11_05_15_40_16.bag' , 
                                data_locations.fly_bottle_data + '/noyeast/noyeast_4_2016_11_07_15_11_51.bag' , 
                                data_locations.fly_bottle_data + '/noyeast/noyeast_5_2016_11_07_15_13_32.bag' , 
                                data_locations.fly_bottle_data + '/noyeast/noyeast_6_2016_11_07_15_15_20.bag' ,
                    ]
    
    
    bottle_values_noyeast = []
    for b, bottle in enumerate(qnoyeast_bottle_files):
        t, f = load_bottle_data_and_align_peaks(bottle)
        #ax.vlines( np.max(f), 0, model.get_val(np.max(f)), color='brown', linewidth=0.5 )
        figurefirst.deprecated_regenerate.mpl('vlines', layout, figure, axis, fifidatafile, 'no yeast vlines '+str(b), [],
                                  np.max(f), 0, model.get_val(np.max(f)), color='brown', linewidth=0.5 )
        #ax.hlines( model.get_val(np.max(f)), 0,  np.max(f), color='brown', linewidth=0.5)
        figurefirst.deprecated_regenerate.mpl('hlines', layout, figure, axis, fifidatafile, 'no yeast hlines '+str(b), [],
                                  model.get_val(np.max(f)), 0,  np.max(f), color='brown', linewidth=0.5 )
        bottle_values_noyeast.append(model.get_val(np.max(f)))                
                    
    bottle_values_200mg = []
    for b, bottle in enumerate(q200mg_bottle_files):
        t, f = load_bottle_data_and_align_peaks(bottle)
        #ax.vlines( np.max(f), 0, model.get_val(np.max(f)), color='blue', linewidth=0.5 )
        figurefirst.deprecated_regenerate.mpl('vlines', layout, figure, axis, fifidatafile, '200 mg vlines '+str(b), [],
                                  np.max(f), 0, model.get_val(np.max(f)), color='blue', linewidth=0.5 )
        #ax.hlines( model.get_val(np.max(f)), 0,  np.max(f), color='blue', linewidth=0.5)
        figurefirst.deprecated_regenerate.mpl('hlines', layout, figure, axis, fifidatafile, '200 mg hlines '+str(b), [],
                                  model.get_val(np.max(f)), 0,  np.max(f), color='blue', linewidth=0.5 )
        bottle_values_200mg.append(model.get_val(np.max(f)))
        
    bottle_values_400mg = []
    for b, bottle in enumerate(q400mg_bottle_files):
        t, f = load_bottle_data_and_align_peaks(bottle)
        #ax.vlines( np.max(f), 0, model.get_val(np.max(f)), color='green', linewidth=0.5 )
        figurefirst.deprecated_regenerate.mpl('vlines', layout, figure, axis, fifidatafile, '400 mg vlines '+str(b), [],
                                  np.max(f), 0, model.get_val(np.max(f)), color='green', linewidth=0.5 )
        #ax.hlines( model.get_val(np.max(f)), 0,  np.max(f), color='green', linewidth=0.5)
        figurefirst.deprecated_regenerate.mpl('hlines', layout, figure, axis, fifidatafile, '400 mg hlines '+str(b), [],
                                  model.get_val(np.max(f)), 0,  np.max(f), color='green', linewidth=0.5 )
        bottle_values_400mg.append(model.get_val(np.max(f)))
        
    bottle_values_200mg_flies = []
    for b, bottle in enumerate(q200mg_bottle_flies_files):
        t, f = load_bottle_data_and_align_peaks(bottle)
        #ax.vlines( np.max(f), 0, model.get_val(np.max(f)), color='red', linewidth=0.5 )
        figurefirst.deprecated_regenerate.mpl('vlines', layout, figure, axis, fifidatafile, '200 mg flies vlines '+str(b), [],
                                  np.max(f), 0, model.get_val(np.max(f)), color='red', linewidth=0.5 )
        #ax.hlines( model.get_val(np.max(f)), 0,  np.max(f), color='red', linewidth=0.5)
        figurefirst.deprecated_regenerate.mpl('hlines', layout, figure, axis, fifidatafile, '200 mg flies hlines '+str(b), [],
                                  model.get_val(np.max(f)), 0,  np.max(f), color='red', linewidth=0.5 )
        bottle_values_200mg_flies.append(model.get_val(np.max(f)))
    
    yticks = [400, 2000, 10000]
    xticks = [0, 1000, 2000]
    #figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], yticks=yticks, xticks=xticks, linewidth=0.5, spine_locations={'left': 5, 'bottom': 5})
    #ax.tick_params(length=2.5)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust spines', [], 
                                  ['left', 'bottom'], yticks=yticks, xticks=xticks, linewidth=0.5, spine_locations={'left': 5, 'bottom': 5}, tick_length=2.5)


    #ax.set_xlabel('Measured CO2 concentration, ppm')
    #ax.set_ylabel('CO2 concentration in bottle, ppm')
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, figure, axis, fifidatafile, 'set fontsize', [], 
                                  6)

    #fig.savefig('Calibration_curve.pdf', format='pdf')

    return np.array(bottle_values_noyeast), np.array(bottle_values_200mg), np.array(bottle_values_400mg), np.array(bottle_values_200mg_flies)

def get_paper_layout():
    svg = figure_template_locations.figure1_ferment_trapassay_flybottles
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout

def make_plot():
    layout = get_paper_layout()
    figure = 'bottles'
    axis = 'fly_bottles'
    
    bottle_values_noyeast, bottle_values_200mg, bottle_values_400mg, bottle_values_200mg_flies = plot_calibration_values(layout, 'co2_calibration', 'calibration_curve') # add ax if you want to direct the plot
    black = (0.0001, 0.0001, 0.0001)
    #fpl.scatter_box(layout.axes_groups['bottles']['fly_bottles'], 0, bottle_values_noyeast, xwidth=0.4, color=black, use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_box',
                                 layout, figure, axis, fifidatafile, 'no yeast',
                                {'x': 'index',
                                 'CO2': 'CO2 in bottle with no yeast'},
                                 0, bottle_values_noyeast, 
                                 xwidth=0.4, color=black, use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)

    #fpl.scatter_box(layout.axes_groups['bottles']['fly_bottles'], 1, bottle_values_200mg, xwidth=0.4, color=black, use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_box',
                                 layout, figure, axis, fifidatafile, '200 mg',
                                ['index',
                                 'CO2 in bottle with 200 mg'],
                                 1, bottle_values_200mg, 
                                 xwidth=0.4, color=black, use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)

    #fpl.scatter_box(layout.axes_groups['bottles']['fly_bottles'], 2, bottle_values_400mg, xwidth=0.4, color=black, use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_box',
                                 layout, figure, axis, fifidatafile, '400 mg',
                                ['index',
                                 'CO2 in bottle with 400 mg yeast'],
                                 2, bottle_values_400mg, 
                                 xwidth=0.4, color=black, use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)

    #fpl.scatter_box(layout.axes_groups['bottles']['fly_bottles'], 3, bottle_values_200mg_flies, xwidth=0.4, color=black, use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_box',
                                 layout, figure, axis, fifidatafile, '200mg flies',
                                ['index',
                                 'CO2 in bottle with 200mg yeast and flies'],
                                 3, bottle_values_200mg_flies, 
                                 xwidth=0.4, color=black, use='mean', optimize_scatter_distance=True, optimize_scatter_distance_y_scale=0.001, markersize=2)

    #layout.axes_groups['bottles']['fly_bottles'].hlines([400], [-1,5], colors=['red'], linestyles=['dashed'])
    x = np.linspace(-1,4)
    y = np.ones_like(x)*400
    #layout.axes_groups['bottles']['fly_bottles'].plot(x,y,linestyle='--',color='red',linewidth=1,zorder=-100, dashes = [2,2,2,2])
    figurefirst.deprecated_regenerate.mpl('plot',
                                 layout, figure, axis, fifidatafile, 'atmospheric line',
                                ['index',
                                 '400 ppm atmospheric level'],
                                 x,y,linestyle='--',color='red',linewidth=1,zorder=-100, dashes = [2,2,2,2])

    yticks = [0, 5000, 10000, 15000]
    xticks = [0,1,2,3]
    #layout.axes_groups['bottles']['fly_bottles'].set_ylim(0,15000)
    figurefirst.deprecated_regenerate.mpl('set_ylim', layout, figure, axis, fifidatafile, 'y lim', [], 
                              0,15000)

    #figurefirst.mpl_functions.adjust_spines(layout.axes_groups['bottles']['fly_bottles'], ['left', 'bottom'], yticks=yticks, xticks=xticks, spine_locations={'left': 5, 'bottom': 5}, linewidth=0.5)
    #layout.axes_groups['bottles']['fly_bottles'].tick_params(length=2.5)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust spines', [], 
                                 ['left', 'bottom'], yticks=yticks, xticks=xticks, spine_locations={'left': 5, 'bottom': 5}, linewidth=0.5, tick_length=2.5)


    
    #layout.axes_groups['bottles']['fly_bottles'].set_yticklabels(['0', '0.5', '1', '1.5'])
    figurefirst.deprecated_regenerate.mpl('set_yticklabels', layout, figure, axis, fifidatafile, 'y labels', [], 
                              ['0', '0.5', '1', '1.5'])

    #layout.axes_groups['bottles']['fly_bottles'].set_xticklabels([])
    figurefirst.deprecated_regenerate.mpl('set_xticklabels', layout, figure, axis, fifidatafile, 'x labels', [], 
                              [])

    #flytext.set_fontsize(layout.axes_groups['bottles']['fly_bottles'].figure, 6)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, figure, axis, fifidatafile, 'font size', [], 6)


    
    layout.append_figure_to_layer(layout.figures['bottles'], 'bottles', cleartarget=True)
    layout.write_svg( figure_template_locations.figure1_ferment_trapassay_flybottles)
    
    
    
    
    
    
## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    make_plot()
    ## Read data #############################################################
    '''
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
                 topics=['/phidgets_daq/co2'])            
    
    f = h5py.File(output_fname, 'r')  
            
    t = f['phidgets_daq']['co2']['t']
    force = f['phidgets_daq']['co2']['value']
    
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
    '''
