from multi_tracker_analysis import bag2hdf5
import os
import h5py
import numpy as np
import multi_tracker_analysis as mta
import data_fit
import figurefirst
import fly_plot_lib.text as flytext

import co2_paper_locations

def load_calibrations():
    calibration_2000ppm = os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_2000ppm_calibration_2017_04_05_17_35_16.bag')
    calibration_0ppm = os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_0ppm_calibration_2017_04_05_17_36_37.bag')

    co2_2000ppm_mean, co2_2000ppm_std = load_raw_co2_mean_and_std(calibration_2000ppm)
    co2_0ppm_mean, co2_0ppm_std = load_raw_co2_mean_and_std(calibration_0ppm)

    return co2_0ppm_mean, co2_2000ppm_mean

def load_raw_co2_mean_and_std(filename):
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
    co2 = f['phidgets_daq']['co2']['value']

    return np.mean(co2), np.std(co2)

def load_calibrated_co2_mean_and_std(filename):

    co2_0ppm_mean, co2_2000ppm_mean = load_calibrations()

    co2_mean, co2_std = load_raw_co2_mean_and_std(filename)

    co2_mean -= co2_0ppm_mean
    co2_mean /= co2_2000ppm_mean
    co2_mean *= 2000 # in ppm

    co2_std /= co2_2000ppm_mean
    co2_std *= 2000 # in ppm

    return co2_mean, co2_std

def get_co2_series_at_distance(distance=2, concentration=5):

    distance = '_d' + str(distance) + '_'
    concentration = '_c' + str(concentration) + '_'
    filenames = mta.read_hdf5_file_to_pandas.get_filenames(co2_paper_locations.data_locations.co2_calibration_windtunnel, distance)

    means = {}
    stds = {}
    for filename in filenames:
        if concentration not in filename:
            continue

        co2_mean, co2_std = load_calibrated_co2_mean_and_std(filename)


        elevation = filename.split(distance)[-1].split('_')[0]
        if 'm' in elevation:
            elevation = float(elevation[1:-1])*-1
        else:
            elevation = float(elevation[1:])

        means[elevation] = co2_mean
        stds[elevation] = co2_std

    return means, stds

def plot_co2_series_at_distance_on_ax(ax, sccm=60, distance=2, concentration=5):
    
    model = plot_concentration_profile(None)

    means, stds = get_co2_series_at_distance(distance, concentration)
    elevations = means.keys()
    elevations.sort()
    
    mean_vals = np.array([means[e] for e in elevations])
    std_vals = np.array([stds[e] for e in elevations])

    intercept = model.parameters['intercept']
    slope_2cm = model.parameters['slope']

    if distance == 10:
        co2_mean, co2_std = load_calibrated_co2_mean_and_std(os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_c5_d10_e0_2017_04_05_14_13_01.bag'))
        slope_10cm = (co2_mean - intercept) / 5.
    elif distance == 0:
        co2_mean, co2_std = load_calibrated_co2_mean_and_std(os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_c1_dcenter_e0_2017_04_05_17_32_12.bag'))
        slope_0cm = (co2_mean - intercept) / 1.

    if distance == 2:
        model.parameters['slope'] = slope_2cm
    elif distance == 10:
        model.parameters['slope'] = slope_10cm
    elif distance == 0:
        mode.parameters['slope'] = slope_0cm

    percent_co2_at_0elevation = (model.parameters['slope']*sccm + model.parameters['intercept'])/float(1e6)*100.

    mean_vals = mean_vals / np.max(mean_vals) * percent_co2_at_0elevation

    #ax.fill_betweenx(elevations, mean_vals-std_vals, mean_vals+std_vals, 
    #    facecolor='red', edgecolor='none', alpha=0.3)

    if distance == 2:
        ax.plot(mean_vals, elevations, color='red')
    elif distance == 10:
        ax.plot(mean_vals, elevations, color='purple')

def plot_concentration_profile(ax, distance=2, elevation=0):

    distance = '_d' + str(distance) + '_'
    elevation = '_e' + str(elevation) + '_'
    filenames = mta.read_hdf5_file_to_pandas.get_filenames(co2_paper_locations.data_locations.co2_calibration_windtunnel, distance)

    means = {}
    stds = {}
    for filename in filenames:
        if elevation not in filename:
            continue    

        co2_mean, co2_std = load_calibrated_co2_mean_and_std(filename)

        concentration = filename.split('windtunnel_')[-1].split('_')[0]
        concentration = float(concentration[1:])

        means[concentration] = co2_mean
        stds[concentration] = co2_std

    concentrations = means.keys()
    concentrations.sort()

    mean_vals = [means[c] for c in concentrations]
    std_vals = [stds[c] for c in concentrations]

    if ax is not None:
        ax.plot(concentrations, mean_vals, 'o', color='red', markersize=3, zorder=100)

    # root model
    if 0:
        model = data_fit.models.RootModel()
        model.fit(np.array(mean_vals), np.array(concentrations))
        x = np.arange(0,70, 0.5)
        y = model.get_val(x)
        ax.plot(x,y, color='red')
        co2_15 = model.get_val(15)
        co2_60 = model.get_val(60)
        ax.plot(15, co2_15, 'o', color='red')
        ax.plot(60, co2_60, 'o', color='red')


    # linear model
    model = data_fit.models.LinearModel()
    model.fit(np.array(mean_vals), np.array(concentrations))

    if ax is not None:
        x = np.arange(0,70, 0.5)
        y = model.get_val(x)
        ax.plot(x,y, color='red', zorder=0)
        co2_15 = model.get_val(15)
        co2_60 = model.get_val(60)
        #ax.plot(15, co2_15, 'o', color='blue', markersize=3, zorder=10)
        #ax.plot(60, co2_60, 'o', color='blue', markersize=3, zorder=10)



    #ax.plot(x,x*1000,'--',color='black', zorder=-100)




    return model

def get_paper_layout():
    layout = figurefirst.svg_to_axes.FigureLayout(co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_simple) 
    layout.make_mplfigures()
    return layout

def write_layout(layout):
    layout.write_svg(co2_paper_locations.figure_template_locations.figure3_windtunnel_walking_simple)

def plot_co2_concentration():
    layout = get_paper_layout()

    ax = layout.axes[('co2_concentration', 'co2_concentration')]
    plot_co2_series_at_distance_on_ax(ax, distance=2, concentration=5)
    ax.set_ylim(-2,2)
    ax.set_xlim(0,3)
    figurefirst.mpl_functions.adjust_spines(ax, [])#['left', 'bottom'], xticks=xticks, yticks=yticks)

    #ax = layout.axes[('co2_concentration', 'co2_concentration')]
    plot_co2_series_at_distance_on_ax(ax, distance=10, concentration=5)
    ax.set_ylim(-2,2)
    ax.set_xlim(0,3)
    xticks = [0, 1.5, 3]
    yticks = [-2, 0, 2]
    figurefirst.mpl_functions.adjust_spines(ax, ['right', 'bottom'], xticks=xticks, yticks=yticks, spine_locations={'left': 5, 'bottom': 5}, linewidth=0.5)
    ax.set_xticklabels(['0', '1.5', '3'])
    ax.tick_params(length=2.5)

    flytext.set_fontsize(ax.figure, 6)
    layout.append_figure_to_layer(layout.figures['co2_concentration'], 'co2_concentration', cleartarget=True)
    write_layout(layout)

def plot_co2_calibration():
    layout = figurefirst.svg_to_axes.FigureLayout(co2_paper_locations.figure_template_locations.figure3_co2_calibration) 
    layout.make_mplfigures()

    

    def plot_it(ax, xticks, yticks):
        model = plot_concentration_profile(ax)
        
        
        figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], spine_locations={'left': 5, 'bottom': 5}, xticks=xticks, yticks=yticks, linewidth=0.5)
        ytickvals = np.array(yticks)/float(1e6)*100.
        ax.set_yticklabels(ytickvals)

        #layout.svgitems['slope'].style['font-size'] = 6
        #s = str( int(model.parameters['slope']*100)/100. )
        #layout.svgitems['slope'].text = 'slope=' + str(model.parameters['slope'])

        intercept = model.parameters['intercept']
        slope_2cm = model.parameters['slope']
        xs = np.arange(0,70)

        co2_mean, co2_std = load_calibrated_co2_mean_and_std(os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_c5_d10_e0_2017_04_05_14_13_01.bag'))
        slope_10cm = (co2_mean - intercept) / 5.
        ax.plot(5, co2_mean, 'o', color='purple', markersize=3)
        ax.plot(xs, xs*slope_10cm+intercept, color='purple')

        co2_mean, co2_std = load_calibrated_co2_mean_and_std(os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_c1_dcenter_e0_2017_04_05_17_32_12.bag'))
        slope_0cm = (co2_mean - intercept) / 1.
        ax.plot(1, co2_mean, 'o', color='green', markersize=3)
        ax.plot(xs, xs*slope_0cm+intercept, color='green')
        ax.tick_params(length=2.5)

    ax = layout.axes[('co2_calibration', 'co2_calibration')]
    ax.set_ylim(0,30000)
    ax.set_xlim(0,70)
    plot_it(ax, [0,5,15,60], [0, 10000, 20000, 30000])

    ax = layout.axes[('co2_calibration', 'co2_calibration_zoom')]
    ax.set_ylim(0,3000)
    ax.set_xlim(0,6)
    plot_it(ax, [0, 1, 5], [0, 1000, 2000, 3000])

    '''
    co2_mean, co2_std = load_calibrated_co2_mean_and_std(os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_c5_d2_e0_2017_04_05_17_26_34.bag'))
    slope_adjustment = co2_mean / (model.parameters['slope']*5)
    print '******'
    print slope_adjustment
    print '******'
    text_2_cm_5sccm = '5 sccm: ' + str(int(100*(model.parameters['slope']*5/float(1e6)*100.))/100.) + '%'
    text_2_cm_15sccm = '15 sccm: ' + str(int(100*(model.parameters['slope']*15/float(1e6)*100.))/100.) + '%'
    text_2_cm_60sccm = '60 sccm: ' + str(int(100*(model.parameters['slope']*60/float(1e6)*100.))/100.) + '%'

    co2_mean, co2_std = load_calibrated_co2_mean_and_std(os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_c5_d10_e1.5m_2017_04_05_14_28_40.bag'))
    slope_adjustment = co2_mean / (model.parameters['slope']*5)
    text_10_cm_5sccm = str(int(100*(model.parameters['slope']*5*slope_adjustment/float(1e6)*100.))/100.) + '%'
    text_10_cm_15sccm = str(int(100*(model.parameters['slope']*15*slope_adjustment/float(1e6)*100.))/100.) + '%'
    text_10_cm_60sccm = str(int(100*(model.parameters['slope']*60*slope_adjustment/float(1e6)*100.))/100.) + '%'


    co2_mean, co2_std = load_calibrated_co2_mean_and_std(os.path.join(co2_paper_locations.data_locations.co2_calibration_windtunnel, 'windtunnel_c1_dcenter_e0_2017_04_05_17_32_12.bag'))
    slope_adjustment = co2_mean / (model.parameters['slope']*1)
    #layout.svgitems['max_1'].text = '1 sccm: ' + str(int(100*(model.parameters['slope']*1*slope_adjustment/float(1e6)*100.))/100.) + '%'
    text_0_cm_5sccm = '5 sccm: ' + str(int(100*(model.parameters['slope']*5*slope_adjustment/float(1e6)*100.))/100.) + '%'
    text_0_cm_15sccm = '15 sccm: ' + str(int(100*(model.parameters['slope']*15*slope_adjustment/float(1e6)*100.))/100.) + '%'
    text_0_cm_60sccm = '60 sccm: ' + str(int(100*(model.parameters['slope']*60*slope_adjustment/float(1e6)*100.))/100.) + '%'
    '''


    flytext.set_fontsize(ax.figure, 6)
    
    '''
    svgitems = ['slope', 'peak_5', 'peak_15', 'peak_60', 'max_5', 'max_15', 'max_60', 'peak_10_5', 'peak_10_15', 'peak_10_60']
    for svgitem in svgitems:
        layout.svgitems[svgitem].style['font-size'] = 6

    layout.apply_svg_attrs(['slope', 'peak_5', 'peak_15', 'peak_60', 'max_5', 'max_15', 'max_60', 'peak_10_5', 'peak_10_15', 'peak_10_60'])
    '''

    layout.append_figure_to_layer(layout.figures['co2_calibration'], 'co2_calibration', cleartarget=True)
    layout.write_svg(co2_paper_locations.figure_template_locations.figure3_co2_calibration)





if __name__ == '__main__':

    plot_co2_concentration()
    plot_co2_calibration()