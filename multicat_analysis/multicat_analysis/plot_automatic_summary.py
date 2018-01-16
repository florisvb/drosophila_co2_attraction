import plot_n_flies
import plot_circadian
import figurefirst

import numpy as np
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
from optparse import OptionParser



def get_paper_layout():
    svg = 'automatic_summary.svg'
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout
    
def make_automatic_summary(directory):
    layout = get_paper_layout()
    
    plot_circadian.plot_circadian(layout.axes[(u'circadian', u'circadian')], directory)
    layout.append_figure_to_layer(layout.figures['circadian'], 'circadian', cleartarget=True)

    paths = mta.read_hdf5_file_to_pandas.get_filenames(directory, contains='day')
    pd_data_sets = []
    for path in paths:
        pd_filename = mta.read_hdf5_file_to_pandas.get_filename(path, 'pd_data.pickle')
        pd_data = pandas.read_pickle(pd_filename)
        pd_data_sets.append(pd_data)
    pd_data = pandas.concat(pd_data_sets)
    
    flowrates = pd_data.flowrate.unique()
    flowrates.sort()
    flowrate_axes = ['flowrate_a', 'flowrate_b', 'flowrate_c', 'flowrate_d']
    portions_of_day = ['afternoon', 'dusk', 'night', 'morning']
    config = mta.read_hdf5_file_to_pandas.load_config_from_path(paths[0])
    for i, flowrate in enumerate(flowrates):
        for portion_of_day in portions_of_day:
            ax = layout.axes[flowrate_axes[i], portion_of_day]
            localtimerange = config.portion_of_day_to_local_time_range[portion_of_day]
            try:
                plot_n_flies.plot_n_flies_in_control_and_odor_regions(ax, paths, localtimerange, flowrate, average_within_paths=False)
            except:
                print 'FAILED PLOTTING: ', paths, localtimerange, flowrate
        svgitem = 'text_' + flowrate_axes[i]
        layout.svgitems[svgitem].style['font-size'] = 8
        layout.svgitems[svgitem].text = "%0.1f" % flowrate
        layout.append_figure_to_layer(layout.figures[flowrate_axes[i]], flowrate_axes[i], cleartarget=True)
    
    figure_output = os.path.join(directory, 'automatic_figure_output.svg')
    layout.apply_svg_attrs()
    layout.write_svg(figure_output)
    

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("--path", type="str", dest="path", default='',
                        help="path to data")
    (options, args) = parser.parse_args()  
   
    make_automatic_summary(options.path)
