import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas
import multi_tracker_analysis as mta
import matplotlib
import fly_plot_lib.plot as fpl
import figurefirst
from optparse import OptionParser

get_filename = mta.read_hdf5_file_to_pandas.get_filename
get_filenames = mta.read_hdf5_file_to_pandas.get_filenames

def get_layout(layout_filename=None):
	if layout_filename is None:
		layout_filename = 'trajectory_plot.svg'
	layout = figurefirst.svg_to_axes.FigureLayout(layout_filename)
	layout.make_mplfigures()
	return layout

def get_color_for_odor_stimulus(config, time_epoch):
	color = []
	for t in time_epoch:
		odor_status = config.get_odor_status(t)
		if odor_status == 'off':
			color.append(-1)
		else:
			color.append(1)
	return np.array(color)

def make_trajectory_plot_for_data_selection_path(path,  destination=None, 
														layout_filename=None, 
														axis='trajectory_plot', 
														figure='trajectory_plot',
														frame_for_odor_stimulus=-1):
	'''
	Loads pickled pandas dataframe found in path, and plots the trajectories as a rasterized scatter plot on top of the img.
	Red = after odor
	Black = before odor 

	'''
	pd, config = mta.read_hdf5_file_to_pandas.load_data_selection_from_path(path)
	layout = get_layout(layout_filename=layout_filename)

	image_sequence_dir = get_filename(path, 'image_sequence')
	bgimg_filenames = get_filenames(os.path.join(path, image_sequence_dir), '.png')
	bgimg = plt.imread(bgimg_filenames[0])

	# equalize bgimgs to match
	r = 0.34/np.mean(bgimg)
	bgimg *= r

	ax = layout.axes[(figure, axis)]

	ax.set_aspect('equal')
	

	#ax.imshow(bgimg, cmap='gray', zorder=1)
	mids = 1- np.abs(bgimg - 0.5)
	darks = 1- np.abs(0.2-bgimg) 
	ax.imshow(bgimg + mids*2 - darks*2, cmap='gray', zorder=-1)

	#ax.plot(pd.position_x.values, pd.position_y.values, '.', color='black', zorder=10)
	color = get_color_for_odor_stimulus(config, pd.time_epoch)
	fpl.scatter(ax, pd.position_x.values, pd.position_y.values, color=color, colormap='seismic', cmap=None, edgecolor='none', radius=1, colornorm=[-1.25,1.25], alpha=1, zorder=10)

	ax.set_ylim(0,bgimg.shape[0])
	ax.set_xlim(0,bgimg.shape[1])

	fpl.adjust_spines(ax, [])
	ax.set_frame_on(False)

	config.draw_mpl(ax, pd.time_epoch.values[frame_for_odor_stimulus] )

	ax.set_rasterization_zorder(1000)

	if destination is None:
		destination = os.path.join(path, 'trajectory_plot.svg')

	layout.append_figure_to_layer(layout.figures[figure], figure, cleartarget=True)
	layout.write_svg(destination)

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('--path', type=str, default='none', help="option: path that points to standard named filename, background image, dvbag, config. If using 'path', no need to provide filename, bgimg, dvbag, and config. Note")
    (options, args) = parser.parse_args()

    make_trajectory_plot_for_data_selection_path(options.path)