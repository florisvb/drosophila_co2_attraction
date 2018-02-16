

if __name__ == '__main__':

    #import make_time_on_pad_plots_multiconcentration
    #make_time_on_pad_plots_multiconcentration.plot_all_lengths(version='new')

    #import plot_co2_concentration
    #plot_co2_concentration.plot_co2_concentration()
    #plot_co2_concentration.plot_co2_calibration()

    import plot_trajectories
    plot_trajectories.plot_all_trajectories_on_figure(version='new')


    if 0:
	    data_dict_lengths = make_time_on_pad_plots_multiconcentration.ann.load_lengths()
	    data_dict_entries = make_time_on_pad_plots_multiconcentration.ann.load_lengths(load_data='odor_entries')
	    data_dict_distances = make_time_on_pad_plots_multiconcentration.ann.load_lengths(load_data='distances')
	    data_dict_speeds = make_time_on_pad_plots_multiconcentration.ann.load_lengths(load_data='speeds')
	    data = [data_dict_lengths, data_dict_entries, data_dict_distances, data_dict_speeds]

	    make_time_on_pad_plots_multiconcentration.plot_lengths_supplemental('h2o', data)
	    make_time_on_pad_plots_multiconcentration.plot_lengths_supplemental('co2', data)
	    make_time_on_pad_plots_multiconcentration.plot_lengths_supplemental('eth', data)
	    make_time_on_pad_plots_multiconcentration.plot_lengths_supplemental('vinegar', data)