import figurefirst

if __name__ == '__main__':

    figurefirst.regenerate.clear_fifidata('supp_fig_2_walking_data.dillpickle', 'all')

    # confirmed
    from multicat_analysis import plot_co2_concentration
    plot_co2_concentration.plot_co2_calibration() # this is the supplement

    # confirmed
    from multicat_analysis import make_time_on_pad_plots_multiconcentration
    data_dict_lengths = make_time_on_pad_plots_multiconcentration.ann.load_lengths()
    data_dict_entries = make_time_on_pad_plots_multiconcentration.ann.load_lengths(load_data='odor_entries')
    data_dict_distances = make_time_on_pad_plots_multiconcentration.ann.load_lengths(load_data='distances')
    data_dict_speeds = make_time_on_pad_plots_multiconcentration.ann.load_lengths(load_data='speeds')
    data = [data_dict_lengths, data_dict_entries, data_dict_distances, data_dict_speeds]

    make_time_on_pad_plots_multiconcentration.plot_lengths_supplemental('h2o', data)
    make_time_on_pad_plots_multiconcentration.plot_lengths_supplemental('co2', data)
    make_time_on_pad_plots_multiconcentration.plot_lengths_supplemental('eth', data)
    make_time_on_pad_plots_multiconcentration.plot_lengths_supplemental('vinegar', data)

    from multicat_analysis import windtunnel_pad_compared_to_walking_data
    windtunnel_pad_compared_to_walking_data.update_all_and_save()   