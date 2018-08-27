import figurefirst

if __name__ == '__main__':

    figurefirst.regenerate.clear_fifidata('fig_2_walking_data.dillpickle', 'all')

    # confirmed
    from multicat_analysis import make_time_on_pad_plots_multiconcentration
    make_time_on_pad_plots_multiconcentration.plot_all_lengths(version='new')

    # confirmed
    from multicat_analysis import plot_co2_concentration
    plot_co2_concentration.plot_co2_concentration()

    # confirmed
    from multicat_analysis import plot_trajectories
    plot_trajectories.plot_all_trajectories_on_figure(version='new')

