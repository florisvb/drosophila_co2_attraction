if __name__ == '__main__':
    import os
    import co2_paper_locations
    from multicat_analysis import make_speed_sorted_figure
    import figurefirst

    figures = ['single_flies_h_c_s_0', 'single_flies_h_c_s_1', 'hcs_00', 'hcs_04', 'hcs_12', 'hcs_34']
    for figure in figures:
        figurefirst.regenerate.clear_fifidata('supp_fig_4_data.dillpickle', (figure, 'high_speed'))
        figurefirst.regenerate.clear_fifidata('supp_fig_4_data.dillpickle', (figure, 'low_speed'))

    directory = os.path.join(co2_paper_locations.data_locations.walking_arena_symmetric_stim_genetics, 'single_flies_h_c_s')
    for flowrate in [0, 1]:
        make_speed_sorted_figure.make_cross_directory_figure(directory, 
                                                             'single_flies_h_c_s', 
                                                             flowrate=flowrate, 
                                                             use_speed_intercept=make_speed_sorted_figure.USE_SPEED_INTERCEPT, 
                                                             figure='single_concentration_model', 
                                                             add_flowrate_to_label=True)
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='concentration')

    print '#'*125
    print 'To complete remake of Supp Fig 4 run the following notebooks: '
    print 'N_approaches.ipynb'
    print 'model.ipynb'
