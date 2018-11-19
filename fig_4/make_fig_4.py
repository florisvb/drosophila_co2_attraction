import figurefirst

if __name__ == '__main__':

    from multicat_analysis import make_speed_sorted_figure
    import co2_paper_locations
    import os
    import multicat_analysis

    if 1:
        figurefirst.regenerate.clear_fifidata('fig_4_mutants_data.dillpickle', 'all')

        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2
        make_speed_sorted_figure.make_colormaps(svg)

        # hcs full figure
        make_speed_sorted_figure.make_hcs_full()

        # primary genetics figure
        make_speed_sorted_figure.make_summary_of_genetic_experiments()


        
        multicat_analysis.plot_genetic_summary_table.plot_attraction_aversion_array_on_layout(odor='co2')
        svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2
        layout = figurefirst.svg_to_axes.FigureLayout(svg)
        figurefirst.regenerate.clear_fifidata('figure5_speed_sorted_walking_genetics_only_CO2_new_data.dillpickle', 'Supplemental Data')
        multicat_analysis.plot_genetic_summary_table.save_statistical_data(layout)
   
    #multicat_analysis.plot_genetic_summary_table.write_in_stats(odor='co2', bonferoni_N=11)
    #multicat_analysis.plot_genetic_summary_table.write_in_values(odor='co2')

