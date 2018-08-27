import figurefirst

if __name__ == '__main__':

    from multicat_analysis import make_speed_sorted_figure
    import co2_paper_locations
    import os

    figurefirst.regenerate.clear_fifidata('supp_fig_6_co2control_data.dillpickle', 'all')

    svg = co2_paper_locations.figure_template_locations.figure5_walking_arena_genetics_only_CO2_control
    make_speed_sorted_figure.make_colormaps(svg)

    # hcs full figure
    make_speed_sorted_figure.make_hcs_full(figure='co2-control')

    # primary genetics figure
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='co2-control')


