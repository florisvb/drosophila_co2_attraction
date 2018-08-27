import figurefirst

if __name__ == '__main__':

    from multicat_analysis import make_speed_sorted_figure
    import co2_paper_locations
    import os

    figurefirst.regenerate.clear_fifidata('supp_fig_8_eth_vinegar_data.dillpickle', 'all')

    svg = co2_paper_locations.figure_template_locations.supp_fig_8_eth_vinegar
    make_speed_sorted_figure.make_colormaps(svg)

    # hcs full figure
    make_speed_sorted_figure.make_hcs_full(odor='ethanol', figure='ethanol')

    # primary genetics figure
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='ethanol')
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='ethanol-control')

    # hcs full figure
    make_speed_sorted_figure.make_hcs_full(odor='vinegar', figure='vinegar')

    # primary genetics figure
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='vinegar')
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='vinegar-control')


