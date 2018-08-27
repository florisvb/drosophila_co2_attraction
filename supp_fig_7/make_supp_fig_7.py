import figurefirst

if __name__ == '__main__':

    from multicat_analysis import make_speed_sorted_figure
    import co2_paper_locations
    import os

    figurefirst.regenerate.clear_fifidata('supp_fig_7_ir25a_ir40a_data.dillpickle', 'all')

    # primary genetics figure
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='ir25a_ir40a-co2')
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='ir25a_ir40a-control')

