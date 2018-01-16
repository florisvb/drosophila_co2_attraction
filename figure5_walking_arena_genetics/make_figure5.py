
if __name__ == '__main__':

    import make_speed_sorted_figure
    import co2_paper_locations


    # hcs full figure
    make_speed_sorted_figure.make_hcs_full()
    make_speed_sorted_figure.make_hcs_full(odor='ethanol')
    make_speed_sorted_figure.make_hcs_full(odor='vinegar')

    # primary genetics figure
    make_speed_sorted_figure.make_summary_of_genetic_experiments()
    make_speed_sorted_figure.make_summary_of_genetic_experiments('ethanol')
    make_speed_sorted_figure.make_summary_of_genetic_experiments('vinegar')

    # primary figure control
    make_speed_sorted_figure.make_hcs_full(figure='co2-control')
    make_speed_sorted_figure.make_summary_of_genetic_experiments('co2-control')
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='ethanol-control')
    make_speed_sorted_figure.make_summary_of_genetic_experiments(figure='vinegar-control')

    
