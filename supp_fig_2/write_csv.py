from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'supp_fig_2_walking_data.dillpickle'
    data_filename = regenerate.compress(data_filename, 300)

    header = ("# Source data for Extended Data Figure 2\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")

    panel_id_to_figure =      {'a': ['co2_calibration'],
                               'c': ['h2o_time', 'co2_time', 'eth_time', 'vinegar_time',
                                     'h2o_distances', 'co2_distances', 'eth_distances', 'vinegar_distances',
                                     'h2o_approaches', 'co2_approaches', 'eth_approaches', 'vinegar_approaches',
                                     'h2o_speeds', 'co2_speeds', 'eth_speeds', 'vinegar_speeds'],
                               'd': ['wind_co2'],
                               'e': ['wind_walking_co2'],
                               'f': ['wind_eth'],
                               'g': ['wind_walking_eth']}

    data = regenerate.load_data_file(data_filename)

    panel_id_to_layout_keys = {}
    for panel_id, figure_names in panel_id_to_figure.items():
        layout_keys = []
        for figure_name in figure_names:
            for layout_key in data.keys():
                if unicode(figure_name) == layout_key[0]:
                    layout_keys.append(layout_key)
        panel_id_to_layout_keys[panel_id] = layout_keys

    
    regenerate.write_to_csv(data_filename, 'Extended Data 2', 
                            panel_id_to_layout_keys, header)




