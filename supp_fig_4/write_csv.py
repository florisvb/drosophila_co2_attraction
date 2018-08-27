from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'supp_fig_4_data.dillpickle'
    data_filename = regenerate.compress(data_filename, 90)

    header = ("# Source data for Extended Data Figure 4\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")

    panel_id_to_figure = {'a': ['single_flies_h_c_s_0', 'single_flies_h_c_s_1'],
                          'b': ['hcs_00', 'hcs_04', 'hcs_12', 'hcs_35'],
                          'd': ['model'],
                          'e': ['n_approaches_stats'],
                          'f': ['n_approaches_00', 'n_approaches_04', 'n_approaches_12', 'n_approaches_35']}

    
    data = regenerate.load_data_file(data_filename)

    panel_id_to_layout_keys = {}
    for panel_id, figure_names in panel_id_to_figure.items():
        layout_keys = []
        for figure_name in figure_names:
            for layout_key in data.keys():
                if unicode(figure_name) == layout_key[0]:
                    layout_keys.append(layout_key)
        panel_id_to_layout_keys[panel_id] = layout_keys


    string_replacements = {'h_c_s_0': 'HCS_control',
                           'h_c_s_1': 'HCS_CO2',
                           '_00': ' 0 percent CO2',
                           '_04': ' 1.7 percent CO2',
                           '_12': ' 5 percent CO2',
                           '_35': ' 15 percent CO2',
                           }


    regenerate.write_to_csv(data_filename, 'Extended Data 4', panel_id_to_layout_keys, header, decimals=3, string_replacements=string_replacements)