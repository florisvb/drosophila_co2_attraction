from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'supp_fig_10_shakefly_data.dillpickle'
    data_filename = regenerate.compress(data_filename, 90)

    header = ("# Source data for Extended Data Figure 8\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")

    panel_id_to_figure = {'a': ['shake']}

    
    data = regenerate.load_data_file(data_filename)

    panel_id_to_layout_keys = {}
    for panel_id, figure_names in panel_id_to_figure.items():
        layout_keys = []
        for figure_name in figure_names:
            for layout_key in data.keys():
                if unicode(figure_name) == layout_key[0]:
                    layout_keys.append(layout_key)
        panel_id_to_layout_keys[panel_id] = layout_keys


    string_replacements = {}


    regenerate.write_to_csv(data_filename, 'Extended Data 10', panel_id_to_layout_keys, header, decimals=3, string_replacements=string_replacements)