from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'supp_fig_5_temp_data.dillpickle'
    data_filename = regenerate.compress(data_filename, 300)

    header = ("# Source data for Extended Figure 5\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")

    panel_id_to_layout_keys = {'a': [('10min_thermistor', '10min_timeseries')],
                               'b': [('10min_thermistor', '10min_co2')],
                               'c': [('10min_thermistor', '10min_air')]}


    string_replacements = {}


    regenerate.write_to_csv(data_filename, '4', panel_id_to_layout_keys, header, decimals=3, string_replacements=string_replacements)