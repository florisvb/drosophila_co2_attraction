from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'supp_fig_8_eth_vinegar_data.dillpickle'
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

    panel_id_to_figure = {'a': ['eth_hcs_full'],
                          'b': ['eth_hcs', 'eth_anosmic', 'eth_antennaless', 'eth_orco', 'eth_ir8aM120', 'eth_M37ir25a2', 'eth_ir25a2_BAC', 'eth_M106'],
                          'c': ['control_eth_hcs', 'control_eth_anosmic', 'control_eth_antennaless', 'control_eth_orco', 'control_eth_ir8aM120', 'control_eth_M37ir25a2', 'control_eth_ir25a2_BAC', 'control_eth_M106'],
                          'd': ['eth_hcs_full'],
                          'e': ['vinegar_hcs', 'vinegar_anosmic', 'vinegar_antennaless', 'vinegar_orco', 'vinegar_ir8aM120', 'vinegar_M37ir25a2', 'vinegar_ir25a2_BAC', 'vinegar_M106'],
                          'f': ['control_vinegar_hcs', 'control_vinegar_anosmic', 'control_vinegar_orco', 'control_vinegar_ir8aM120', 'control_vinegar_M37ir25a2', 'control_vinegar_M106'],
                          }

    
    data = regenerate.load_data_file(data_filename)

    panel_id_to_layout_keys = {}
    for panel_id, figure_names in panel_id_to_figure.items():
        layout_keys = []
        for figure_name in figure_names:
            for layout_key in data.keys():
                if unicode(figure_name) == layout_key[0]:
                    layout_keys.append(layout_key)
        panel_id_to_layout_keys[panel_id] = layout_keys


    string_replacements = {'M106': 'orco_ir8a_gr63a',
                           'M37ir25a2': 'ir25a',
                           'double_gr63_ir64': 'gr63a_ir64a',
                           'ir25a2_BAC': 'ir25a_BAC',
                           'ir8aM120': 'ir8a',
                           'hcs_full': 'hcs',
                           'nflies': 'Preference index',
                           'high_speed': 'High Activity (speed > 1.0 mm/sec)',
                           }


    regenerate.write_to_csv(data_filename, 'Extended Data 8', panel_id_to_layout_keys, header, decimals=3, string_replacements=string_replacements)