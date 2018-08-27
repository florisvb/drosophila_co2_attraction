from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'supp_fig_7_ir25a_ir40a_data.dillpickle'
    data_filename = regenerate.compress(data_filename, 90)

    header = ("# Source data for Extended Data Figure 7\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")

    panel_id_to_figure = {'a': ['ir25a_nothot', 'ir25a_hot', 'M37ir25a2_nothot', 'M37ir25a2_hot', 'ir25a2_BAC_nothot', 'ir25a2_BAC_hot'],
                          'b': ['control_ir25a_nothot', 'control_ir25a_hot', 'control_M37ir25a2_nothot', 'control_M37ir25a2_hot', 'control_ir25a2_BAC_nothot', 'control_ir25a2_BAC_hot'],
                          'c': ['ir40a_nothot',  'ir40a_hot'],
                          'd': ['control_ir40a_nothot',  'control_ir40a_hot'],
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
                           'M37ir25a2': 'ir25a2',
                           'double_gr63_ir64': 'gr63a_ir64a',
                           'ir8aM120': 'ir8a',
                           'hcs_full': 'hcs',
                           'nflies': 'Preference index',
                           'high_speed': 'High Activity (speed > 2.3 mm/sec)',
                           'low_speed': 'Low Activity (speed < 2.3 mm/sec)',
                           }


    regenerate.write_to_csv(data_filename, 'Extended Data Figure 7', panel_id_to_layout_keys, header, decimals=3, string_replacements=string_replacements)