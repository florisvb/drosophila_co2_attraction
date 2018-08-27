from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'fig_4_mutants_data.dillpickle'
    data_filename = regenerate.compress(data_filename, 90)

    header = ("# Source data for Figure 4\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")

    panel_id_to_figure = {'a': ['hcs_full'],
                          'b': ['ir64a', 'gr63a', 'double_gr63_ir64'],
                          'c': ['hcs', 'anosmic', 'antennaless', 'orco', 'ir8aM120', 'M37ir25a2', 'ir25a2_BAC', 'M106', 'ir40a'],
                          'd': ['summary_array'],
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
                           'high_speed': 'High Activity (speed > 2.3 mm/sec)',
                           'low_speed': 'Low Activity (speed < 2.3 mm/sec)',
                           }


    regenerate.write_to_csv(data_filename, '4', panel_id_to_layout_keys, header, decimals=3, string_replacements=string_replacements)