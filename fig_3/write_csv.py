from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'fig_3_activity_data.dillpickle'
    data_filename = regenerate.compress(data_filename, 90)


    header = ("# Source data for Figure 3\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")
    
    panel_id_to_figure = {'b': ['singlepulse_2min', 'singlepulse_2min_nflies'],
                          'c': ['singlepulse_2hrs', 'singlepulse_2hrs_nflies'],
                          'd': ['control_circadian_colors', 'control_circadian', 'control_afternoon', 'control_dusk', 'control_night', 'control_morning',
                                'hcs_circadian_colors', 'hcs_circadian', 'hcs_afternoon', 'hcs_dusk', 'hcs_night', 'hcs_morning',
                                'starved_circadian_colors', 'starved_circadian', 'low_flow_co2_starved_afternoon', 'low_flow_co2_starved_dusk', 'low_flow_co2_starved_night', 'low_flow_co2_starved_morning',
                                'ethanol_circadian_colors', 'ethanol_circadian', 'ethanol_starved_afternoon', 'ethanol_starved_dusk', 'ethanol_starved_night', 'ethanol_starved_morning',],
                          'e': ['high_flow_co2_dusk', 'high_flow_arista_dusk', 'high_flow_warm_dusk', 'high_flow_ethanol_dusk'],
                          'f': ['activity_correlation_relative'],
                          }

    data = regenerate.load_data_file(data_filename)

    panel_id_to_layout_keys = {}
    for panel_id, figure_names in panel_id_to_figure.items():
        layout_keys = []
        for figure_name in figure_names:
            for layout_key in data.keys():
                if unicode(figure_name) in layout_key[0]:
                    layout_keys.append(layout_key)
        panel_id_to_layout_keys[panel_id] = layout_keys


    regenerate.write_to_csv(data_filename, '3', panel_id_to_layout_keys, header)