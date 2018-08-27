from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'fig_2_walking_data.dillpickle'
    # data_filename = regenerate.compress(data_filename, 90)

    header = ("# Source data for Figure 2\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")
    panel_id_to_layout_keys = { 'c': [(u'co2_concentration', u'co2_concentration')],
                                'd': [(u'trajectories', u'h2o'), (u'trajectories', u'co2'), (u'trajectories', u'ethanol')],
                                'e': [(u'distances', u'h2o_60sccm'), (u'distances', u'co2_60sccm'), (u'distances', u'eth_60sccm'), (u'distances', u'eth-co2_60_15'), (u'distances', u'vinegar_60sccm'), 
                                (u'time', u'h2o_60sccm'), (u'time', u'co2_60sccm'), (u'time', u'eth_60sccm'), (u'time', u'eth-co2_60_15'), (u'time', u'vinegar_60sccm'), 
                                (u'speeds', u'h2o_60sccm'), (u'speeds', u'co2_60sccm'), (u'speeds', u'eth_60sccm'), (u'speeds', u'eth-co2_60_15'), (u'speeds', u'vinegar_60sccm'), 
                                (u'approaches', u'h2o_60sccm'), (u'approaches', u'co2_60sccm'), (u'approaches', u'eth_60sccm'), (u'approaches', u'eth-co2_60_15'), (u'approaches', u'vinegar_60sccm')],
                                }

    string_replacements = {'singlepulse_2min': 'singlepulse_10min'}

    regenerate.write_to_csv(data_filename, '2', 
                            panel_id_to_layout_keys, header, decimals=3, string_replacements=string_replacements)

