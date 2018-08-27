from figurefirst import regenerate

if __name__ == '__main__':

    data_filename = 'supp_fig_1_trapassay_data.dillpickle'
    data_filename = regenerate.compress(data_filename, 90)

    header = ("# Source data for Extended Data Figure 1\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")

    panel_id_to_layout_keys = {'a': [(u'mpl', u'xticks'), (u'mpl', u'gravity'), (u'mpl', u'co2')],
                               'c': [(u'mpl', u'pref_index')],
                               'd': [(u'bottles', u'fly_bottles')],
                               'f': [(u'co2_calibration', u'calibration_timecourse')],
                               'g': [(u'co2_calibration', u'calibration_curve')],
                                }

    
    regenerate.write_to_csv(data_filename, 'Extended Data 1', 
                            panel_id_to_layout_keys, header)




