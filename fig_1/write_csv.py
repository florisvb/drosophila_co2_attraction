from figurefirst import regenerate

if __name__ == '__main__':

    header = ("# Source data for Figure 1\n"
              "-----------------------------------------------------------------\n"
              "## Paper:\n"
              "### Distinct activity-gated pathways mediate attraction and aversion to CO2 in Drosophila\n"
              "-----------------------------------------------------------------\n"
              "## Authors: \n"
              "### Floris van Breugel, Ainul Huda, and Michael H Dickinson\n"
              "-----------------------------------------------------------------\n"
              "#####  Data and code: https://github.com/florisvb/drosophila_co2_attraction\n")

    panel_id_to_layout_keys = { 'c': [(u'co2_control', u'heatmap_xz'), (u'co2_control', u'heatmap_xy'), (u'co2_odor', u'heatmap_xz'), (u'co2_odor', u'heatmap_xy')],
                                'd': [(u'eth_control', u'heatmap_xz'), (u'eth_control', u'heatmap_xy'), (u'eth_odor', u'heatmap_xz'), (u'eth_odor', u'heatmap_xy')],
                                'e': [(u'windtunnel_quanitified', u'approached_pad'), (u'windtunnel_quanitified', u'landing_pad'), (u'windtunnel_quanitified', u'approached_blob')],}

    regenerate.write_to_csv('fig_1_freeflight_data.dillpickle', '1', 
                            panel_id_to_layout_keys, header)




