
if __name__ == '__main__':

    import make_summary_figure

    make_summary_figure.make_single_pulse_plots()

    make_summary_figure.plot_circadian()
    make_summary_figure.plot_attraction_index_vs_speed()

    high_flow_panels = ['high_flow_ethanol_dusk', 'high_flow_warm_dusk', 'high_flow_arista_dusk', 'high_flow_co2_dusk']
    make_summary_figure.update_panels(high_flow_panels)

    if 1:
        panel_basenames = ['control', 'low_flow_co2', 'low_flow_co2_starved', 'ethanol_starved']
        for panel_basename in panel_basenames:
            make_summary_figure.update_panels_for_full_day(panel_basename)

    import windtunnel_pad_compared_to_walking_data
    windtunnel_pad_compared_to_walking_data.update_all_and_save()