import sys
import figurefirst

if __name__ == '__main__':

    figurefirst.regenerate.clear_fifidata('supp_fig_1_trapassay_data.dillpickle', 'all')
    
    # confirmed
    sys.path.append('fly_bottle_analysis')
    import analyze_co2_results
    analyze_co2_results.make_plot()
    analyze_co2_results.plot_calibration_supplemental_figure()
    
    # confirmed
    sys.path.append('fermentation')
    import fermentation_analysis
    fermentation_analysis.make_figure()
