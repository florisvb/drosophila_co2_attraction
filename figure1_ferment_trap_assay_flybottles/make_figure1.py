import sys


if __name__ == '__main__':

    sys.path.append('fly_bottle_analysis')
    import analyze_co2_results
    analyze_co2_results.make_plot()
    
    sys.path.append('fermentation')
    import fermentation_analysis
    fermentation_analysis.make_figure()
