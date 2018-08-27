import imp
import figurefirst

if __name__ == '__main__':

    figurefirst.regenerate.clear_fifidata('fig_1_freeflight_data.dillpickle', 'all')

    # confirmed
    freeflight_analysis = imp.load_source('freeflight_analysis', 'freeflight_analysis.py')
    freeflight_analysis.make_paper_figure()
    freeflight_analysis.plot_colorbar()

    # confirmed 
    freeflight_analysis.make_statistics_plot()

    
