import imp


if __name__ == '__main__':

    freeflight_analysis = imp.load_source('freeflight_analysis', 'freeflight_analysis.py')
    freeflight_analysis.make_paper_figure()


    freeflight_analysis.make_statistics_plot()

    freeflight_analysis.plot_colorbar()
