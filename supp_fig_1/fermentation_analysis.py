import figurefirst
import numpy as np
import fly_plot_lib.text as flytext
import fly_plot_lib.plot as fpl
import scipy.stats

import co2_paper_locations

fifidatafile = 'figure_data.pickle'

def get_paper_layout():
    svg = co2_paper_locations.figure_template_locations.figure1_ferment_trapassay_flybottles
    #svg = 'simon_foundation.svg'
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout
    
def specific_gravity_to_abv(og, fg):   
    # Ritchie Products Ltd, (Zymurgy, Summer 1995, vol. 18, no. 2) 
    # Michael L. Hall's article Brew by the Numbers: Add Up What's in Your Beer, and Designing Great Beers by Daniels
    ABV =(76.08 * (og-fg) / (1.775-og)) * (fg / 0.794) / 100.
    # in percent / 100 
    return ABV
    
def plot_preference_index(layout, figure, axis, t):
    
    # 2 - day
    control_male = np.array([10,97,60,16,29,31,11,17,17,18,27,35])
    control_female = np.array([9,15,13,17,29,25,7,8,20,24,34,30])
    test_male = np.array([18,75,113,27,43,36,13,23,19,23,31,27])
    test_female = np.array([33,28,42,45,35,46,10,22,18,7,21,23])
    control = control_male + control_female
    test = test_male + test_female
    pi_2day = (test-control) / (test.astype(float)+control.astype(float))
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_box',
                                 layout, figure, axis, fifidatafile, 'pref index 2day',
                                 ['Time (hrs)',
                                  'Pref Index for 2 days'],
                                 48, pi_2day, 
                                 xwidth=16, ywidth=2, color=(0.001,0.001,0.001), edgecolor=(0.001,0.001,0.001), flipxy=False, shading='95conf', markersize=2, linewidth=1, use='mean', optimize_scatter_distance=True)
        
    # 6 - day
    control_male = np.array([46,59,15,4,19,18,57,57,63,55,54,110])
    control_female = np.array([21,38,15,5,10,42,25,31,93,45,35,105])
    test_male = np.array([21,79,26,7,30,26,24,55,63,20,21,66])
    test_female = np.array([14,67,21,4,21,27,17,27,67,19,23,34])
    control = control_male + control_female
    test = test_male + test_female
    pi_6day = (test-control) / (test.astype(float)+control.astype(float))
    #fpl.scatter_box(ax, 168, pi_6day, xwidth=16, ywidth=2, color=(0.001,0.001,0.001), edgecolor=(0.001,0.001,0.001), flipxy=False, shading='95conf', markersize=2, linewidth=1, use='mean', optimize_scatter_distance=True)
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_box',
                                 layout, figure, axis, fifidatafile, 'pref index 6day',
                                ['Time (hrs)',
                                 'Pref Index for 6 days'],
                                 168, pi_6day,
                                 xwidth=16, ywidth=2, color=(0.001,0.001,0.001), edgecolor=(0.001,0.001,0.001), flipxy=False, shading='95conf', markersize=2, linewidth=1, use='mean', optimize_scatter_distance=True)
    

    # 12 - day
    control_male = np.array([17,45,28,21,63,17,53,41,22,70,8,21])
    control_female = np.array([16,29,15,15,27,8,23,41,12,47,10,12])
    test_male = np.array([8,25,36,14,34,15,24,46,9,14,10,15])
    test_female = np.array([17,6,19,10,14,10,6,24,2,5,10,12])
    control = control_male + control_female
    test = test_male + test_female
    pi_12day = (test-control) / (test.astype(float)+control.astype(float))
    #fpl.scatter_box(ax, 288, pi_12day, xwidth=16, ywidth=2, color=(0.001,0.001,0.001), edgecolor=(0.001,0.001,0.001), flipxy=False, shading='95conf', markersize=2, linewidth=1, use='mean', optimize_scatter_distance=True)
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_box',
                                 layout, figure, axis, fifidatafile, 'pref index 12day',
                                ['Time (hrs)',
                                 'Pref Index for 12 days'],
                                 288, pi_12day,
                                 xwidth=16, ywidth=2, color=(0.001,0.001,0.001), edgecolor=(0.001,0.001,0.001), flipxy=False, shading='95conf', markersize=2, linewidth=1, use='mean', optimize_scatter_distance=True)
    


    x = [48]*len(pi_2day) + [168]*len(pi_6day) + [288]*len(pi_12day)
    y = []
    y.extend(pi_2day)
    y.extend(pi_6day)
    y.extend(pi_12day)
    slope, intercept, r, p, stderr = scipy.stats.linregress(x,y)
    
    #ax.plot(np.array(x), np.array(x)*slope+intercept, color='red')
    figurefirst.deprecated_regenerate.mpl('plot',
                                 layout, figure, axis, fifidatafile, 'correlation',
                                ['Time (hrs)',
                                 'Pref Index correlation'],
                                 np.array(x), np.array(x)*slope+intercept, color='red')


    print 'Trap stats: ', 'slope, intercept, r, p, stderr'
    print slope, intercept, r, p, stderr

    #ax.set_ylim(-0.5,0.5)
    figurefirst.deprecated_regenerate.mpl('set_ylim', layout, figure, axis, fifidatafile, 'ylim', [],
                              -0.5,0.5)
    #ax.set_xlim(0,14.5*24)
    figurefirst.deprecated_regenerate.mpl('set_xlim', layout, figure, axis, fifidatafile, 'xlim', [],
                              0,14.5*24)
    xticks = t
    yticks = [-0.5,0,0.5]

    #figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], yticks=yticks, xticks=[2*24,7*24,12*24], spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5)
    #ax.tick_params(length=2.5)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, figure, axis, fifidatafile, 'adjust PI spines', [], 
                                 ['left', 'bottom'], yticks=yticks, xticks=[2*24,7*24,12*24], spine_locations={'left': 5, 'bottom': 5}, linewidth=0.5, tick_length=2.5)


    #ax.set_xticklabels([2,7,12])
    figurefirst.deprecated_regenerate.mpl('set_xticklabels', layout, figure, axis, fifidatafile, 'x labels', [],
                              [2,7,12])

    #flytext.set_fontsize(ax.figure, 6)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, figure, axis, fifidatafile, 'set fontsize', [], 
                                 6)
    #ax.set_yticklabels([])
    
    return r, p

def make_figure():
    

    A = np.array([1.09,  1.088,  1.069, 1.052, 1.04, 1.028, 1.019, 1.012, 1.006, 1.002, 0.999, 0.996, 0.995, 0.994, 0.994])
    r = scipy.stats.uniform(-0.002, 0.002) # note differences betwen for measurments, adding a little noise so that envelope is at least visible in derivative.  
    A1 = A + r.rvs(len(A))
    A2 = A + r.rvs(len(A))
    A3 = A + r.rvs(len(A))
    t = np.array([24*i for i in range(len(A))])
    ABV1 = np.abs([specific_gravity_to_abv( A1[0], a ) for a in A1])
    ABV2 = np.abs([specific_gravity_to_abv( A2[0], a ) for a in A2])
    ABV3 = np.abs([specific_gravity_to_abv( A3[0], a ) for a in A3])

    
    
    
    #alcohol_volume = np.mean(alcohol_volume, axis=0)
    ABV = np.abs([specific_gravity_to_abv( A[0], a ) for a in A])

    volume = 0.130 # L
    alcohol_volume_1 = ABV1*volume
    alcohol_volume_2 = ABV2*volume
    alcohol_volume_3 = ABV3*volume
    alcohol_volume_0 = ABV*volume
    alcohol_volume = np.vstack((alcohol_volume_0, alcohol_volume_1, alcohol_volume_2, alcohol_volume_3))

    #### continue working here

    alcohol_density = 789 # g / L = kg / m^3
    alcohol_molar_mass = 46.068 # g/mol

    alcohol_mols = (alcohol_volume*alcohol_density) / alcohol_molar_mass
    total_CO2_mols = alcohol_mols

    d_alcohol_mols_dt = np.diff(alcohol_mols) / np.diff(t)
    d_CO2_dt = d_alcohol_mols_dt

    mean_volume_of_CO2_per_minute = (np.mean(d_CO2_dt, axis=0) /60.) * 0.08206 * 298 / 1.
    volume_of_CO2_per_minute = (d_CO2_dt /60.) * 0.08206 * 298 / 1.
    volume_of_CO2_per_minute = [v for v in volume_of_CO2_per_minute]

    volume_of_air_per_min = np.array([5., 15., 60., 200.]) # mL / min
    volume_of_air_per_min /= 1000.
    volume_of_air_per_hour = volume_of_air_per_min*60.


    R = 0.08206 # L atm mol-1 K-1
    T = 295. # K
    P = 1. # atm
    # PV = nRT

    mols_air_per_hour = P*volume_of_air_per_hour / (R*T)

    peak_CO2_concentration = np.max( np.mean(d_CO2_dt, axis=0) ) / mols_air_per_hour




    # make figure

    layout = get_paper_layout()

    ax_gravity = layout.axes[('mpl','gravity')]
    ax_co2 = layout.axes[('mpl','co2')]
    ax_xticks = layout.axes[('mpl','xticks')]

    #A_yadj = (np.array(A) - 0.9)
    #ax_gravity.plot(t, A_yadj, color='gray')
    #ax_gravity.plot(t, A_yadj, 'o', markersize=3, markeredgecolor='gray', markerfacecolor='gray')

    #ax_gravity.plot(t, ABV, color='blue')
    figurefirst.deprecated_regenerate.mpl('plot', layout, 'mpl', 'gravity', fifidatafile, 'Alcohol by volume',
                              ['Time',
                               'Alcohol by volume'],
                               t, ABV, color='blue')


    #fpl.scatter_line(ax_gravity, t, [ABV, ABV1, ABV2, ABV3], color='blue', shading='95conf', show_lines=False, use='mean')
    figurefirst.deprecated_regenerate.custom('fly_plot_lib', 'plot.scatter_line',
                                 layout, 'mpl', 'gravity', fifidatafile, 'Alcohol by volume confidence interval',
                                ['Time',
                                 'List of alcohol by volume'],
                                 t, [ABV, ABV1, ABV2, ABV3],
                                 color='blue', shading='95conf', show_lines=False, use='mean')


    #ax_gravity.plot(t, ABV, 'o', markersize=3, markeredgecolor='blue', markerfacecolor='blue')
    figurefirst.deprecated_regenerate.mpl('plot', layout, 'mpl', 'gravity', fifidatafile, 'Alcohol by volume points',
                              [],
                               t, ABV, 'o', markersize=3, markeredgecolor='blue', markerfacecolor='blue')

    t_dt = t[1:] - 12
    #ax_co2.plot(t_dt, mean_volume_of_CO2_per_minute, color='green')
    figurefirst.deprecated_regenerate.mpl('plot', layout, 'mpl', 'co2', fifidatafile, 'CO2',
                              ['Time',
                               'Mean volume of CO2 per minute'],
                               t_dt, mean_volume_of_CO2_per_minute, color='green')

    #ax_co2.plot(t_dt, mean_volume_of_CO2_per_minute, 'o', markersize=3, markeredgecolor='green', markerfacecolor='green')
    figurefirst.deprecated_regenerate.mpl('plot', layout, 'mpl', 'co2', fifidatafile, 'CO2 points',
                              [],
                               t_dt, mean_volume_of_CO2_per_minute, 'o', markersize=3, markeredgecolor='green', markerfacecolor='green')

    #fpl.scatter_line(ax_co2, t_dt, volume_of_CO2_per_minute, color='green', shading='95conf', show_lines=False, use='mean')
    figurefirst.deprecated_regenerate.custom( 'fly_plot_lib', 'plot.scatter_line', layout, 'mpl', 'co2', fifidatafile, 'CO2 confidence interval',
                                  ['Time',
                                   'List of volumes of CO2 per minute'],
                                   t_dt, volume_of_CO2_per_minute, color='green', shading='95conf', show_lines=False, use='mean')

    #ax_gravity.set_ylim(0, 0.2)
    figurefirst.deprecated_regenerate.mpl('set_ylim', layout, 'mpl', 'gravity', fifidatafile, 'ABV ylim', [], 
                              0, 0.2)

    #ax_gravity.set_xlim(0,14.5*24)
    figurefirst.deprecated_regenerate.mpl('set_xlim', layout, 'mpl', 'gravity', fifidatafile, 'ABV xlim', [], 
                              0,14.5*24)
    
    #ax_co2.set_xlim(0,14.5*24)
    figurefirst.deprecated_regenerate.mpl('set_xlim', layout, 'mpl', 'co2', fifidatafile, 'CO2 xlim', [], 
                              0,14.5*24)

    #xticks = [2*24,7*24,12*24]
    xticks = [0*24,5*24,10*24, 15*24]
    yticks = [0, 0.1, 0.2]

    #ax_xticks.set_xlim(0,15)
    figurefirst.deprecated_regenerate.mpl('set_xlim', layout, 'mpl', 'xticks', fifidatafile, 'shared x ticks', [], 
                              0,15)

    #figurefirst.mpl_functions.adjust_spines(ax_xticks, ['bottom'], xticks=xticks, spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5)
    #ax_xticks.tick_params(length=2.5, color='black')
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, 'mpl', 'xticks', fifidatafile, 'adjust x spines', [], 
                              ['bottom'], xticks=xticks, spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5, tick_length=2.5, color='black')

    #ax_xticks.set_xticklabels([0,5,10,15])
    figurefirst.deprecated_regenerate.mpl('set_xticklabels', layout, 'mpl', 'xticks', fifidatafile, 'x tick labels', [], 
                              [0,5,10,15])

    #figurefirst.mpl_functions.adjust_spines(ax_gravity, ['left'], xticks=xticks, yticks=yticks, spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5)
    #ax_gravity.tick_params(length=2.5, color='blue')
    #ax_gravity.spines['left'].set_color('blue')
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, 'mpl', 'gravity', fifidatafile, 'adjust ABV spines', [], 
                              ['left',], xticks=xticks, yticks=yticks, spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5, tick_length=2.5, color='blue')


    #ax_gravity.set_yticklabels([0, 10, 20], color='blue')
    figurefirst.deprecated_regenerate.mpl('set_yticklabels', layout, 'mpl', 'gravity', fifidatafile, 'gravity tick labels', [], 
                              [0, 10, 20], color='blue')

    #ax_co2.set_ylim(0, 0.0012)
    figurefirst.deprecated_regenerate.mpl('set_ylim', layout, 'mpl', 'co2', fifidatafile, 'co2 y lim', [], 
                              0, 0.0012)

    yticks = [0, 0.0006, 0.0012]
    #figurefirst.mpl_functions.adjust_spines(ax_co2, ['right'], spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5, yticks=yticks)
    #ax_co2.tick_params(length=2.5, color='green')
    #ax_co2.spines['right'].set_color('green')
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.adjust_spines', layout, 'mpl', 'co2', fifidatafile, 'adjust CO2 spines', [], 
                              ['right',], xticks=xticks, yticks=yticks, spine_locations={'bottom': 5, 'right': 5, 'left': 5}, linewidth=0.5, tick_length=2.5, color='green')

    #ax_co2.set_yticklabels([str(tick*1000) for tick in yticks], color='green')
    figurefirst.deprecated_regenerate.mpl('set_yticklabels', layout, 'mpl', 'co2', fifidatafile, 'co2 y labels', [], 
                              [str(tick*1000) for tick in yticks], color='green')

    #flytext.set_fontsize(ax_gravity.figure, 6)
    figurefirst.deprecated_regenerate.custom( 'figurefirst', 'mpl_functions.set_fontsize', layout, 'mpl', 'gravity', fifidatafile, 'set fontsize', [], 
                                  6)

    layout.axes[('mpl', 'co2')]._set_fontsize([], 6)
    

    if 0:
        ax_percent_co2 = layout.axes[('mpl', 'percent_co2')]
        ax_percent_co2.plot(volume_of_air_per_min, peak_CO2_concentration, color='black')
        ax_percent_co2.plot(volume_of_air_per_min, peak_CO2_concentration, 'o', markersize=3, markeredgecolor='black', markerfacecolor='black')
        ax_percent_co2.set_ylim(0, 0.25)
        figurefirst.mpl_functions.adjust_spines(ax_percent_co2, ['left', 'bottom'], xticks=volume_of_air_per_min, yticks=[0, 0.05, 0.2])
        yticklabels = ['0%', '5%', '20%']
        xticklabels = [int(v) for v in volume_of_air_per_min*1000]
        #ax_percent_co2.set_xticklabels(xticklabels)
        #ax_percent_co2.set_yticklabels(yticklabels)
        
        ax_percent_co2.set_xticklabels(['' for i in xticklabels])
        ax_percent_co2.set_yticklabels(['' for i in yticklabels])
    
    
    
    r, p = plot_preference_index(layout, 'mpl', 'pref_index', t)
    x = np.linspace(1*24,13*24)
    y = np.zeros_like(x)
    #layout.axes[('mpl','pref_index')].plot(x,y,linestyle='--',color='black',linewidth=1,zorder=-100, dashes=[2,2,2,2])
    figurefirst.deprecated_regenerate.mpl('plot', layout, 'mpl', 'pref_index', fifidatafile, 'pref index dashes', [], 
                              x,y,linestyle='--',color='black',linewidth=1,zorder=-100, dashes=[2,2,2,2])
    
    
    if p<0.001:
        p = 0.001
    regression_text = 'p<' + '{0:.3f}'.format(p) + ', ' + 'r2=' + '{0:.2f}'.format(r**2)
    layout.svgitems['regression'].style['font-size'] = 6
    layout.svgitems['regression'].text = regression_text
    
    
    layout.append_figure_to_layer(layout.figures['mpl'], 'mpl', cleartarget=True)
    
    layout.apply_svg_attrs()
    
    
    layout.write_svg(co2_paper_locations.figure_template_locations.figure1_ferment_trapassay_flybottles )
    #layout.write_svg('simon_foundation.svg' )





if __name__ == '__main__':

    make_figure()