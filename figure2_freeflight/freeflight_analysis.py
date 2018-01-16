import figurefirst
import flydra_pandas_analysis as fpa
import numpy as np
import multi_tracker_analysis as mta
import matplotlib.pyplot as plt
import fly_plot_lib.plot as fpl
import flystat
import fly_plot_lib.text as flytext
import orchard.annotated_time_on_pad_analysis as ann

from co2_paper_locations import data_locations, figure_template_locations

def get_paper_layout():
    layout = figurefirst.svg_to_axes.FigureLayout(figure_template_locations.figure2_freeflight) 
    layout.make_mplfigures()
    return layout
    
def load_co2_datasets(orco=False):
    if not orco:
        pd_co2, config_co2 = fpa.flydra_pandas_dataset.load_pickled_pandas_dataframes_as_one_dataframe(data_locations.freeflight_rotpad_hcs+'/60sccm_co2_dry_hcs')
        pd_co2_odor = fpa.flydra_pandas_dataset.get_pd_for_odor_on(pd_co2, config_co2, start_minute=0, end_minute=30, odor_presentations=None)
        pd_co2_ctrl = fpa.flydra_pandas_dataset.get_pd_for_odor_off(pd_co2, config_co2, start_minute=20, end_minute=50, odor_presentations=None)
        
        datasets = {'pd_co2_odor': pd_co2_odor,
                    'pd_co2_ctrl': pd_co2_ctrl,
                    }
        return datasets
    
    else:
        pd_co2, config_co2 = fpa.flydra_pandas_dataset.load_pickled_pandas_dataframes_as_one_dataframe(data_locations.freeflight_rotpad_orco+'/orco_15percentco2_dry')
        pd_co2_odor = fpa.flydra_pandas_dataset.get_pd_for_odor_on(pd_co2, config_co2, start_minute=0, end_minute=30, odor_presentations=None)
        pd_co2_ctrl = fpa.flydra_pandas_dataset.get_pd_for_odor_off(pd_co2, config_co2, start_minute=20, end_minute=50, odor_presentations=None)
        
        datasets = {'pd_co2_odor': pd_co2_odor,
                    'pd_co2_ctrl': pd_co2_ctrl,
                    }
        return datasets


def load_eth_datasets(orco=False):
    if not orco:
        pd_eth, config_eth = fpa.flydra_pandas_dataset.load_pickled_pandas_dataframes_as_one_dataframe(data_locations.freeflight_rotpad_hcs+'/hcs_60sccm_ethanol_dry')
        pd_eth_odor = fpa.flydra_pandas_dataset.get_pd_for_odor_on(pd_eth, config_eth, start_minute=0, end_minute=30, odor_presentations=None)
        pd_eth_ctrl = fpa.flydra_pandas_dataset.get_pd_for_odor_off(pd_eth, config_eth, start_minute=20, end_minute=50, odor_presentations=None)
        
        datasets = {'pd_eth_odor': pd_eth_odor,
                    'pd_eth_ctrl': pd_eth_ctrl,
                    }
        return datasets
    else:
        pd_eth, config_eth = fpa.flydra_pandas_dataset.load_pickled_pandas_dataframes_as_one_dataframe(data_locations.freeflight_rotpad_orco+'/orco_15sccm_ethanol_dry')
        pd_eth_odor = fpa.flydra_pandas_dataset.get_pd_for_odor_on(pd_eth, config_eth, start_minute=0, end_minute=30, odor_presentations=None)
        pd_eth_ctrl = fpa.flydra_pandas_dataset.get_pd_for_odor_off(pd_eth, config_eth, start_minute=20, end_minute=50, odor_presentations=None)
        
        datasets = {'pd_eth_odor': pd_eth_odor,
                    'pd_eth_ctrl': pd_eth_ctrl,
                    }
        return datasets
    
def __add_approaches_landings_blob_data_to_panel__(label, analysis, odor_state, ax, x):
    data_to_analyze = {'hcs_co2_15sccm' : data_locations.freeflight_rotpad_hcs+'/15percentco2_dry',
                       'hcs_co2_60sccm' : data_locations.freeflight_rotpad_hcs+'/60sccm_co2_dry_hcs',
                       'hcs_eth_15sccm' : data_locations.freeflight_rotpad_hcs+'/hcs_15sccm_ethanol_dry',
                       'hcs_eth_60sccm' : data_locations.freeflight_rotpad_hcs+'/hcs_60sccm_ethanol_dry',
                       'orco_co2_15sccm': data_locations.freeflight_rotpad_orco+'/orco_15percentco2_dry',
                       'orco_eth_15sccm': data_locations.freeflight_rotpad_orco+'/orco_15sccm_ethanol_dry',
                       'orco_eth_60sccm': data_locations.freeflight_rotpad_orco+'/orco_60sccm_ethanol_dry',
                       }
    directory = data_to_analyze[label]
    if odor_state == 'on':
        if 'co2' in directory:
            color = 'magenta'
        elif 'ethanol' in directory:
            color = (0.001, 0.001, 0.001)
        else:
            color = 'green'
    else:
        color = 'blue'
    analysis_on = fpa.approaches_landings.load_stats(directory, analysis=analysis, odor_state=odor_state)
    analysis_means_on = flystat.resampling.bootstrap_medians_from_data(analysis_on, iterations=1000, use='mean')
    confidence_interval_95 = flystat.resampling.bootstrap_confidence_intervals_for_medians(analysis_means_on, confidence_lo=0.025, confidence_hi=0.975)
    fpl.plot_confidence_interval(ax, x, np.mean(analysis_on), confidence_interval_95, confidence_interval_50=None, width=0.2, color=color, linewidth=1, alpha95=0.3, alpha50=0)

def make_approaches_landings_blob_panel(ax, analysis, orco=False):
    if not orco:
        __add_approaches_landings_blob_data_to_panel__('hcs_co2_15sccm', analysis, 'off', ax, -.3)
        __add_approaches_landings_blob_data_to_panel__('hcs_co2_60sccm', analysis, 'off', ax, -.1)
        __add_approaches_landings_blob_data_to_panel__('hcs_co2_15sccm', analysis, 'on', ax, .1)
        __add_approaches_landings_blob_data_to_panel__('hcs_co2_60sccm', analysis, 'on', ax, .3)
        
        __add_approaches_landings_blob_data_to_panel__('hcs_eth_15sccm', analysis, 'off', ax, .7)
        __add_approaches_landings_blob_data_to_panel__('hcs_eth_60sccm', analysis, 'off', ax, .9)
        __add_approaches_landings_blob_data_to_panel__('hcs_eth_15sccm', analysis, 'on', ax, 1.1)
        __add_approaches_landings_blob_data_to_panel__('hcs_eth_60sccm', analysis, 'on', ax, 1.3)
    
    else:
        __add_approaches_landings_blob_data_to_panel__('orco_co2_15sccm', analysis, 'off', ax, -.3)
        __add_approaches_landings_blob_data_to_panel__('orco_co2_15sccm', analysis, 'on', ax, .1)
        
        __add_approaches_landings_blob_data_to_panel__('orco_eth_15sccm', analysis, 'off', ax, .7)
        __add_approaches_landings_blob_data_to_panel__('orco_eth_15sccm', analysis, 'on', ax, 1.1)
    
    xticks = [0, 1]
    yticks = [0, 0.15, 0.3]
    ax.set_xlim(-.5,1.5)
    ax.set_ylim(0,.3)
    if analysis == 'approached':
        figurefirst.mpl_functions.adjust_spines(ax, ['left'], yticks=yticks, spine_locations={'bottom': 5, 'bottom': 5}, linewidth=0.5)
        ax.tick_params(length=2.5)
        ax.set_yticklabels([0,15,30])
        ax.set_xticklabels([])#'CO2', 'EtOH'])
        flytext.set_fontsize(ax.figure, 6)
    else:
        figurefirst.mpl_functions.adjust_spines(ax, [])#['bottom'], xticks=xticks)
        ax.set_xticklabels([])

    
def make_xz_plot(ax, pd, aspect_ratio, draw_boxes=False):
    resolution = 0.003
    vmax = None
    binsx = np.arange(-.5,1.25,resolution)
    binsy = np.arange(-.15,.15,resolution)
    if vmax is None:
        vmax = pd.shape[0]*1e-5
        print 'vmax: ', vmax
    ax = mta.plot.plot_heatmap(pd, binsx, binsy, ax=ax, vmax=vmax, vmin=0, position_x='position_x', position_y='position_z', position_z='position_y', position_z_slice=[-.1,.1])
    xmin = -.15
    xmax = .75
    ax.set_xlim(xmin, xmax)
    x_range = xmax - xmin
    y_range = x_range / aspect_ratio
    print 'xz plot y range: ', y_range
    ax.set_ylim(-1*y_range/2.,y_range/2.)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    ax.set_axis_off()
    
    if draw_boxes:
        from matplotlib import patches
        rect_downwind = patches.Rectangle( (0.06, 0), .1-0.06, 0.03, fill=False, color=(1,.4235,1), linewidth=1)
        rect_blob = patches.Rectangle( (0.3, -.15), .1, .05, fill=False, color='lime', linewidth=1)
        rect_onpad = patches.Rectangle( (-.03, -.01), .07, .01+.003, fill=False, color='cyan', linewidth=1)
        rect_control = patches.Rectangle( (0.31, -.03), .2, .06+.03, fill=False, color='white', linewidth=1)
        
        ax.add_artist(rect_downwind)
        ax.add_artist(rect_blob)
        ax.add_artist(rect_onpad)
        ax.add_artist(rect_control)
    '''
        
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.06,.1], 
                                                                  [-.07, .07], 
                                                                  [0,.03], 
                                                                  region_name='downwind')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [-.03,.04], 
                                                                  [-.01, .066], 
                                                                  [-.003,.01], 
                                                                  region_name='onpad')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.3,.4], 
                                                                  [-.07, .07], 
                                                                  [-.15,-.1], 
                                                                  region_name='blob')
    mta.data_slicing.calc_frames_with_object_in_rectangular_region(  pd, 
                                                                  [.31,.51], 
                                                                  [-.07, .07], 
                                                                  [-.03,.06], 
                                                                  region_name='control')
                                                                  
    '''                                                          

def make_xy_plot(ax, pd, aspect_ratio):
    resolution = 0.003
    vmax = None
    binsx = np.arange(-.5,1.25,resolution)
    binsy = np.arange(-.2,.2,resolution)
    zslice=[-.01, .02]
    if vmax is None:
        vmax = pd.shape[0]*3e-4
        print 'vmax: ', vmax
    ax = mta.plot.plot_heatmap(pd, binsx, binsy, ax=ax, vmax=vmax, vmin=0, position_x='position_x', position_y='position_y', position_z='position_z', position_z_slice=zslice)
    xmin = -.05
    xmax = .05
    ax.set_xlim(xmin, xmax)
    x_range = xmax - xmin
    #y_range = aspect_ratio*x_range 
    y_range = x_range / aspect_ratio
    ax.set_ylim(-1*y_range/2.,y_range/2.)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    ax.set_axis_off()
    
def plot_alicat(ax):
    alicat_bag = data_locations.freeflight_rotpad_hcs+'/hcs_60sccm_ethanol_dry/20151111_164955_Nflydra_alicat.bag'
    
    def fill_half_hour(ax, t, color):
        ax.fill_between([t, t+1800], [70, 70], [0, 0], facecolor=color, edgecolor='none', alpha=0.3)
    
    
    config = fpa.flydra_pandas_dataset.load_config(data_locations.freeflight_rotpad_hcs+'/hcs_60sccm_ethanol_dry' )
    
    for row in config.odor_control:
        if row[1] > 0:
            fill_half_hour(ax, row[0], 'orange')
        elif row[1] == 0:
            fill_half_hour(ax, row[0]+20*60, (0.001, 0.001, 0.001))
    
    fpa.plot_alicat.plot_alicat_data(alicat_bag, ax)
    
    xticks = [0, 3*3600]
    #xticklabels = [-3, 0]
    yticks = [0,60]
    #figurefirst.mpl_functions.adjust_spines(ax, ['left'], xticks=xticks, yticks=yticks)
    figurefirst.mpl_functions.adjust_spines(ax, [])
    #ax.set_xticklabels(xticklabels)
    
    
def plot_xy_and_xz_heatmaps(layout, group, dataset, draw_boxes=False):
    aspect_ratio_xz = layout.axes_groups[group]['heatmap_xz'].w / layout.axes_groups[group]['heatmap_xz'].h
    aspect_ratio_xy = layout.axes_groups[group]['heatmap_xy'].w / layout.axes_groups[group]['heatmap_xy'].h
    make_xz_plot(layout.axes_groups[group]['heatmap_xz']['axis'], dataset, 
                 aspect_ratio_xz, draw_boxes=draw_boxes)
    make_xy_plot(layout.axes_groups[group]['heatmap_xy']['axis'], dataset, 
                 aspect_ratio_xy)
                 
def plot_colorbar():
    layout = get_paper_layout()
    ax = layout.axes_groups['colorbar']['colorbar']
    fpl.colorbar(ax, ticks=None, ticklabels=None, colormap='hot', aspect='auto', orientation='horizontal', filename=None, flipspine=False)
    layout.append_figure_to_layer(layout.figures['colorbar'], 'colorbar', cleartarget=True)
    layout.write_svg(figure_template_locations.figure2_freeflight )
    
    
def print_n():
    co2_datasets = load_co2_datasets()
    print 'CO2'
    print 'N cohorts: ', len(co2_datasets['pd_co2_odor'].dataset_code.unique())
    print 'N trajecs total odor: ', len(co2_datasets['pd_co2_odor'].objid_int.unique())
    print 'N trajecs total NO odor: ', len(co2_datasets['pd_co2_ctrl'].objid_int.unique())
    print 
    
    eth_datasets = load_eth_datasets()
    print 'ETHANOL'
    print 'N cohorts: ', len(eth_datasets['pd_eth_odor'].dataset_code.unique())
    print 'N trajecs total odor: ', len(eth_datasets['pd_eth_odor'].objid_int.unique())
    print 'N trajecs total NO odor: ', len(eth_datasets['pd_eth_ctrl'].objid_int.unique())
    
    ######################
    print
    print 'Quant data: '
    
    
    def print_N_for_quant(label, analysis, odor_state):
        data_to_analyze = {'hcs_co2_15sccm' : data_locations.freeflight_rotpad_hcs+'/15percentco2_dry',
                           'hcs_co2_60sccm' : data_locations.freeflight_rotpad_hcs+'/60sccm_co2_dry_hcs',
                           'hcs_eth_15sccm' : data_locations.freeflight_rotpad_hcs+'/hcs_15sccm_ethanol_dry',
                           'hcs_eth_60sccm' : data_locations.freeflight_rotpad_hcs+'/hcs_60sccm_ethanol_dry',
                           'orco_co2_15sccm': data_locations.freeflight_rotpad_orco+'/orco_15percentco2_dry',
                           'orco_eth_15sccm': data_locations.freeflight_rotpad_orco+'/orco_15sccm_ethanol_dry',
                           'orco_eth_60sccm': data_locations.freeflight_rotpad_orco+'/orco_60sccm_ethanol_dry',
                           }
                           
                           
        directory = data_to_analyze[label]
        analysis_on = fpa.approaches_landings.load_stats(directory, analysis=analysis, odor_state=odor_state)
            
        print
        print label, analysis, odor_state
        print len(analysis_on)
        
    
    for label in ['hcs_co2_15sccm', 'hcs_co2_60sccm', 'hcs_eth_15sccm', 'hcs_eth_60sccm']:
        for analysis in ['approached', 'landed', 'blob']:
            for odor_state in ['off', 'on']:
                print_N_for_quant(label, analysis, odor_state)
    
def make_paper_figure(co2_datasets=None, eth_datasets=None, orco=False):
    layout = get_paper_layout()
    
    
    #plot_alicat(layout.axes['mass_flow_timeseries']['axis'])
    
    if 1:
        if co2_datasets is None:
            co2_datasets = load_co2_datasets(orco=orco)
        if 1:
            plot_xy_and_xz_heatmaps(layout, 'co2_control', co2_datasets['pd_co2_ctrl'], draw_boxes=True)
            plot_xy_and_xz_heatmaps(layout, 'co2_odor', co2_datasets['pd_co2_odor'])
         
    if 1:
        if eth_datasets is None:
            eth_datasets = load_eth_datasets(orco=orco)
        if 1:
            plot_xy_and_xz_heatmaps(layout, 'eth_control', eth_datasets['pd_eth_ctrl'])
            plot_xy_and_xz_heatmaps(layout, 'eth_odor', eth_datasets['pd_eth_odor'])
    
    
    
    
    
    layout.append_figure_to_layer(layout.figures['co2_control'], 'co2_control', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['co2_odor'], 'co2_odor', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['eth_control'], 'eth_control', cleartarget=True)
    layout.append_figure_to_layer(layout.figures['eth_odor'], 'eth_odor', cleartarget=True)
    
    layout.write_svg(figure_template_locations.figure2_freeflight )

def make_statistics_plot(orco=False):
    layout = get_paper_layout()
    make_approaches_landings_blob_panel(layout.axes[('windtunnel_quanitified', 'approached_pad')]['axis'], 'approached', orco=orco)
    make_approaches_landings_blob_panel(layout.axes[('windtunnel_quanitified', 'landing_pad')]['axis'], 'landed', orco=orco)
    make_approaches_landings_blob_panel(layout.axes[('windtunnel_quanitified', 'approached_blob')]['axis'], 'blob', orco=orco)
    layout.append_figure_to_layer(layout.figures['windtunnel_quanitified'], 'windtunnel_quanitified', cleartarget=True)
    layout.write_svg(figure_template_locations.figure2_freeflight )
    

if __name__ == '__main__':
    
    make_paper_figure()
    
    plot_colorbar()
