import os

# for figure 1
fly_bottle_data = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/fly_bottle_co2_measurements/bottles'
fly_bottle_calibration = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/fly_bottle_co2_measurements/calibrations'

# for figure 2
freeflight_rotpad_hcs = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/windtunnel_freeflight_rotpad_hcs'
freeflight_rotpad_orco = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/windtunnel_freeflight_rotpad_orco'

# primarily for figure3
basepath_for_windtunnel_walking = '/media/caveman/CO2_hdf5_uw/singlepad_windtunnel_blacktop/'
windtunnel_walking = {#  'orco_co260': basepath_for_windtunnl_walking+'60sccm_co2_ORCO',
                         'co2_5sccm':  basepath_for_windtunnel_walking+'CO2/5sccm_CO2',
                         'co2_60sccm': basepath_for_windtunnel_walking+'CO2/60sccm_co2_mag',
                         'eth_60sccm': basepath_for_windtunnel_walking+'ETHANOL/60sccm_ethanol',
                         'co2_200sccm': basepath_for_windtunnel_walking+'CO2/200sccm_co2_mag',
                         #'co2_hungry': basepath_for_windtunnel_walking+'co2_hungry_60sccm',
                         'eth_200sccm': basepath_for_windtunnel_walking+'ETHANOL/200sccm_ethanol',
                         'vinegar_15sccm': basepath_for_windtunnel_walking+'VINEGAR/15sccm_vinegar_hcs',
                         'vinegar_60sccm': basepath_for_windtunnel_walking+'VINEGAR/60sccm_vinegar',
                         'vinegar_200sccm': basepath_for_windtunnel_walking+'VINEGAR/200sccm_vinegar',
                         #'vinegar_60sccm_15co2': basepath_for_windtunnel_walking+'60sccm_vinegar_15sccm_co2',
                         'co2_15sccm': basepath_for_windtunnel_walking+'CO2/co2_15sccm',
                         'h2o_60sccm': basepath_for_windtunnel_walking+'H2O/h20_60sccm',
                         'ethylacetate_15sccm': basepath_for_windtunnel_walking+'ETHYLACETATE/ethylacetate_200ul_100ml_h20_15sccm',
                         'ethylacetate_15_50': basepath_for_windtunnel_walking+'ETHYLACETATE/ethacet_15sccm_50ul_in100ml_h2o',
                         'eth_15sccm': basepath_for_windtunnel_walking+'ETHANOL/ethanol_15sccm',
                         'ethylacetate_60sccm':  basepath_for_windtunnel_walking+'ETHYLACETATE/ethylacetate_200ul_100ml_h20_60sccm',
                         'eth-co2_60_15': basepath_for_windtunnel_walking+'MIXTURES/ethanol_60sccm_15sccm_co2',
                         'h2o_200sccm': basepath_for_windtunnel_walking+'H2O/h2o_200sccm',
                         'h2o_15sccm': basepath_for_windtunnel_walking+'H2O/h2o_15sccm',
                         }
                         
# primarily for figure 4
basepath_arena_activity = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/'
walking_arena_activity = {  'control_afternoon': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_air_slow_delayed_ssr_lowflow_constant_flowrate_change', 
                            'control_dusk': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_air_slow_delayed_ssr_lowflow_constant_flowrate_change',
                            'control_night': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_air_slow_delayed_ssr_lowflow_constant_flowrate_change', 
                            'control_morning': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_air_slow_delayed_ssr_lowflow_constant_flowrate_change',
                            
                            'low_flow_co2_afternoon': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1sccm', 
                            'low_flow_co2_dusk': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1sccm',
                            'low_flow_co2_night': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1sccm', 
                            'low_flow_co2_morning': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1sccm',
                            
                            'example': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1sccm',

                            'balanced_flow_0': basepath_arena_activity+'low_flow_walking_arena_B/hcs/symmetric_0_1_hcs_noonstarved',
                            'balanced_flow_1': basepath_arena_activity+'low_flow_walking_arena_B/hcs/symmetric_0_1_hcs_noonstarved',

                            'singlepulse_2hrs_nflies': basepath_arena_activity+'low_flow_walking_arena_B_short_pulses/2hr_exps',
                            #'singlepulse_2min_nflies': '/media/Orchard4/tmaze_style_exps/short_exps',
                            'singlepulse_2min_nflies':  basepath_arena_activity+'low_flow_walking_arena_B_short_pulses/10min_exps',

                            'singlepulse_2hrs': basepath_arena_activity+'low_flow_walking_arena_B_short_pulses/'+'2hr_exps',
                            'singlepulse_2min':  basepath_arena_activity+'low_flow_walking_arena_B_short_pulses/'+'10min_exps',


                            'low_flow_co2_starved_afternoon': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1ssm_24hrstarved', 
                            'low_flow_co2_starved_dusk': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1ssm_24hrstarved',
                            'low_flow_co2_starved_night': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1ssm_24hrstarved', 
                            'low_flow_co2_starved_morning': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_co2_1ssm_24hrstarved',
                            
                            'ethanol_afternoon': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_ethanol_1sccm', 
                            'ethanol_dusk': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_ethanol_1sccm',
                            'ethanol_night': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_ethanol_1sccm', 
                            'ethanol_morning': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_ethanol_1sccm',
                            
                            'ethanol_starved_afternoon': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_ethanol_1sccm_24hrstarved', 
                            'ethanol_starved_dusk': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_ethanol_1sccm_24hrstarved',
                            'ethanol_starved_night': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_ethanol_1sccm_24hrstarved', 
                            'ethanol_starved_morning': basepath_arena_activity+'low_flow_walking_arena_A/'+'multicat_ethanol_1sccm_24hrstarved',
                            
                            'low_flow_warm_dusk': basepath_arena_activity+'low_flow_walking_arena_A/multicat_co2_ssr_delayed_lowflow_1sccmco2_hot',
                            'high_flow_co2_dusk': basepath_arena_activity+'high_flow_walking_arena/'+'holyminiolfacto_co2_5percent_HCS_noonstarved_noonstarted',
                            'high_flow_warm_dusk': basepath_arena_activity+'high_flow_walking_arena/'+'hot_100sccm',
                            'high_flow_arista_dusk': basepath_arena_activity+'high_flow_walking_arena/'+'aristae',
                            
                            'high_flow_ethanol_dusk': basepath_arena_activity+'high_flow_walking_arena/'+'eth_concentrations_100sccm_1_5_15',
                            }

# for figure 5
walking_arena_hcs_concentration = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/low_flow_walking_arena_B_HCS_concentration'
walking_arena_symmetric_stim_genetics = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/low_flow_walking_arena_B'
walking_arena_symmetric_stim_genetics_ethanol = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/low_flow_walking_arena_B_ETHANOL'
walking_arena_symmetric_stim_genetics_vinegar = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/low_flow_walking_arena_B_VINEGAR'
walking_arena_symmetric_stim_genetics_ethanol_hcs = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/low_flow_walking_arena_B_ETHANOL/eth_hcs'
walking_arena_symmetric_stim_genetics_vinegar_hcs = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/low_flow_walking_arena_B_VINEGAR/vinegar_hcs'
walking_arena_tripleparadigm_hcs = '/media/caveman/CO2_hdf5_uw/CO2_data_notrajecs/low_flow_walking_arena_B/hcs'

# for figure 3
co2_calibration_windtunnel = '/home/caveman/Dropbox_flyranch/2017_CO2/drosophila_co2_attraction/fig_2/co2_concentration_wind_tunnel' 

# for supplemental figure
shake_fly_mosquito_data = '/home/caveman/Dropbox_flyranch/2017_CO2/drosophila_co2_attraction/supp_fig_10/shake_fly_mosquito_data'