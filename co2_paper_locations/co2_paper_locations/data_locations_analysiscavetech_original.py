import os


# for figure 1
fly_bottle_data = '/media/Orchard3/fly_bottle_data/bottles'
fly_bottle_calibration = '/media/Orchard3/fly_bottle_data/calibrations'

# for figure 2
freeflight_rotpad_hcs = '/media/caveman/Orchard/rotpad'
freeflight_rotpad_orco = '/media/Orchard2/rotpad2'

# primarily for figure3
windtunnel_walking = {#  'orco_co260': '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_co2_ORCO',
                         'co2_5sccm':  '/media/Orchard2/singlepad_windtunnel_blacktop/5sccm_CO2',
                         'co2_60sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_co2_mag',
                         'eth_60sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_ethanol',
                         'co2_200sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/200sccm_co2_mag',
                         #'co2_hungry': '/media/Orchard2/singlepad_windtunnel_blacktop/co2_hungry_60sccm',
                         'eth_200sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/200sccm_ethanol',
                         'vinegar_15sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/15sccm_vinegar_hcs',
                         'vinegar_60sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_vinegar',
                         'vinegar_200sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/200sccm_vinegar',
                         'vinegar_60sccm_15co2': '/media/Orchard2/singlepad_windtunnel_blacktop/60sccm_vinegar_15sccm_co2',
                         'co2_15sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/co2_15sccm',
                         'h2o_60sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/h20_60sccm',
                         'ethylacetate_15sccm': '/media/Orchard2/singlepad_windtunnel_blacktop/ethylacetate_200ul_100ml_h20_15sccm',
                         'ethylacetate_15_50': '/media/Orchard2/singlepad_windtunnel_blacktop/ethacet_15sccm_50ul_in100ml_h2o',
                         'eth_15sccm': '/media/Orchard3/more_single_pad_blacktop_windtunnel/ethanol_15sccm',
                         'ethylacetate_60sccm':  '/media/Orchard2/singlepad_windtunnel_blacktop/ethylacetate_200ul_100ml_h20_60sccm',
                         'eth-co2_60_15': '/media/Orchard3/more_single_pad_blacktop_windtunnel/ethanol_60sccm_15sccm_co2',
                         'h2o_200sccm': '/media/Orchard3/more_single_pad_blacktop_windtunnel/h2o_200sccm',
                         'h2o_15sccm': '/media/Orchard3/more_single_pad_blacktop_windtunnel/h2o_15sccm',
                         }
                         
# primarily for figure 4
walking_arena_activity = {  'control_afternoon': '/media/Orchard3/CONTROLED_5massflow/multicat_air_slow_delayed_ssr_lowflow_constant_flowrate_change', 
                            'control_dusk': '/media/Orchard3/CONTROLED_5massflow/multicat_air_slow_delayed_ssr_lowflow_constant_flowrate_change',
                            'control_night': '/media/Orchard3/CONTROLED_5massflow/multicat_air_slow_delayed_ssr_lowflow_constant_flowrate_change', 
                            'control_morning': '/media/Orchard3/CONTROLED_5massflow/multicat_air_slow_delayed_ssr_lowflow_constant_flowrate_change',
                            
                            'low_flow_co2_afternoon': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1sccm', 
                            'low_flow_co2_dusk': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1sccm',
                            'low_flow_co2_night': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1sccm', 
                            'low_flow_co2_morning': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1sccm',
                            
                            'example': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1sccm',

                            'balanced_flow_0': '/media/Orchard4/SYMMETRIC_stim/hcs/symmetric_0_1_hcs_noonstarved',
                            'balanced_flow_1': '/media/Orchard4/SYMMETRIC_stim/hcs/symmetric_0_1_hcs_noonstarved',

                            'singlepulse_2hrs_nflies': '/media/Orchard4/tmaze_style_exps/2hr_exps',
                            #'singlepulse_2min_nflies': '/media/Orchard4/tmaze_style_exps/short_exps',
                            'singlepulse_2min_nflies':  '/media/Orchard4/tmaze_style_exps/10min_exps',

                            'singlepulse_2hrs': '/media/Orchard4/tmaze_style_exps/2hr_exps',
                            #'singlepulse_2min': '/media/Orchard4/tmaze_style_exps/short_exps',
                            'singlepulse_2min':  '/media/Orchard4/tmaze_style_exps/10min_exps',


                            'low_flow_co2_starved_afternoon': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1ssm_24hrstarved', 
                            'low_flow_co2_starved_dusk': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1ssm_24hrstarved',
                            'low_flow_co2_starved_night': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1ssm_24hrstarved', 
                            'low_flow_co2_starved_morning': '/media/Orchard3/CONTROLED_5massflow/multicat_co2_1ssm_24hrstarved',
                            
                            'ethanol_afternoon': '/media/Orchard3/CONTROLED_5massflow/multicat_ethanol_1sccm', 
                            'ethanol_dusk': '/media/Orchard3/CONTROLED_5massflow/multicat_ethanol_1sccm',
                            'ethanol_night': '/media/Orchard3/CONTROLED_5massflow/multicat_ethanol_1sccm', 
                            'ethanol_morning': '/media/Orchard3/CONTROLED_5massflow/multicat_ethanol_1sccm',
                            
                            'ethanol_starved_afternoon': '/media/Orchard3/CONTROLED_5massflow/multicat_ethanol_1sccm_24hrstarved', 
                            'ethanol_starved_dusk': '/media/Orchard3/CONTROLED_5massflow/multicat_ethanol_1sccm_24hrstarved',
                            'ethanol_starved_night': '/media/Orchard3/CONTROLED_5massflow/multicat_ethanol_1sccm_24hrstarved', 
                            'ethanol_starved_morning': '/media/Orchard3/CONTROLED_5massflow/multicat_ethanol_1sccm_24hrstarved',
                            
                            'low_flow_warm_dusk': '/media/Orchard3/confounded_data/multicat_co2_ssr_delayed_lowflow_1sccmco2_hot',
                            'high_flow_co2_dusk': '/media/caveman/Orchard/miniorchard/co2_5percent/holyminiolfacto_co2_5percent_HCS_noonstarved_noonstarted',
                            'high_flow_warm_dusk': '/media/Orchard2/hot_100sccm',
                            'high_flow_arista_dusk': '/media/Orchard2/aristae',
                            
                            'high_flow_ethanol_dusk': '/media/Orchard2/eth_concentrations_100sccm_1_5_15',
                            }

# for figure 5
walking_arena_hcs_concentration = '/media/Orchard4/SYMMETRIC_HCS_CONCENTRATION'
walking_arena_symmetric_stim_genetics = '/media/Orchard4/SYMMETRIC_stim'

# for figure 3
co2_calibration_windtunnel = '/home/caveman/Dropbox (flyranch)/2017_CO2/analysis/figure3_windtunnel_walking/co2_concentration_wind_tunnel'

# for supplemental figure
shake_fly_mosquito_data = '/home/caveman/Dropbox (flyranch)/2017_CO2/analysis/sup_fig_shake_fly_mosquito/shake_fly_mosquito_data'


