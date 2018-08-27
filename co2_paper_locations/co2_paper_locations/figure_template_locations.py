import os

# before running analysis, set the environment variable co2_paper_locations to something that corresponds to this document
# e.g. export co2_paper_locations=analysiscavetech-original

try:
    co2_paper_locations_profile = os.environ["co2_paper_locations"]
except:
    co2_paper_locations_profile = "default"

if co2_paper_locations_profile == "analysiscavetech-original" or co2_paper_locations_profile == "default":
		basepath = "/home/caveman/Dropbox_flyranch/2017_CO2/drosophila_co2_attraction"
elif co2_paper_locations_profile == "analysiscavetech-organized":
		basepath = "/home/caveman/Dropbox_flyranch/2017_CO2/drosophila_co2_attraction"
elif co2_paper_locations_profile == "analysiscavetech-portable":
		basepath = "/home/caveman/Dropbox_flyranch/2017_CO2/drosophila_co2_attraction"
elif co2_paper_locations_profile == "uw":
        basepath = "/home/caveman/Dropbox_flyranch/2017_CO2/drosophila_co2_attraction"            

figure1_ferment_trapassay_flybottles = os.path.join(basepath, "supp_fig_1/supp_fig_1_trapassay.svg")
sup_figure_co2_calibration = os.path.join(basepath, "supp_fig_1/supp_fig_1_trapassay.svg")

figure2_freeflight = os.path.join(basepath, "fig_1/fig_1_freeflight.svg")

figure3_windtunnel_walking = os.path.join(basepath, "fig_2/fig_2_walking.svg")
figure3_windtunnel_walking_simple = os.path.join(basepath, "fig_2/fig_2_walking.svg")
figure3_windtunnel_walking_trajecs = os.path.join(basepath, "fig_2/fig_2_walking.svg")
figure3_windtunnel_walking_trajecs = os.path.join(basepath, "fig_2/fig_2_walking.svg")

figure3_windtunnel_walking_supp = os.path.join(basepath, "supp_fig_2/supp_fig_2_walking.svg")
figure3_co2_calibration = os.path.join(basepath, "supp_fig_2/supp_fig_2_walking.svg")

figure4_walking_arena_activity = os.path.join(basepath, "fig_3/fig_3_activity.svg")

#figure5_walking_arena_genetics_only = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_genetics_only.svg")
#figure5_walking_arena_hcs_full = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_hcs_full.svg")
#figure5_walking_arena_genetics_speed = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_genetics_speed.svg")
#figure5_walking_arena_genetics_speed_control = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_genetics_speed_control.svg")

figure5_walking_arena_genetics_control = os.path.join(basepath, "supp_fig_6/supp_fig_6_co2control.svg")

#figure5_full_template = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_full_supplement_template.svg")

figure5_walking_arena_genetics_only_CO2 = os.path.join(basepath, "fig_4/fig_4_mutants.svg")
figure5_walking_arena_genetics_only_CO2_control = os.path.join(basepath, "supp_fig_6/supp_fig_6_co2control.svg")

figure5_single_concentration_model = os.path.join(basepath, "supp_fig_4/supp_fig_4.svg")

supp_fig_8_eth_vinegar = os.path.join(basepath, "supp_fig_8/supp_fig_8_eth_vinegar.svg")
supp_fig_7_ir25a_ir40a = os.path.join(basepath, "supp_fig_7/supp_fig_7_ir25a_ir40a.svg")