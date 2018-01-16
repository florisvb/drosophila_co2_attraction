import os

# before running analysis, set the environment variable co2_paper_locations to something that corresponds to this document
# e.g. export co2_paper_locations=analysiscavetech-original

try:
    co2_paper_locations_profile = os.environ["co2_paper_locations"]
except:
    co2_paper_locations_profile = "default"

if co2_paper_locations_profile == "analysiscavetech-original" or co2_paper_locations_profile == "default":
		basepath = "/home/caveman/Dropbox_flyranch/2017_CO2/analysis"
elif co2_paper_locations_profile == "analysiscavetech-organized":
		basepath = "/home/caveman/Dropbox_flyranch/2017_CO2/analysis"
elif co2_paper_locations_profile == "analysiscavetech-portable":
		basepath = "/home/caveman/Dropbox_flyranch/2017_CO2/analysis"
elif co2_paper_locations_profile == "uw":
        basepath = "/home/caveman/Dropbox_flyranch/2017_CO2/analysis"            

figure1_ferment_trapassay_flybottles = os.path.join(basepath, "figure1_ferment_trap_assay_flybottles/figure1_ferment_trapassay_flybottles.svg")
sup_figure_co2_calibration = os.path.join(basepath, "figure1_ferment_trap_assay_flybottles/sup_figure_co2_calibration.svg")

figure2_freeflight = os.path.join(basepath, "figure2_freeflight/figure2_freeflight.svg")

figure3_windtunnel_walking = os.path.join(basepath, "figure3_windtunnel_walking/figure3_windtunnelwalking_multiconcentration.svg")
figure3_windtunnel_walking_simple = os.path.join(basepath, "figure3_windtunnel_walking/figure3_windtunnelwalking_simple.svg")
figure3_windtunnel_walking_trajecs = os.path.join(basepath, "figure3_windtunnel_walking/figure3_windtunnelwalking_multiconcentration_trajecs.svg")
figure3_windtunnel_walking_trajecs = os.path.join(basepath, "figure3_windtunnel_walking/figure3_windtunnelwalking_multiconcentration_trajecs.svg")
figure3_windtunnel_walking_supp = os.path.join(basepath, "figure3_windtunnel_walking/figure3_windtunnelwalking_simple_supp.svg")
figure3_co2_calibration = os.path.join(basepath, "figure3_windtunnel_walking/figure_sup_co2_calibration.svg")

figure4_walking_arena_activity = os.path.join(basepath, "figure4_walking_arena_activity/figure4_walking_arena_activity.svg")

figure5_walking_arena_genetics_only = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_genetics_only.svg")
figure5_walking_arena_hcs_full = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_hcs_full.svg")
figure5_walking_arena_genetics_speed = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_genetics_speed.svg")
figure5_walking_arena_genetics_speed_control = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_genetics_speed_control.svg")

figure5_walking_arena_genetics_control = os.path.join(basepath, "figure5_walking_arena_genetics/figure5_speed_sorted_walking_genetics_control.svg")

