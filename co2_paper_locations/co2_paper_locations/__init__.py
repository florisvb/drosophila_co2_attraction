# before running analysis, set the environment variable co2_paper_locations to something that corresponds to this document
# e.g. export co2_paper_locations=analysiscavetech-original

#import data_locations
import figure_template_locations

import os
co2_paper_locations_profile = os.environ['co2_paper_locations']
print 'USING PROFILE: ', co2_paper_locations_profile

if co2_paper_locations_profile == 'analysiscavetech-original' or co2_paper_locations_profile == 'default':
	import data_locations_analysiscavetech_original as data_locations
elif co2_paper_locations_profile == 'analysiscavetech-portable':
	import data_locations_analysiscavetech_portable as data_locations
elif co2_paper_locations_profile == 'analysiscavetech-organized':
	import data_locations_analysiscavetech_organized as data_locations
elif co2_paper_locations_profile == 'uw':
    import data_locations_analysiscavetech_uw as data_locations