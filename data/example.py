# Instructions for regenerating figures

### Unzip the files

### Install the following:
* https://github.com/FlyRanch/figurefirst
* https://github.com/florisvb/FlyPlotLib
* https://github.com/florisvb/FlyStat

import figurefirst

LAYOUT = path_to_svg_template_file
OUTPUT = path_to_svg_output_file # could be the same as the LAYOUT
DATA = path_to_the_pickle_file_containing_the_data

figurefirst.regenerate.replot(LAYOUT, OUTPUT, DATA)

# Instructions for opening the DATA.dillpickle files

import dill as pickle
f = open(DATA, 'r')
data = pickle.load(f)
f.close()

'''
Data is a dictionary, where each key corresponds to a sub-panel in the associated figure template. 
The value corresponding to that key is a list of the plotting actions used to create the figure,
these actions include calls to matplotlib, https://github.com/florisvb/FlyPlotLib, and figurefirst.
Each action is a dictionary that contains data, plotting instructions, and a short description.
'''

### Example of accessing the data

f = open('fig_1_freeflight_data.dillpickle', 'r')
data = pickle.load(f)
f.close()
plot_action = data[(u'windtunnel_quanitified', u'landing_pad')][0]

print plot_action['title']

for i, arg in enumerate(plot_action['args']):
    print plot_action['args_description'][i]
    print plot_action['args'][i]

print plot_action['function'] 
print plot_action['kwargs']





