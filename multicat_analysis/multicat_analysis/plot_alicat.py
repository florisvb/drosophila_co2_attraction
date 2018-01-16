import os, sys
import imp
import numpy as np
import matplotlib.pyplot as plt
import pandas

def get_filename(path, contains, does_not_contain=[]):
    cmd = 'ls ' + path
    ls = os.popen(cmd).read()
    all_filelist = ls.split('\n')
    try:
        all_filelist.remove('')
    except:
        pass
    filelist = []
    for i, filename in enumerate(all_filelist):
        if contains in filename:
            fileok = True
            for nc in does_not_contain:
                if nc in filename:
                    fileok = False
            if fileok:
                return os.path.join(path, filename)



def plot_alicat(path, unit='C'):
    configuration_filename = get_filename(path, 'config_', does_not_contain=['~', '.pyc'])
    configuration = imp.load_source('configuration', configuration_filename)
    config = configuration.Config(path)
    
    fig = plt.figure()
    
    if 'N1' in config.identifiercode:
        units = ['A', 'B', 'I', 'J', 'K', 'L']
        #ssr_port_0 = config.alicat_data['phidgets_interface_ssr']['output_data']['states_0']
        #ssr_port_1 = config.alicat_data['phidgets_interface_ssr']['output_data']['states_1']
        #ssr_t = config.alicat_data['phidgets_interface_ssr']['output_data']['t'] 
    elif 'N2' in config.identifiercode:
        units = ['D', 'E', 'F', 'G', 'H', 'C']
        #ssr_port_0 = config.alicat_data['phidgets_interface_ssr']['output_data']['states_2']
        #ssr_port_1 = config.alicat_data['phidgets_interface_ssr']['output_data']['states_3']
        #ssr_t = config.alicat_data['phidgets_interface_ssr']['output_data']['t'] 
        
    data = []
    for i, unit in enumerate(units):
        ax = fig.add_subplot(len(units),1,i+1)
        a = config.alicat_data['alicat_flow_rate_'+unit]['data']
        t = config.alicat_data['alicat_flow_rate_'+unit]['t']
        ax.plot(t, a)
        data.append(a)

        if np.max(a) < 20:
            ax.set_ylim(0,4)

        ax.set_ylabel(unit)
    