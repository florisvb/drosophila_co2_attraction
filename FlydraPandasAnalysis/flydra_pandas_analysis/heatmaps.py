import multi_tracker_analysis as mta
import numpy as np
import matplotlib.pyplot as plt
import flydra_pandas_dataset
import os
from optparse import OptionParser
from matplotlib.backends.backend_pdf import PdfPages
import imp
import flydra_pandas_dataset as fpd

def get_heatmap_difference_for_ops(pd, config, iterations=10, resolution=0.01, vmin=-10, vmax=10, figsize=(1.9*4,.5*4), pval=0.05):
    h_ons = []
    h_offs = []
    for f in range(0, np.max(pd.dataset_code)+1):
        pd_f = pd[pd.dataset_code==f]
        
        for op in range(0,8):
            pd_odor_on = fpd.get_pd_for_odor_on(pd_f, config, odor_presentations=[op])
            pd_odor_off = fpd.get_pd_for_odor_off(pd_f, config, 30, 60, odor_presentations=[op])
            
            h_on, h_off = get_heatmap_difference(pd_odor_on, pd_odor_off, resolution=resolution)
            h_ons.append(h_on)
            h_offs.append(h_off)
            
    actual_difference = np.mean(h_ons, axis=0) - np.mean(h_offs, axis=0)
    
    all_hs = []
    all_hs.extend(h_ons)
    all_hs.extend(h_offs)
    
    differences = []
    for i in range(iterations):
        indices_on = np.random.randint(0, len(all_hs), len(h_ons))
        indices_off = np.random.randint(0, len(all_hs), len(h_offs))

        fake_on = np.mean([all_hs[n] for n in indices_on], axis=0)
        fake_off = np.mean([all_hs[n] for n in indices_off], axis=0)
        
        differences.append(fake_on-fake_off)
    differences = np.array(differences)
    
    pvals = np.ones_like(actual_difference)
    for r in range(0, differences.shape[1]):
        for c in range(0, differences.shape[2]):
            q = differences[:,r,c]
            q.sort()
            iq = np.argmin( np.abs(q-actual_difference[r,c]) )
            p = 1 - np.abs((iq - iterations/2.) / (iterations/2.))
            pvals[r,c] = p
            
    indices = np.where(pvals > pval)
    actual_difference[indices] = 0
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0,0,1,1])
    binsx = np.arange(-.5,1.25,resolution)
    binsy = np.arange(-.2,.2,resolution)
    
    ax.imshow(actual_difference.T, cmap=plt.get_cmap('bwr'), vmin=vmin, vmax=vmax, origin='lower', interpolation='nearest', extent=[binsx[0],binsx[-1],binsy[0],binsy[-1]])
    
    ax.set_ylim(-.2,.2)
    ax.set_xlim(-.5,1.25)
    ax.set_xticks([])
    ax.set_yticks([])
    
    save = os.path.join(config.path, 'results/heatmapdiff.pdf')
    fig.savefig(save, format='pdf')
    
    return actual_difference, differences, pvals
        
    
    
    
    
def get_heatmap_difference(pd_odor_on, pd_odor_off, resolution=0.01):
    binsx = np.arange(-.5,1.25,resolution)
    binsy = np.arange(-.2,.2,resolution)
    h_on = mta.plot.get_heatmap(pd_odor_on, binsx, binsy, position_y='position_z', position_z='position_y', position_z_slice=[-.1, .1])
    h_off = mta.plot.get_heatmap(pd_odor_off, binsx, binsy, position_y='position_z', position_z='position_y', position_z_slice=[-.1, .1])
    return h_on, h_off
    
def plot_heatmap_difference(pd_odor_on, pd_odor_off, vmin=-500, vmax=500, resolution=0.01, save='', figsize=(1.9*4,.5*4)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0,0,1,1])
    binsx = np.arange(-.5,1.25,resolution)
    binsy = np.arange(-.2,.2,resolution)
    h_on = mta.plot.get_heatmap(pd_odor_on, binsx, binsy, position_y='position_z', position_z='position_y', position_z_slice=[-.1, .1])
    h_off = mta.plot.get_heatmap(pd_odor_off, binsx, binsy, position_y='position_z', position_z='position_y', position_z_slice=[-.1, .1])
    h_diff = h_on - h_off
    
    if vmax is None:
        vmax = pd_odor_on.shape[0]*1e-5
        print 'vmax: ', vmax
        
    ax.imshow(h_diff.T, cmap=plt.get_cmap('bwr'), vmin=vmin, vmax=vmax, origin='lower', interpolation='nearest', extent=[binsx[0],binsx[-1],binsy[0],binsy[-1]])
    
    ax.set_ylim(-.2,.2)
    ax.set_xlim(-.5,1.25)
    ax.set_xticks([])
    ax.set_yticks([])
    if len(save) > 0:
        ax.figure.savefig(save, format='pdf')
    return fig
    
def plot_xz_heatmap(pd, vmax=None, resolution=0.003, save='', figsize=(1.9*4,.5*4)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0,0,1,1])
    binsx = np.arange(-.5,1.25,resolution)
    binsy = np.arange(-.2,.2,resolution)
    if vmax is None:
        vmax = pd.shape[0]*1e-5
        print 'vmax: ', vmax
    ax = mta.plot.plot_heatmap(pd, binsx, binsy, ax=ax, vmax=vmax, vmin=0, position_x='position_x', position_y='position_z', position_z='position_y', position_z_slice=[-.1,.1])
    ax.set_ylim(-.2,.2)
    ax.set_xlim(-.5,1.25)
    ax.set_xticks([])
    ax.set_yticks([])
    if len(save) > 0:
        ax.figure.savefig(save, format='pdf')
    return fig
    
def plot_xy_heatmap(pd, vmax=None, resolution=0.001, save='', zslice=[-.01, .02]):
    fig = plt.figure(figsize=(1.9*4,.5*4))
    ax = fig.add_axes([0,0,1,1])
    binsx = np.arange(-.5,1.25,resolution)
    binsy = np.arange(-.2,.2,resolution)
    if vmax is None:
        vmax = pd.shape[0]*3e-5
    ax = mta.plot.plot_heatmap(pd, binsx, binsy, ax=ax, vmax=vmax, vmin=0, position_x='position_x', position_y='position_y', position_z='position_z', position_z_slice=zslice)
    ax.set_ylim(-.2,.2)
    ax.set_xlim(-.5,1.25)
    ax.set_xticks([])
    ax.set_yticks([])
    if len(save) > 0:
        ax.figure.savefig(save, format='pdf')
    return fig
    
def plot_yz_heatmap(pd, vmax=None, resolution=0.001, save='', zslice=[-0.1, .4]):
    fig = plt.figure(figsize=(0.5*4,.5*4))
    ax = fig.add_axes([0,0,1,1])
    binsx = np.arange(-.2,.2,resolution)
    binsy = np.arange(-.2,.2,resolution)
    if vmax is None:
        vmax = pd.shape[0]*3e-5
    ax = mta.plot.plot_heatmap(pd, binsx, binsy, ax=ax, vmax=vmax, vmin=0, position_x='position_y', position_y='position_z', position_z='position_x', position_z_slice=zslice)
    ax.set_ylim(-.2,.2)
    ax.set_xlim(-.2,.2)
    ax.set_xticks([])
    ax.set_yticks([])
    if len(save) > 0:
        ax.figure.savefig(save, format='pdf')
    return fig
    
def plot_heatmaps_for_on_and_off(pd, config, savepath=''):
    pd_odor_on = flydra_pandas_dataset.get_pd_for_odor_on(pd, config)
    pd_odor_off = flydra_pandas_dataset.get_pd_for_odor_off(pd, config, start_minute=20, end_minute=50)
    
    save = os.path.join(savepath, 'xy_odor_on.pdf')
    plot_xy_heatmap(pd_odor_on, save=save)
    
    save = os.path.join(savepath, 'xz_odor_on.pdf')
    plot_xz_heatmap(pd_odor_on, save=save)
    
    save = os.path.join(savepath, 'xy_odor_off.pdf')
    plot_xy_heatmap(pd_odor_off, save=save)
    
    save = os.path.join(savepath, 'xz_odor_off.pdf')
    plot_xz_heatmap(pd_odor_off, save=save)
    
def plot_xz_heatmap_book(path, odor=True):
    s = 'results/heatmap_book' + str(odor) + '.pdf'
    figure_path = os.path.join(path, s)
    pp = PdfPages(figure_path)
    

    contains = ['trackedobjects.pickle']
    filenames = fpd.get_filename(path, contains)
    
    print filenames
    
    config_filename = fpd.get_filename(path, 'config.py')[0]  
    Config = imp.load_source('Config', config_filename)
    config = Config.Config(path)
    
    pds = []
    odor_controls = []
    for f, filename in enumerate(filenames):
        identifiercode = os.path.basename(filename).split('_trackedobjects')[0]
        print identifiercode
        print 'ODOR: ', odor
        pd, odor_control = fpd.load_and_align_timestamps_for_pickled_pandas_dataframe(filename, identifiercode)
        side = config.sides[identifiercode]
        if side == 'left':
            print identifiercode, 'left', 'flipped'
            pd.position_y = pd.position_y*-1
        pd['dataset_code'] = f
        config.odor_control = odor_control
            
        pd_odor_on = fpd.get_pd_for_odor_on(pd, config)
        pd_odor_off = fpd.get_pd_for_odor_off(pd, config)
    
        if odor:
            fig = plot_xz_heatmap(pd_odor_on, vmax=None, resolution=0.003, save='', figsize=(1.9*4,.5*4))
        else:
            fig = plot_xz_heatmap(pd_odor_off, vmax=None, resolution=0.003, save='', figsize=(1.9*4,.5*4))
        pp.savefig()
        plt.close('all')
        
    pp.close()    
    
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("--path", type="str", dest="path", default='',
                        help="path data")
    parser.add_option("--save", type="str", dest="save", default='',
                        help="path to save to")
    (options, args) = parser.parse_args()
    
    pd, config = flydra_pandas_dataset.load_pickled_pandas_dataframes_as_one_dataframe(options.path)
    #plot_heatmaps_for_on_and_off(pd, config, savepath=options.save)
    
    pd_odor_on = flydra_pandas_dataset.get_pd_for_odor_on(pd, config)
    pd_odor_off = flydra_pandas_dataset.get_pd_for_odor_off(pd, config, 30, 60)
    
    
    a, d, p = get_heatmap_difference_for_ops(pd, config, iterations=1000, resolution=0.01, vmax=10, vmin=-10, pval=0.05)
    
    
    plot_xz_heatmap_book(options.path, odor=False)
    plot_xz_heatmap_book(options.path, odor=True)
    
    
    
