import matplotlib.pyplot as plt
import random
from matplotlib import colors
import numpy as np
import pandas
import itertools
import warnings
import numpy as np
import itertools
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
import scipy.signal as signal
import os 
from collections import defaultdict
import math
#%%Functions
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
#%%

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
    #%

def crosscorrelate(sua1, sua2, lag=None, n_pred=1, predictor=None,
                   display=False, kwargs={}):
    #Runs on samples, not spike times!
    assert predictor is 'shuffle' or predictor is None, "predictor must be \
    either None or 'shuffle'. Other predictors are not yet implemented."
    #Check whether sua1 and sua2 are SpikeTrains or arrays
    sua = []
    for x in (sua1, sua2):
        #if isinstance(x, SpikeTrain):
        if hasattr(x, 'spike_times'):
            sua.append(x.spike_times)
        elif x.ndim == 1:
            sua.append(x)
        elif x.ndim == 2 and (x.shape[0] == 1 or x.shape[1] == 1):
            sua.append(x.ravel())
        else:
            raise TypeError("sua1 and sua2 must be either instances of the" \
                            "SpikeTrain class or column/row vectors")
    sua1 = sua[0]
    sua2 = sua[1]
    if sua1.size < sua2.size:
        if lag is None:
            lag = np.ceil(10*np.mean(np.diff(sua1)))
        reverse = False
    else:
        if lag is None:
            lag = np.ceil(20*np.mean(np.diff(sua2)))
        sua1, sua2 = sua2, sua1
        reverse = True
    #construct predictor
    if predictor is 'shuffle':
        isi = np.diff(sua2)
        sua2_ = np.array([])
        for ni in range(1, n_pred+1):
            idx = np.random.permutation(isi.size-1)
            sua2_ = np.append(sua2_, np.add(np.insert(
                (np.cumsum(isi[idx])), 0, 0), sua2.min() + (
                np.random.exponential(isi.mean()))))
    #calculate cross differences in spike times
    differences = np.array([])
    pred = np.array([])
    for k in np.arange(0, sua1.size): #changed xrange() for np.arange()
        differences = np.append(differences, sua1[k] - sua2[np.nonzero(
            (sua2 > sua1[k] - lag) & (sua2 < sua1[k] + lag))])
    if predictor == 'shuffle':
        for k in np.arange(0, sua1.size): #changed xrange() for np.arange()
            pred = np.append(pred, sua1[k] - sua2_[np.nonzero(
                (sua2_ > sua1[k] - lag) & (sua2_ < sua1[k] + lag))])
    if reverse is True:
        differences = -differences
        pred = -pred
    norm = np.sqrt(sua1.size * sua2.size)
    return differences, pred, norm
#%Spike detection
def spike_detection_npx(npx_voltage, Thr):

    #Len_Thr = 0.5 * Thr
#
    b = npx_voltage - Thr
#
    ThrCrosses = np.diff(np.sign(b))
    #Len_ThrCrosses = np.diff(np.sign(b+Len_Thr))

    ThrCrossStart = []
    ThrCrossEnd = []
    ThrCrossStart.append(np.where(ThrCrosses <= -2)[0])
    ThrCrossEnd.append(np.where(ThrCrosses >= 2)[0])
#
    spike_times = []
    spike_end = np.empty((len(ThrCrossStart[0]),1), dtype = int)
    for spike_start in range(len(ThrCrossStart[0])):
        spike_end[spike_start] = ThrCrossStart[0][spike_start] + 10
    
    for i in range(len(ThrCrossStart[0])):
        if ThrCrossStart[0][i] + 10< len(npx_voltage):
            spike_times.append( (ThrCrossStart[0][i] + np.argmin(b[ThrCrossStart[0][i]:spike_end[i][0]] )) )
        else:
            pass
    #Len_ThrCrossStart = []
    #Len_ThrCrossEnd = []
    #Len_ThrCrossStart.append(np.where(Len_ThrCrosses == -2))
    #Len_ThrCrossEnd.append(np.where(Len_ThrCrosses == 2))
    #plt.scatter(spike_times, b[spike_times]+thr, c = 'r')
    #plt.plot(b+thr)
    #plt.hlines(thr,0,30000)
#
    #ThrCrossDurs = [ThrCrossEnd[0][0][i] - ThrCrossStart[0][0][i] for i in range(len(ThrCrossEnd[0][0]))]
    #Len_ThrCrossDurs = [Len_ThrCrossEnd[0][0][i] - Len_ThrCrossStart[0][0][i] for i in range(len(Len_ThrCrossStart[0][0]))]

    #spike_times = []

    #for i in range(len(ThrCrossEnd[0][0])):
     #   spike_times.append( (ThrCrossStart[0][0][i] + np.argmin(b[ThrCrossStart[0][0][i]:ThrCrossEnd[0][0][i]] )) )

    spike_times = np.asarray(spike_times)
#
    return spike_times

#%
def crosscorrel_on_probe(histogram_array, rms_array, cell_id, probe_map, width, height, rms_min, rms_max, central_chan, num_rows_side = 7):
    import matplotlib.gridspec as gridspec
    from matplotlib import colors
    import matplotlib.cm as cmx
#test loads#%%#%%
   #histogram_array = crosscorrels_30chans
   #rms_array = rms[cells_to_analyse.index(cell), :]
   #cell_id = cells_to_analyse.index(cell)
   #probe_map = npx_mapcsv
   #width = 5.4
   #height = 9.31
   #rms_min = np.min(rms[cells_to_analyse.index(cell),:])
   #rms_max = np.max(rms[cells_to_analyse.index(cell),:])
   #central_chan = npx_chan_peaks[0,0]
   #num_rows_side = 7 
   #
#
    rc_channel = np.where(probe_map == central_chan)
    sub_map = probe_map[rc_channel[0][0]-num_rows_side:rc_channel[0][0]+num_rows_side+1,:]
    electrodes = sub_map[sub_map > 0 ]

    cm= plt.get_cmap('cool')

    
    cNorm=colors.Normalize(vmin=rms_min, vmax= rms_max )#np.max(electrodes))
    scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)
    #
    fig, axs= plt.subplots(15,4,figsize=cm2inch(width,height), sharex=True, sharey=True)
    #

    for chan in range(len(electrodes)):
        colorVal=scalarMap.to_rgba( rms_array[chan] )
        r = np.where(sub_map == electrodes[chan])[0][0]
        c = np.where(sub_map == electrodes[chan])[1][0]#replaced bins 101 by 21
        cc = axs[r, c].hist(crosscorrels_30chans[chan][0]/30.0, bins = 101,  align = 'mid', normed = False, color = 'k')
        #plt.subplot2grid((15,4), (rc_channel[0][0], rc_channel[1][0]))
        plt.subplots_adjust(top = 0.98, bottom = 0.02, left = 0.1, right = 0.9, wspace = 0.1, hspace = 0.0)
        axs[r, c].axis('off')
        plt.xlim(-50,50)#-50 50
        if electrodes[chan] == central_chan:
            axs[r, c].set_title('%s'%(electrodes[chan]), fontsize = 5, y = 0.35, alpha = 1, loc='right', fontweight = 'bold')
        else:
            axs[r, c].set_title('%s'%(electrodes[chan]), fontsize = 5, y = 0.35, alpha = 0.6, loc='right')
    scalebar_loc = np.where(sub_map[:,0]==-1)[0]
    axs[scalebar_loc[0],0].vlines(-50,0, int(axs[0,0].get_ylim()[1]), lw = 2)
    axs[scalebar_loc[0],0].text(-45, 0,'%s'%(int(axs[0,0].get_ylim()[1])), fontsize = 8)
    
    scalebar_loc = np.where(sub_map[:,0]==-1)[0]
    axs[scalebar_loc[-1], 0].hlines(int(axs[0,0].get_ylim()[1])/10, -50, 0, lw = 2)
    axs[scalebar_loc[-1], 0]. text (-45, int(axs[0,0].get_ylim()[1])/5, '50 ms', fontsize = '8')

    #colorbar_loc = np.where(sub_map[:,0]==-1)[0]
    #axs[colorbar_loc[-2], 0]. text (-45, int(axs[0,0].get_ylim()[1])/5, '[%s %s uV]' %( round(rms_min,1), round(rms_max,1) ), fontsize = '6')
    
    for row in range(15):
        for col in range(4):
            axs[row, col].axis('off')
            

    

#cbar.ax.set_xticklabels(['0', '%s' %(100*round(histmax/2.0,2))+'%', '%s' %(100* round(histmax,2) )+'%'])
    
    #

#%
def return_nearest_power10(number):
    if number/1000.0 > 1:
        return int(np.ceil(number/1000.0)*1000)
    elif number/100.0 > 1:
        return int(np.ceil(number/100.0)*100)
    elif number/10.0 > 1:
        return int(np.ceil(number/10.0)*10)
    else:
        return number

#%%Data frames and structures
listcells = os.listdir('E:/code/for analysis')
suffix = '_aligned_analysis'

data_summary = pandas.read_excel('C:/Users/Andre Marques-Smith/Dropbox/Data Summary.xlsx')

cells_above_10uV = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 10].tolist()
#excluded_cells = [2, 6, 9, 10, 12, 21, 23]
excluded_cells = [2, 5, 9, 12, 21, 23]
cells_to_analyse = [ cells_above_10uV[i] for i in range(len(cells_above_10uV)) if i not in excluded_cells]

zoomed_histo = [[]] * len(cells_to_analyse)
row_zoom = 0
rms = np.empty((len(cells_to_analyse),30), dtype = float)
Thr = np.empty((len(cells_to_analyse),30), dtype = float)

#%%Batch Analysis for every cell out of 21 further analysed
#output_dir = 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/repro/fig5/'

for cell in cells_to_analyse[20:21]:
    paths = defaultdict(list)
    cell_idx = listcells.index(data_summary.loc[cell]['Cell']+suffix)

    aligned_directory = 'E:/code/for analysis/'+listcells[cell_idx]
    
    for file in os.listdir(aligned_directory):
        print file
        if file.endswith('meta.npy'):
            paths['meta'].append(aligned_directory + '/' + file)
        elif file.endswith('patch_preprocessed.npy'):
            paths['patch_v'].append(aligned_directory + '/' + file)
        elif file.endswith('sta_wc.npy'):
            paths['patch_sta'].append(aligned_directory + '/' + file)
        elif file.endswith('wc_spike_samples.npy'):
            paths['patch_spikes'].append(aligned_directory + '/' + file)
        elif file.endswith('npx_preprocessed.bin'):
            paths['npx_v'].append(aligned_directory + '/' + file)
        elif file.endswith('sta_windows_chan_array.npy'):
            paths['npx_sta_array'].append(aligned_directory + '/' + file)
        elif file.endswith('sta_np_by_channel.npy'):
            paths['npx_sta_mean'].append(aligned_directory + '/' + file)
        elif file.endswith('channels_by_pk.npy'):
            paths['npx_chan_peaks'].append(aligned_directory + '/' + file)
        else:
            pass

    expt_meta = np.load(paths['meta'][0]).item()
    
    patch_v = np.load(paths['patch_v'][0])#, dtype = 'float64')
    patch_sta = np.load(paths['patch_sta'][0])#, dtype = 'float64')
    patch_spikes = np.load(paths['patch_spikes'][0])#, dtype = 'int16')
    #
    npx_voltage = np.memmap( paths['npx_v'][0], mode = 'c', dtype=np.int16 )
    npx_voltage = npx_voltage.reshape((expt_meta['npx'][0][0], expt_meta['npx'][0][1]))
    #
    npx_sta_array = np.load(paths['npx_sta_array'][0])#, dtype = 'int16')
    npx_sta_mean = np.load(paths['npx_sta_mean'][0])#, dtype = 'float64')
    #npx_chan_peaks = np.load(paths['npx_chan_peaks'][0])#, dtype = 'float64')
    npx_chan_peaks = [np.max(npx_sta_mean[i]) - np.min(npx_sta_mean[i]) for i in range(384)]
    
    central_chan = np.where(npx_chan_peaks ==np.max(npx_chan_peaks))[0][0]
    #
    m = len(patch_v)/float(len(npx_voltage[0,:]))
    patch_spikes_npx = np.asarray([int(patch_spikes[i] / m) for i in range(len(patch_spikes))])
    nc = cells_to_analyse.index(cell)
#%%In Development
dot_prod = np.empty((len(patch_spikes),3), dtype = 'float')

non_spikes = []

rand_times = random.sample(range(30, len(npx_voltage[0])-31), len(patch_spikes))

chan = central_chan-10

while len(non_spikes)< len(patch_spikes):
    r = random.randint(30, len(npx_voltage[0])-31)
    if r not in patch_spikes_npx:
        if r not in non_spikes:
            non_spikes.append(r)

for spike in range(len(patch_spikes)):
    dot_prod[spike,0] = np.dot(unit_vector(npx_sta_array[chan,30:90,spike]), unit_vector(npx_sta_mean[chan, 30:90]))
    dot_prod[spike,1] = np.dot(unit_vector(npx_voltage[chan,non_spikes[spike]-30:non_spikes[spike]+30]), unit_vector(npx_sta_mean[chan, 30:90]))
    dot_prod[spike,2] = np.dot(unit_vector(npx_voltage[chan,rand_times[spike]-30:rand_times[spike]+30]), unit_vector(npx_sta_mean[chan, 30:90]))
#
weights = np.ones_like(dot_prod[:,0])/float(len(dot_prod[:,0]))
plt.hist(dot_prod[:,0], bins = 100, alpha = 0.5, color = 'b', weights = weights)    
plt.hist(dot_prod[:,1], bins = 100, alpha = 0.5, color = 'r', weights = weights)

#plt.hist(dot_prod[:,2], bins = np.arange(0,1.05,0.01), alpha = 0.1, color = 'm')
    
    
    
    #%%
    #central_chan = int(npx_chan_peaks[0,0])
    num_rows_side = 7
    #
    rc_channel = np.where(npx_mapcsv == central_chan)
    sub_map = npx_mapcsv[rc_channel[0][0]-num_rows_side:rc_channel[0][0]+num_rows_side+1,:]
    electrodes = sub_map[sub_map > 0 ]
    
#Detect spikes I - calculate threshold as 7* the median absolute deviation for each ell and each of the 30 channels closest to that cell
    for channel in range(len(electrodes)):
        Thr[cells_to_analyse.index(cell), channel] = -7*np.median(np.abs(npx_voltage[electrodes[channel],:] - np.median(npx_voltage[electrodes[channel],:] )))#-5*np.median(np.absolute(npx_voltage[channel])/0.6745)

    spikes_30chans = [[]]*30
#Detect spikes II - find spikes as crossings of the threshold Thr for each of the 30 channels in every cell
    for channel in range(len(electrodes)):
        spikes_30chans[channel] = spike_detection_npx(npx_voltage[electrodes[channel]], Thr[cells_to_analyse.index(cell), channel])
#Run cross-correlation function and find lag between patch spike and every spike found for 30 channels in every cell that was +- 50 ms of the patch spike peak
    crosscorrels_30chans = [[]]*30
    
    for channel in range(len(electrodes)):
        crosscorrels_30chans[channel] = crosscorrelate(spikes_30chans[channel], patch_spikes_npx, lag = 1500)
        rms[cells_to_analyse.index(cell), channel] = np.sqrt(np.nanmean(npx_voltage[channel,:]**2))
    #
    #
    #yscale = crosscorrel_on_probe(crosscorrels_30chans ,npx_mapcsv, width = 5.4, height = 9.31, color_probe = 'b', central_chan = npx_chan_peaks[0,0], num_rows_side = 7)
    crosscorrel_on_probe(crosscorrels_30chans, rms[cells_to_analyse.index(cell), :], cells_to_analyse.index(cell), npx_mapcsv, width = 5.4, height = 9.31, rms_min = np.min(rms[cells_to_analyse.index(cell),:]), rms_max = np.max(rms[cells_to_analyse.index(cell),:]), central_chan = npx_chan_peaks[0,0], num_rows_side = 7)
#
    plt.savefig( output_dir + 'Fig 5_'+data_summary.loc[cell]['Cell']+'_correl_array_mad7.png', dpi = 1200)
    plt.close()

    npx_voltage = None
    
    
    zoomed_histo[cells_to_analyse.index(cell)] = crosscorrels_30chans[15][0]
    
 
    
#%%Figure 5D - Plot the peri-event time histogram of extracellular vs patch spikes fpr the closest channel to the soma for each of the 21 cells

    fig, axs= plt.subplots(3,7,figsize=cm2inch(17.59,5.72), sharex='col', sharey=False)
    axs = axs.ravel()
    #

    for cell in range(len(cells_to_analyse)):
        #r = np.where(sub_map == electrodes[chan])[0][0]
        #c = np.where(sub_map == electrodes[chan])[1][0] bins replaced by 100, removed range as not needed
        axs[cell].hist(zoomed_histo[cell]/30.0, bins = 101, align = 'mid', normed = False, color = 'k')
        ymax = return_nearest_power10(np.max(np.histogram(zoomed_histo[cell]/30.0, bins = 101, range=(-50,50))[0]))
        axs[cell].set_yticks([0, ymax])
        axs[cell].set_title(data_summary.loc[cells_to_analyse[cell]]['Cell'], fontsize = 8, y = 0.9)
        axs[cell].tick_params(labelsize=6)
        plt.subplots_adjust(hspace = 0.7, wspace = 0.7, bottom = 0.1, top = 0.91, left = 0.07, right = 0.97)

        plt.xlim(-50,50)
        
    #%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 4_materials/' + 'Fig 4_zoom_all_cells.png', dpi = 1200)
#%%Figure 5E - Find which time bin (relative to patch spike) is the time bin where ground-truth extracellular spikes are occurring and compute what proportion that is,
#compared to the total number of patch spike times. IE, does it detect more spikes than the cell fired or less?
total_patch_spikes = np.empty((len(cells_to_analyse),1), dtype = float)
total_npx_spikes = np.empty((len(cells_to_analyse),1), dtype = float)
negative_pks = np.empty((len(cells_to_analyse),1), dtype = float)
proportion = np.empty((len(cells_to_analyse),1), dtype = float)
    
for cell in cells_to_analyse:
    aligned_directory = 'E:/code/for analysis/'+data_summary.loc[cell]['Cell']+suffix
    
    for file in os.listdir(aligned_directory):
        if file.endswith('sta_np_by_channel.npy'):
            npx_sta = np.load(aligned_directory + '/' + file)
            negative_pks[cells_to_analyse.index(cell)] = np.min(npx_sta)
        elif file.endswith('wc_spike_samples.npy'):
            patch_spikes = np.load(aligned_directory+ '/' + file)
            total_patch_spikes[cells_to_analyse.index(cell)] = len(patch_spikes)
        else:            
            pass
    which_bin = 48+np.argmax(np.histogram(zoomed_histo[ cells_to_analyse.index(cell) ]/30.0, bins = 101, range=(-50,50))[0][48:51] )
    total_npx_spikes[cells_to_analyse.index(cell)] = np.histogram(zoomed_histo[ cells_to_analyse.index(cell) ]/30.0, bins = 101, range=(-50,50))[0][which_bin]

#
for cell in range(21):
    proportion[cell] = total_npx_spikes[cell]/total_patch_spikes[cell]
#Detectability plot, Figure 5E

fig = plt.figure(figsize = (5,5))
for cell in range(21):
    detectability = plt.scatter( np.log2(negative_pks[cell]/np.median(Thr[cell,:])), proportion[cell], cmap = 'viridis', c = -negative_pks[cell], vmin = 0, vmax = 200 )
for i, txt in enumerate(data_summary.loc[cells_to_analyse]['Cell']):
    plt.annotate(txt, (np.log2(negative_pks[i]/np.median(Thr[i,:])), proportion[i]), xytext=((-1)**i+3,10+i), textcoords='offset pixels', horizontalalignment='left' )
plt.xlabel('Negative peak/Spike threshold (log2)', fontweight = 'bold')
plt.ylabel('Extracellular divided by\npatch spike count', fontweight = 'bold')
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25])
plt.ylim(0.0,1.25)
plt.xlim(-2,4)
col = plt.colorbar(detectability, orientation = 'horizontal', shrink = 0.5, ticks = range(0,250,50), pad = 0.15)
col.ax.set_xlabel('Amplitude of Spike Negative Peak (uV)')
plt.vlines(0, 0.0, 1.5, linestyle = '--', lw = 0.5)
plt.hlines(0.5, -2, 4, linestyle = '--', lw = 0.5)
plt.subplots_adjust(top = 0.9, bottom = 0.1, left = 0.2, right = 0.8)
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 4_materials/' + 'Fig 4_detection.png', dpi = 1200)
#%%Re-plot Figure 5D using color code of 5E

    fig, axs= plt.subplots(3,7,figsize=cm2inch(17.59,5.72), sharex='col', sharey=False)
    axs = axs.ravel()
    #
    cm = plt.get_cmap('viridis')
    cNorm=colors.Normalize(vmin=0, vmax= 200 )#np.max(electrodes))
    scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)

    

    for cell in range(len(cells_to_analyse)):
        colorVal=scalarMap.to_rgba( -negative_pks[cell] )
        #r = np.where(sub_map == electrodes[chan])[0][0]
        #c = np.where(sub_map == electrodes[chan])[1][0] bins replaced by 100, removed range as not needed
        axs[cell].hist(zoomed_histo[cell]/30.0, bins = 101, align = 'mid', normed = False, color = colorVal)
        ymax = return_nearest_power10(np.max(np.histogram(zoomed_histo[cell]/30.0, bins = 101, range=(-50,50))[0]))
        axs[cell].set_yticks([0, ymax])
        axs[cell].set_title(data_summary.loc[cells_to_analyse[cell]]['Cell'], fontsize = 8, y = 0.9)
        axs[cell].tick_params(labelsize=6)
        plt.subplots_adjust(hspace = 0.7, wspace = 0.7, bottom = 0.1, top = 0.91, left = 0.07, right = 0.97)

        plt.xlim(-50,50)
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 4_materials/' + 'Fig 4_zoom_all_cells.png', dpi = 1200)