# -*- coding: utf-8 -*-
#imports
import pandas
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict
import random
import matplotlib.gridspec as gridspec
from matplotlib import colors
import matplotlib.cm as cmx
from scipy.signal import butter, lfilter, freqz
import math
#%%
listcells = os.listdir('E:/code/for analysis')
suffix = '_aligned_analysis'

data_summary = pandas.read_excel('C:/Users/Andre Marques-Smith/Dropbox/Data Summary.xlsx')
cells_above_10uV = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 10].tolist()
#excluded_cells = [2, 6, 9, 10, 12, 21, 23]
excluded_cells = [2, 5, 9, 12, 21, 23]
#cells_to_analyse = [ cells_above_10uV[i] for i in range(len(cells_above_10uV)) if i not in excluded_cells]
cells_to_analyse = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 50].tolist()
n = len(cells_to_analyse)

DV_dist = np.zeros((n, 384,3), dtype = 'float')
#%%
for cell in cells_to_analyse:
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
    #%
    expt_meta = np.load(paths['meta'][0]).item()
    
    patch_v = np.load(paths['patch_v'][0])#, dtype = 'float64')
    patch_sta = np.load(paths['patch_sta'][0])#, dtype = 'float64')
    patch_spikes = np.load(paths['patch_spikes'][0])#, dtype = 'int16')
    #
    #npx_voltage = np.memmap( paths['npx_v'][0], mode = 'c', dtype=np.int16 )
    #npx_voltage = npx_voltage.reshape((expt_meta['npx'][0][0], expt_meta['npx'][0][1]))
    #
    npx_sta_array = np.load(paths['npx_sta_array'][0])#, dtype = 'int16')
    npx_sta_mean = np.load(paths['npx_sta_mean'][0])#, dtype = 'float64')
    #npx_chan_peaks = np.load(paths['npx_chan_peaks'][0])#, dtype = 'float64')
    
    npx_chan_peaks = [np.max(npx_sta_mean[i]) - np.min(npx_sta_mean[i]) for i in range(384)]
    
    central_chan = np.where(npx_chan_peaks ==np.max(npx_chan_peaks))[0][0]
#
    nc = cells_to_analyse.index(cell)
    x0,y0 = chanMapcsv[np.where(chanMapcsv[:,0] == central_chan),1:3][0][0]
    #z = data_summary.loc[data_summary['Cell'] == listcells[cell][:-17]]['Distance']
    neg_pk = np.argmin(npx_sta_mean[central_chan])
    central_pk2pk = np.max(npx_sta_mean[central_chan]) - np.min(npx_sta_mean[central_chan]) 
    
    for chan in range(384):
        DV_dist[nc, chan,1] = (np.max(npx_sta_mean[chan]) - np.min(npx_sta_mean[chan])) / central_pk2pk #relative pk2pk amplitude
        DV_dist[nc, chan,2] = 1000/30000.0 * (np.argmin(npx_sta_mean[chan]) - neg_pk) # latency between channel negative peak and somatic negative peak
        x1,y1 = chanMapcsv[np.where(chanMapcsv[:,0] == chan),1:3][0][0]
        DV_dist[nc, chan,0] = int(math.sqrt( (x0 - x1)**2 + (y0-y1)**2) ) # Distance between channel and somatic channel
        if y1 < y0:
            DV_dist[nc, chan,0] *= -1
    #%%
DV_dist_align = np.zeros((n, 200,3), dtype = 'float')
for c in range(n):
    centre = np.where(DV_dist[c,:,0] ==0)[0][0]
    for col in range(3):
        DV_dist_align[c,:,col] = DV_dist[c,centre-100:centre+100,col]
#%%
    #for cell in range(21):
    for c in range(n):
        plt.scatter(DV_dist_align[c,:,0], DV_dist_align[c,:,1], s = 3)
        plt.xlim(-1000,1000)
        #plt.vlines(0,0,1, linestyle = '--', lw = 0.7)
        print data_summary.loc[cells_to_analyse[c]]['Cell']
    plt.plot( np.linspace(-1000,1000,200), np.average(DV_dist_align[:,:,1], axis = 0))
    plt.hlines(0.12,-1000,1000,colors='r',linestyles='--')
#%%

#%%
sig = 0.12
for cell in range(n):
    #colorVal=scalarMap.to_rgba( n ) 
    sig_pks = np.where(DV_dist_align[cell,:,1] >= sig)[0]
    #for chan in sig_pks:
    plt.scatter(DV_dist_align[cell,sig_pks,2], DV_dist_align[cell,sig_pks,0],  s = 3, alpha = 0.5)
    #plt.xlim(-1000,1000)
#%%
for cell in range(21):
    sig_pks = np.where(DV_dist_align[cell,:,1] >= sig)[0]
    plt.plot(DV_dist_align[cell,:,0], DV_dist_align[cell,:,2])
    plt.xlim(-1000,1000)
#%%
d = 50
for chan in range(centre-d,centre+d):
        plt.plot( (npx_sta_mean[chan]/ -np.min(npx_sta_mean[chan])) + chan - (centre-d) , alpha = 0.3, lw = 0.7, c = 'b' )
        plt.scatter(np.argmin(npx_sta_mean[chan]), np.min(npx_sta_mean[chan]/-np.min(npx_sta_mean[chan])) + chan - (centre-d), c = 'r', s = 3)
        
#%%
pos_pks = np.where(DV_dist[:,2] >= 2)[0]
scale =  (np.max(npx_sta_mean[central_chan]) - np.min(npx_sta_mean[central_chan])) / (np.max(npx_sta_mean[pos_pks]) - np.min(npx_sta_mean[pos_pks]))


for chan in pos_pks:
    plt.plot(npx_sta_mean[chan,:], lw = 0.7)
fig2 = plt.plot(npx_sta_mean[central_chan,:]/scale, lw = 0.7, c = 'r')
plt.vlines(neg_pk, np.min(npx_sta_mean[central_chan]/scale), np.max(npx_sta_mean[pos_pks]), linestyle = '--', color = 'grey' )

##plot derivative of central chan: earliest positive peak is the point at which somatic negative peak starts decelerating!
#%%
for cell in range(0,5):
    sig_pks = np.where(DV_dist[cell,:,1] >= 0.2)[0]
    for chan in sig_pks:
        plt.scatter( DV_dist[cell,chan,0],  DV_dist[cell,chan,2],s = 3, )
#%%
for cell in range(0,5):
    sig_pks = np.where(DV_dist[cell,:,1] >= 0.2)[0]
    for chan in sig_pks:
        plt.scatter(   abs(DV_dist[cell,chan,2] /DV_dist[cell,chan,0]) ,DV_dist[cell,chan,0], s = 3, )
