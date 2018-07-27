# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 11:37:29 2018

@author: Andre Marques-Smith
"""
#%%Imports
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
 
#%%Data frames and structures
listcells = os.listdir('E:/code/for analysis')
suffix = '_aligned_analysis'

data_summary = pandas.read_excel('C:/Users/Andre Marques-Smith/Dropbox/Patch-Clamp Summary 2.0.xlsx')
cells_above_10uV = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 10].tolist()
#excluded_cells = [2, 6, 9, 10, 12, 21, 23]
excluded_cells = [2, 5, 9, 12, 21, 23]
cells_to_analyse = [ cells_above_10uV[i] for i in range(len(cells_above_10uV)) if i not in excluded_cells]

#%%Code for Joy division plots in Fig 7A and 7B - will generate for every cell.
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
    npx_chan_peaks = np.load(paths['npx_chan_peaks'][0])#, dtype = 'float64')
    #
    m = len(patch_v)/float(len(npx_voltage[0,:]))
    patch_spikes_npx = np.asarray([int(patch_spikes[i] / m) for i in range(len(patch_spikes))])
#added int
    central_chan = int(npx_chan_peaks[0,0])  

#These plots only use one probe column, therefore we begin by finding which column a cell was, and which channels are in that column in order to plot them.
    #What we're plotting is the first derivative of voltage over time, normalised to the channel closest to the soma/
    cellcol = np.where(npx_mapcsv == central_chan)[1][0]
    if cellcol == 0:
        channels = np.arange(0,384,4)
    elif cellcol == 1:
        channels = np.arange(2,385,4)
    elif cellcol == 2:
        channels = np.arange(1,385,4)
    elif cellcol == 3:
        channels = np.arange(3,385,4)
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])
    shader = plt.subplot(gs[0])
    shader.imshow(np.diff(npx_sta_mean[channels,:])/abs(np.min(np.diff(npx_sta_mean[channels,:]) )),cmap='seismic_r', interpolation='none', aspect = 3,vmin = -1, vmax = 1)#, norm = MidpointNormalize(midpoint=0.))# vmin = -1, vmax = 1)
    #ax.set_title('%s'%(cell_id ), y = 1.0, fontweight = 'bold' )
#        cbar = fig.colorbar(shader, orientation = 'horizontal', shrink = 0.75)
    #plt.text(3,96 - np.where(npx_mapcsv == central_chan)[0][0]/2,'cell' )
    cellrow = 96 - np.where(npx_mapcsv == central_chan)[0][0]/2 
    shader.set_yticks([cellrow-15, cellrow-1, cellrow +15])
    shader.set_yticklabels(['%s-->'%(central_chan - 64), '%s'%(central_chan), '%s-->'%(central_chan+64)])
    #plt.hlines(96 - np.where(npx_mapcsv == central_chan)[0][0]/2 ,0,120 )

    #ax.set_yticklabels(['384', '360', '336', '312', '288', '264', '240', '216', '192', '168', '144', '120', '96', '72', '48', '24', '0'])
    shader.set_ylabel('Channel', fontweight='bold')
    shader.set_xlabel('Time (ms)', fontweight = 'bold')
    shader.set_xticks([0, 30, 60, 90, 119])
    shader.set_xticklabels(['-2.0', '-1.0',  '0',  '1.0',  '2.0'])
    shader.spines['left'].set_linewidth(1)
    shader.spines['right'].set_linewidth(1)
    shader.spines['top'].set_linewidth(1)
    shader.spines['bottom'].set_linewidth(1)
    shader.tick_params(axis='both', direction='in')
    shader.hlines(cellrow-15,0,119,linestyle = 'dashed', lw = 0.4  )
    shader.hlines(cellrow+15,0,119,linestyle = 'dashed', lw = 0.4  )
    plt.gca().invert_yaxis()
    #
    joydivision = plt.subplot(gs[1])
    chanlist = np.arange(central_chan - 64,central_chan+64,4)
    range_joy = np.min(npx_sta_mean[:,:]), np.max(npx_sta_mean[:,:])
    for chan in np.arange(central_chan - 64,central_chan+64,4):
        if range_joy[0] < (np.min(npx_sta_mean[chanlist[0],:]) +chanlist[0] - central_chan):
        
            if chan == central_chan:
                joydivision.plot( npx_sta_mean[chan,:]+ npx_chan_peaks[0,1]/64*(chan - central_chan) , color = 'm', lw=2, alpha = 0.9)
            else:
                joydivision.plot( npx_sta_mean[chan,:]+ npx_chan_peaks[0,1]/64*(chan - central_chan) , color = 'w', lw=1)
                joydivision.set_facecolor('k')
                joydivision.set_xticks([])
                joydivision.set_xticklabels([''])
                joydivision.set_yticks([npx_chan_peaks[0,1]/64*(chanlist[0] - central_chan), npx_chan_peaks[0,1]/64*(chanlist[-1] - central_chan)  ])
                #joydivision.set_yticks([ np.median(npx_sta_mean[chanlist[0],:] -64  ), np.median(npx_sta_mean[central_chan,:]), np.median(npx_sta_mean[chanlist[-1],:]+62  ) ])
                joydivision.set_yticklabels(['-->', '-->' ])
                joydivision.tick_params(axis='both', direction='in', color = 'w')
                joydivision.set_xticks([0, 30, 60, 90, 120])
                joydivision.set_xticklabels(['-2.0', '-1.0',  '0',  '1.0',  '2.0'])
                joydivision.spines['left'].set_linewidth(0.1)
                joydivision.spines['right'].set_linewidth(0.1)
                joydivision.spines['top'].set_linewidth(0.1)
                joydivision.spines['bottom'].set_linewidth(0.1)
                joydivision.set_xlabel('Time (ms)', fontweight = 'bold')
        else:
            if chan == central_chan:
                joydivision.plot( npx_sta_mean[chan,:]+ chan - central_chan , color = 'm', lw=1, alpha = 0.8)
            else:
                joydivision.plot( npx_sta_mean[chan,:]+ chan - central_chan , color = 'w', lw=1)
                joydivision.set_facecolor('k')
                joydivision.set_xticks([])
                joydivision.set_xticklabels([''])
                joydivision.set_yticks([ np.median(npx_sta_mean[chanlist[0],:] -64  ), np.median(npx_sta_mean[central_chan,:]), np.median(npx_sta_mean[chanlist[-1],:]+62  ) ])
                joydivision.set_yticklabels(['-->', '', '-->' ])
                joydivision.tick_params(axis='both', direction='in', color = 'w')
                joydivision.set_xticks([0, 30, 60, 90, 120])
                joydivision.set_xticklabels(['-2.0', '-1.0',  '0',  '1.0',  '2.0'])
                joydivision.spines['left'].set_linewidth(0.1)
                joydivision.spines['right'].set_linewidth(0.1)
                joydivision.spines['top'].set_linewidth(0.1)
                joydivision.spines['bottom'].set_linewidth(0.1)
                joydivision.set_xlabel('Time (ms)', fontweight = 'bold')
            

    fig.suptitle(data_summary.loc[cell]['Cell'], fontsize=14, fontweight = 'bold')    
    plt.subplots_adjust(wspace=0.16)        
    #
    plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 7_materials/' + 'Fig 7_shaderjoy_%s.png' %(data_summary.loc[cell]['Cell']), dpi = 600)
    plt.close()
#%%
