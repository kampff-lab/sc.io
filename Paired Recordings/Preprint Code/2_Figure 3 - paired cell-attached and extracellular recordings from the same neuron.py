# -*- coding: utf-8 -*-
"""
Created on Fri Mar 02 10:37:23 2018

@author: Andre Marques-Smith
"""
#%%Imports
import numpy as np
import stfio
import matplotlib.pyplot as plt
#from wc_spike_extract_functions import *
import os
from collections import defaultdict
#%%
##Figure 3: Paired cell-attached and extra-cellular recordings from the same neuron.
#%%
listcells = os.listdir('E:/code/for analysis')
suffix = '_aligned_analysis'
#%%
paths = defaultdict(list)
cell = listcells.index('c46'+suffix)

aligned_directory = 'E:/code/for analysis/'+listcells[cell]
#%%
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
#%%%
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


#%%
color_patch = '#000000'
color_probe = '#0433FF'

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
#%%Figure 2B - find 500 random patch spikes and corresponding neuropixel waveforms on closest channel, plot them and overlay mean trace.
#patch traces
import random
trials_to_plot = 500
trials = random.sample(range(0, len(patch_spikes)), trials_to_plot)
#
width = 4.11
height = 2.96
y_offset = 2

patch_timebase = np.arange(-2,2, 4.0/len(patch_sta[0]) )
patch_trial_recordings = patch_sta[trials] 

for i in range(len(trials)):
    patch_trial_recordings[i]-= np.median(patch_trial_recordings[i], axis = 0)

patch_fig = plt.figure(figsize=cm2inch(width, height))
ax = patch_fig.add_subplot(111)


for i in range(len(trials)):    
    ax.plot(patch_timebase, patch_trial_recordings[i], alpha=0.1, color=color_patch, lw=0.3 )
plt.xlim(-2, 2)
plt.xticks([-2, -1, 0, 1, 2])
plt.ylim(-300, 100)
plt.yticks([-300, -200, -100, 0, 100])
ax.tick_params(labelbottom='off', labelleft='off', direction='out')
ax.vlines(0, -300, 100,  linestyles = 'dashed', lw = 0.3,)
ax.vlines(-1.75, -200, -100 )
ax.plot(patch_timebase, np.mean(patch_trial_recordings, axis = 0), lw=0.7, color='#C0C0C0')
plt.subplots_adjust(top = 0.94, bottom = 0.06, left = 0.04, right = 0.96)
#
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 2_materials/Fig_2C_top'+listcells[cell] + 'avg juxta.png', dpi = 1200)
plt.close()
#Npx traces

npx_timebase = np.arange(-2,2, 4.0/len(npx_sta_array[0,:,0]) )
npx_trial_recordings = 2.34375* npx_sta_array[181,:,trials] 
#

npx_fig = plt.figure(figsize=cm2inch(width, height))
ax = npx_fig.add_subplot(111)


for i in range(len(trials)):    
    ax.plot(npx_timebase, npx_trial_recordings[i], alpha=0.1, color=color_probe, lw=0.3 )
plt.xlim(-2, 2)
plt.xticks([-2, -1, 0, 1, 2])
plt.ylim(-250, 150)
plt.yticks([-250, -150, -50, 50, 150])
ax.tick_params(labelbottom='off', labelleft='off', direction='out')
ax.vlines(0, -250, 150,  linestyles = 'dashed', lw = 0.3,)
ax.vlines(-1.75, -150, -50 )
ax.plot(npx_timebase, np.mean(npx_trial_recordings, axis = 0), lw=0.7, color='#76D6FF')
plt.subplots_adjust(top = 0.94, bottom = 0.06, left = 0.04, right = 0.96)
#
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 2_materials/Fig_2C_bottom'+listcells[cell] + 'avg extra.png', dpi = 1200)
plt.close()

#%%Figure 2A - Pick a random 100ms segment of activity in patch-clamp and plot the corresponding segment on the 16 channels closest to the cell on the neuropixel probe.
#figure properties
width = 8.22
height = 2.96
npx_channel=181
length_seconds = 0.1
npx_samprate = 30000

m = len(patch_v)/float(len(npx_voltage[0]))
npx_buffer = int((npx_samprate * length_seconds)/2)
patch_buffer = int(npx_buffer * m)

patch_time = np.arange(0, length_seconds, length_seconds/(2*patch_buffer))
npx_time = np.arange(0, length_seconds, length_seconds/(2*npx_buffer))

#%%pick a spike
from random import randint


rand_spike = randint( np.where(patch_spikes>patch_buffer)[0][0], np.where(patch_spikes < (len(patch_v) - patch_buffer))  [-1][-1])
patch_sample = patch_spikes[rand_spike]
#%%
probe_segment = plt.figure(1, figsize=cm2inch(width, height))
ax = probe_segment.add_subplot(111)

ax.plot( npx_time, 2.34* npx_voltage[npx_channel, int(patch_sample/m) - npx_buffer : int(patch_sample/m) + npx_buffer  ] ,lw = 0.5, color = color_probe)
ax.tick_params(labelbottom='off', labelleft='off', direction='out')
ax.set_ylim([-200,100])
ax.set_yticks([-200, -100, 0, 100])
ax.set_xlim([0,0.1])
ax.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.1])
ax.vlines(0.005, -150,-50 )
plt.subplots_adjust(top = 0.94, bottom = 0.06, left = 0.04, right = 0.96)

#
patch_segment = plt.figure(2, figsize=cm2inch(width, height))
ax1 = patch_segment.add_subplot(111)

ax1.plot( patch_time, patch_v[ patch_sample - patch_buffer: patch_sample + patch_buffer], color=color_patch, lw=0.5)
ax1.tick_params(labelbottom='off', labelleft='off', direction='out')
ax1.set_ylim([-200,100])
ax1.set_yticks([-200, -100, 0, 100])
ax1.set_xlim([0,0.1])
ax1.vlines(0.005, -150,-50 )
plt.subplots_adjust(top = 0.94, bottom = 0.06, left = 0.04, right = 0.96)
ax1.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.1])
#%%
probe_segment.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 2_materials/Fig_2A_'+listcells[cell] + 'paired sample_probe.png', dpi = 1200 )
patch_segment.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 2_materials/Fig_2A_'+listcells[cell] + 'paired sample_patch.png', dpi = 1200 )

#%%Plot that same time segment for the 16 channels closest to the probe, respecting the channel layout; place black markers on samples of patch spike peak
import matplotlib.gridspec as gridspec
from matplotlib import colors
import matplotlib.cm as cmx


npx_samples = range(int(patch_sample/m) - npx_buffer, int(patch_sample/m) + npx_buffer)

npx_spikes = []
for i in patch_spikes:
    if int(i/m) in npx_samples:
        npx_spikes.append(int(i/m)-npx_samples[0])

single_trial_probe = np.empty((384,3000), dtype = float)
for chan in range(384):
    single_trial_probe[chan] = 2.34375*npx_voltage[chan, int(patch_spikes[rand_spike]/m - 1500) : int(patch_spikes[rand_spike]/m+1500)]

num_rows_side = 7

rc_channel = np.where(npx_mapcsv == npx_channel)
sub_map = npx_mapcsv[rc_channel[0][0]-num_rows_side:rc_channel[0][0]+num_rows_side+1,:]
electrodes = sub_map[sub_map > 0 ]

width=8.23
height = 19.27

fig = plt.figure(1, figsize=cm2inch(width,height))

origin_cmap= plt.get_cmap('cool')
cm=origin_cmap
cNorm=colors.Normalize(vmin= npx_channel - 1 - 7*2, vmax= npx_channel + 1+ 7*2 )#np.max(electrodes))
scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)

for chan in electrodes:
    if chan in range(npx_channel - 1 - 7*2, npx_channel + 1+ 7*2):
        colorVal=scalarMap.to_rgba( chan )
        gridspec.GridSpec(num_rows_side*2+1,4)  
        rc_channel = np.where(sub_map == chan)
        plt.subplot2grid((num_rows_side*2+1,4), (rc_channel[0][0], rc_channel[1][0]))
        plt.subplots_adjust(top = 1.0, bottom = 0.0, left = 0.0, right = 1.0, wspace = 0.1, hspace = 0.0)
        plt.axis('off')
        plt.ylim(-200, 100)
        #colorVal=scalarMap.to_rgba( np.max(npx_sta_mean[chan,70:90])-np.min(npx_sta_mean[chan,50:70]) )
        plt.plot(single_trial_probe[chan,:], lw = 0.2, color = colorVal)
        for spike in range(len(npx_spikes)):
            plt.vlines(npx_spikes[spike], 80, 100, color = 'k', lw=0.7)
    else:
# set up subplot grid
    #colorVal=scalarMap.to_rgba( rankpk2pk(npx_sta_mean, electrodes, chan) )
        gridspec.GridSpec(num_rows_side*2+1,4)  
        rc_channel = np.where(sub_map == chan)
        plt.subplot2grid((num_rows_side*2+1,4), (rc_channel[0][0], rc_channel[1][0]))
        plt.subplots_adjust(top = 1.0, bottom = 0.0, left = 0.0, right = 1.0, wspace = 0.1, hspace = 0.0)
        plt.axis('off')
        plt.ylim(-200, 100)
        #colorVal=scalarMap.to_rgba( np.max(npx_sta_mean[chan,70:90])-np.min(npx_sta_mean[chan,50:70]) )
        plt.plot(single_trial_probe[chan,:], lw = 0.2, color = color_probe)
        for spike in range(len(npx_spikes)):
            plt.vlines(npx_spikes[spike], 80, 100, color = 'k', lw=0.7)
gridspec.GridSpec(num_rows_side*2+1,4)
plt.subplot2grid((num_rows_side*2+1,4), (14, 0))
plt.axis('off')
plt.ylim(-200,100)
plt.xlim(0,100)
plt.vlines(25,-50,50)

#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 2_materials/Fig_2C_'+listcells[cell] + 'sample_single_trial_probe.png', dpi = 1200 )
#%%Fig 2D - Plot the average waveform of ground truth spikes for the 16 channels closest to the cell, colour-coded same way as 100ms segments in Fig2C.

def plot_avg_multichan_extra(npx_sta_mean, central_chan, npx_mapcsv, width, height, num_rows_side = 7, yoffset=1):
    def triggerline(x):
        if x is not None:
            ylim = plt.ylim()
            plt.vlines(x, ylim[0], ylim[1], linewidth = 0.3, linestyles='dashed')
        
    rc_channel = np.where(npx_mapcsv == central_chan)
    sub_map = npx_mapcsv[rc_channel[0][0]-num_rows_side:rc_channel[0][0]+num_rows_side+1,:]
    electrodes = sub_map[sub_map > 0 ]
    
    origin_cmap= plt.get_cmap('cool')
    #shrunk_map = shiftedColorMap(origin_cmap, start=0.5, midpoint=0.75, stop=1, name='shrunk')
    cm=origin_cmap
    cNorm=colors.Normalize(vmin=np.min(electrodes), vmax= np.max(electrodes) )#np.max(electrodes))
    scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)

    timebase = np.arange(-2,2, 4.0/len(npx_sta_mean[0]) )
    fig = plt.figure(figsize=cm2inch(width,height))
    ax = fig.add_subplot(111)
    for electrode in electrodes:
        
        colorVal=scalarMap.to_rgba( electrode )
        ax.plot(timebase, npx_sta_mean[ electrode ,:], color=colorVal, lw = 0.8, alpha=0.8)
        plt.subplots_adjust(top = 0.94, bottom = 0.06, left = 0.04, right = 0.96)
        plt.xlim(-2, 2)
        plt.xticks([-2,-1,0,1,2])
        ax.tick_params(labelbottom='off', labelleft='off', direction='out')
        plt.ylim(-200, 100 )
        #plt.axis('off')
        #plt.ylabel('Voltage (uV)', fontsize=20)
        #plt.xlabel('Time (ms)',fontsize=20)
        #plt.xticks(fontsize=15)
        #plt.yticks(fontsize=15)
        #plt.tick_params(right = 'on', top = 'on', bottom = 'on', direction = 'in')
    triggerline(0)
    plt.ylim(-200,100)
    plt.vlines(-1.75, -150,-50)
        
            
#%%
plot_avg_multichan_extra(npx_sta_mean, 181, npx_mapcsv, width = 4.11, height = 6.42)
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 2_materials/'+listcells[cell] + 'multiple_extra_chans.png', dpi = 1200 )
#%%Fig 2E - Same traces as 2D, but plotted according to the probe layout.
def plot_avg_multichan_extra_probe(npx_sta_mean, central_chan, npx_mapcsv, width, height, num_rows_side = 7, yoffset=1):
    rc_channel = np.where(npx_mapcsv == central_chan)
    sub_map = npx_mapcsv[rc_channel[0][0]-num_rows_side:rc_channel[0][0]+num_rows_side+1,:]
    electrodes = sub_map[sub_map > 0 ]

    origin_cmap= plt.get_cmap('cool')
    #shrunk_map = shiftedColorMap(origin_cmap, start=0.5, midpoint=0.65, stop=0.85, name='shrunk')
    cm=origin_cmap
    cNorm=colors.Normalize(vmin=np.min(electrodes), vmax= np.max(electrodes) )#np.max(electrodes))
    scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)
    plt.figure(figsize=cm2inch(width,height))
    
    for chan in electrodes:
    # set up subplot grid
        colorVal=scalarMap.to_rgba( chan )
        gridspec.GridSpec(15,4)  
        rc_channel = np.where(sub_map == chan)
        plt.subplot2grid((15,4), (rc_channel[0][0], rc_channel[1][0]))
        plt.subplots_adjust(top = 1.0, bottom = 0.0, left = 0.0, right = 1.0, wspace = 0.1, hspace = 0.0)
        plt.axis('off')
        plt.ylim(-200, 100 )
        #colorVal=scalarMap.to_rgba( np.max(npx_sta_mean[chan,70:90])-np.min(npx_sta_mean[chan,50:70]) )
        plt.plot(npx_sta_mean[chan,:], lw = 0.6, color = colorVal)
    gridspec.GridSpec(num_rows_side*2+1,4)
    plt.subplot2grid((num_rows_side*2+1,4), (14, 0))
    plt.axis('off')
    plt.ylim(-200,100)
    plt.xlim(0,100)
    plt.vlines(25,-50,50)
#%%
plot_avg_multichan_extra_probe(npx_sta_mean, 181, npx_mapcsv, width = 4.11, height = 12.5)
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 2_materials/'+listcells[cell] + 'multiple_extra_chans_probe.png', dpi = 1200 )