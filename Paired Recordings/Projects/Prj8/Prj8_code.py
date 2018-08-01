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
import math
import matplotlib.gridspec as gridspec
from matplotlib import colors
import matplotlib.cm as cmx 
from scipy.interpolate import interp1d
#%%Data frames and structures
listcells = os.listdir('E:/code/for analysis')
suffix = '_aligned_analysis'

data_summary = pandas.read_excel('C:/Users/Andre Marques-Smith/Dropbox/Data Summary.xlsx')
cells_above_10uV = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 10].tolist()
#excluded_cells = [2, 6, 9, 10, 12, 21, 23]
excluded_cells = [2, 5, 9, 12, 21, 23]
cells_to_analyse = [ cells_above_10uV[i] for i in range(len(cells_above_10uV)) if i not in excluded_cells]
n_cells = len( cells_to_analyse)

backprop = 999*np.ones((n_cells, 49,4), dtype = 'float')
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

    #npx_voltage = np.memmap( paths['npx_v'][0], mode = 'c', dtype=np.int16 )
    #npx_voltage = npx_voltage.reshape((expt_meta['npx'][0][0], expt_meta['npx'][0][1]))
    #
    npx_sta_array = np.load(paths['npx_sta_array'][0])#, dtype = 'int16')
    npx_sta_mean = np.load(paths['npx_sta_mean'][0])#, dtype = 'float64')
    #npx_chan_peaks = np.load(paths['npx_chan_peaks'][0])#, dtype = 'float64')
    npx_chan_peaks = [np.max(npx_sta_mean[i]) - np.min(npx_sta_mean[i]) for i in range(384)]
    
    central_chan = np.where(npx_chan_peaks ==np.max(npx_chan_peaks))[0][0]

    
    nc = cells_to_analyse.index(cell)
    #
    col1 = [43,59]
    col2 = [11, 27]
    if chanMapcsv[central_chan,1] in col1:
        a = np.where(chanMapcsv[:,1] == col1[0])[0]
        b = np.where(chanMapcsv[:,1] == col1[1])[0]
    else:
        a = np.where(chanMapcsv[:,1] == col2[0])[0]
        b = np.where(chanMapcsv[:,1] == col2[1])[0]
    c = np.concatenate((a,b), axis=0)
    c = np.sort(-c)
    col_channels = -c
    #%
    col_distances = np.empty((192), dtype = 'float')
    x0,y0 = chanMapcsv[np.where(chanMapcsv[:,0] == central_chan),1:3][0][0]
    for chan in range(len(col_channels)):
        x1,y1 = chanMapcsv[np.where(chanMapcsv[:,0] == col_channels[chan]),1:3][0][0]
        col_distances[chan] = int(math.sqrt( (x0 - x1)**2 + (y0-y1)**2) )
        if y1 < y0:
            col_distances[chan] *= -1
        
    central_pk_a = np.max(npx_sta_mean[central_chan,40:80]) - np.min(npx_sta_mean[central_chan,40:80])
    central_pk_t = 40+np.argmin(npx_sta_mean[central_chan,40:80])
    
    
    #Create colormap with first derivative of voltage/t for each cell; plot traces by side and find earliest minima
    shader = plt.subplot2grid((25,2), (0,0), rowspan=25)
    shader.imshow(np.diff(npx_sta_mean[col_channels[abs(col_distances)<=240],:])/abs(np.min(np.diff(npx_sta_mean[col_channels[abs(col_distances)<=240],:]) )),cmap='seismic_r', interpolation='none', aspect = 10,vmin = -1, vmax = 1)#, norm = MidpointNormalize(midpoint=0.))# vmin = -1, vmax = 1)
    shader.set_yticks([0, 6, 12, 18, 24])
    shader.set_yticklabels(['+0.24', '+0.12', '0', '-0.12','-0.24'])
    shader.set_ylabel('Distance to soma (mm)', fontweight='bold')
    shader.set_xlabel('Time (ms)', fontweight = 'bold')
    shader.set_xticks([30, 60, 90])
    shader.set_xticklabels(['-1.0',  '0',  '1.0'])
    det_thresh = 10.0
    for chan in col_channels[abs(col_distances)<=480]:#range(len(subset_240)):
        chan_pk2pk = (np.max(npx_sta_mean[chan]) - np.min(npx_sta_mean[chan]))/central_pk_a
        chan_pk_a = np.max(npx_sta_mean[chan,40:80]) - np.min(npx_sta_mean[chan,40:80])
        chan_pk_t = 40+np.argmin(npx_sta_mean[chan,40:80])
        if chan_pk2pk * central_pk_a >= det_thresh:
            color_trc = 'k'
            alph = 1.0
        else:
            color_trc = 'grey'
            alph = 0.75
        if chan in col_channels[abs(col_distances)<=240]:
            joydivision = plt.subplot2grid((25,2), (list(col_channels[abs(col_distances)<=240]).index(chan),1), rowspan=1)
            joydivision.plot(npx_sta_mean[chan], color = color_trc, lw = 0.8, alpha = alph)
            joydivision.axis('off')
    
            if chan_pk2pk * central_pk_a > det_thresh:
                a = npx_sta_mean[chan]- int(0.5*np.min(npx_sta_mean[chan,40:80]))
                b = np.where(np.sign(a)==-1)[0]
                if len(b)>0:
                    b=b[0]
                    c = np.zeros((121), dtype = 'float')
                    c[1:] = np.diff(a)
                    d = b+np.where(c[b:chan_pk_t]>0)[0]
                    if len(d) >1:
                        d = d[0]
                        chan_pk_t = b+np.argmin(a[b:d])
                        #backprop_prof[nc, np.where(col_channels==chan)[0][0],2] = 1000/30000.0 * (chan_pk_t - central_pk_t)
                joydivision.scatter(chan_pk_t , npx_sta_mean[chan, chan_pk_t] ,s=10, alpha = 0.5 , c = 'm')
            backprop[nc, list(col_channels[abs(col_distances)<=480]).index(chan), 1] = chan_pk_a / central_pk_a # Normalised amplitude
            backprop[nc, list(col_channels[abs(col_distances)<=480]).index(chan), 2] = 1000/30000.0 * (chan_pk_t - central_pk_t) # latency between channel negative peak and somatic negative peak
            backprop[nc, list(col_channels[abs(col_distances)<=480]).index(chan), 3] = chan_pk_a
        else:
            backprop[nc, list(col_channels[abs(col_distances)<=480]).index(chan), 1] = chan_pk_a / central_pk_a
            backprop[nc, list(col_channels[abs(col_distances)<=480]).index(chan), 3] = chan_pk_a
            
        plt.tight_layout
    plt.suptitle(data_summary.loc[cells_to_analyse[nc]]['Cell'], fontsize=14, fontweight = 'bold')
    
    plt.savefig( 'E:/repos/sc.io/Paired Recordings/Projects/Prj8/cells/supp_fig_7.1_%s.png' %(data_summary.loc[cells_to_analyse[nc]]['Cell']), dpi = 300)
    plt.close()
    backprop[nc,:,0] = col_distances[abs(col_distances)<=480]
    #%%EAP propagation trajectory figure - plotting only cells where EAP could be tracked at more than 3 points (amplitude > 10uV) and only ±240 um from soma
    origin_cmap= plt.get_cmap('tab20')
    cm=origin_cmap
    cNorm=colors.Normalize(vmin= 0, vmax= 21 )
    scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)
    labels = []
    
    for cell in range(21):
        to_plot = 12+np.where(backprop[cell,12:36,3]>9.5)[0]
        if len(to_plot)>3:
            colorVal=scalarMap.to_rgba( cell )
            plt.plot( backprop[cell, to_plot, 2], backprop[cell, to_plot, 0]/1000.0, color = colorVal)
            labels.append(data_summary.loc[cells_to_analyse[cell]]['Cell'])
        
        plt.xlim(-0.2,0.8)
        plt.yticks([-0.24, -0.120, 0.0, 0.12, 0.24])
        plt.ylim(-0.250,0.250)
        plt.legend(labels, loc=7, frameon=False)
    plt.vlines(0,-0.250,0.250, linestyle='--', lw = 1)
    plt.hlines(0,-0.2,0.6, linestyle='--', lw = 1)
    plt.xlabel('Time relative to soma (ms)', fontweight = 'bold')
    plt.ylabel('Distance to soma (mm)', fontweight = 'bold')
    
    plt.savefig('E:/repos/sc.io/Paired Recordings/Projects/Prj8/Fig_7B_Trajectories.png')
    plt.close()
#%%EAP amplitude propagation profile - interpolate average propagation and find intersection with y = 0.5, 0.25 and 0.12
backprop_interp = np.empty((960,2), dtype = 'float')
backprop_interp[:,0] = np.arange(480,-480,-1)
f = interp1d(np.average(backprop[:,:,0], axis = 0), np.average(backprop[:,:,1], axis = 0) )
backprop_interp[:,1] = f(backprop_interp[:,0])

windows = np.zeros((3,3), dtype = 'float')
windows[0:3,0] = 0.5, 0.25, 0.12
for t in windows[:,0]:
    windows[np.where(windows[:,0]==t),2] = backprop_interp[np.where(abs(backprop_interp[:,1] - t) <=0.005)[0][0],0]
    windows[np.where(windows[:,0]==t),1] = backprop_interp[np.where(abs(backprop_interp[:,1] - t) <=0.005)[0][-1],0]
#
origin_cmap= plt.get_cmap('tab20')
cm=origin_cmap
cNorm=colors.Normalize(vmin= 0, vmax= 21 )
scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)

for cell in range(n_cells):
    colorVal=scalarMap.to_rgba( cell )
    plt.scatter(backprop[cell,:,1], backprop[cell,:,0]/1000.0,  c = colorVal, s = 3, alpha = 0.2)
    
plt.plot(np.average(backprop[:,:,1], axis = 0),np.average(backprop[:,:,0], axis = 0)/1000.0,  c = 'k')
plt.ylabel('Distance to soma (mm)',fontweight='bold')
plt.yticks([-0.48,0.0,0.48])
plt.xticks([0.0, 0.12, 0.25, 0.5, 0.75, 1.0])
plt.xlabel('Normalised Peak-Peak Amplitude', fontweight='bold')
plt.text(0.15,1.0, 'Average (n = 21)',fontweight='bold', color='k')

tmap = plt.get_cmap('cool')
#cm=origin_cmap
cNormt=colors.Normalize(vmin= 0.12, vmax= 0.5 )
scalarMap= cmx.ScalarMappable(norm=cNormt,cmap=tmap)
for t in windows[:,0]:
    colorVal=scalarMap.to_rgba( t )
    plt.hlines(windows[np.where(windows[:,0]==t),1]/1000.0, -0.05, windows[np.where(windows[:,0]==t),0]  , color = colorVal)
    plt.hlines(windows[np.where(windows[:,0]==t),2]/1000.0, -0.05, windows[np.where(windows[:,0]==t),0]  , color = colorVal)
    plt.vlines(t, windows[np.where(windows[:,0]==t),1]/1000.0, windows[np.where(windows[:,0]==t),2]/1000.0   ,  color = colorVal)
    plt.vlines(t, -0.5, windows[np.where(windows[:,0]==t),2]/1000.0   ,  color = colorVal, linestyle='--', alpha = 0.7)
    plt.text( t+0.02,windows[np.where(windows[:,0]==t),2]/1000.0 + 0.01, '%s mm' %( round(windows[np.where(windows[:,0]==t),2][0][0]/1000.0,2) ), style = 'italic',  color = colorVal)
    plt.text( t+0.02,windows[np.where(windows[:,0]==t),1]/1000.0 - 0.03, '%s mm' %( round(windows[np.where(windows[:,0]==t),1][0][0]/1000.0,2) ), style = 'italic',  color = colorVal)
plt.xlim(-0.05,1.05)
plt.ylim(-0.5,0.5)
#%%
plt.savefig('E:/repos/sc.io/Paired Recordings/Projects/Prj8/Fig_7A_EAP_propagation_profile.png')
plt.close()
#%%





###From here below, testing
#%%
origin_cmap= plt.get_cmap('rainbow')
cm=origin_cmap
cNorm=colors.Normalize(vmin= 0, vmax= n_cells )
scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)

for cell in range(n_cells):
    colorVal=scalarMap.to_rgba( cell )
        #if backprop_prof120[cell,6,3] >=50:
    plt.plot(backprop_prof120[cell,:,2], backprop_prof120[cell,:,0], color = colorVal, alpha = 0.7)
    plt.xlim(-1.0,1.0)
    #%%
    backprop_array = np.empty((len(good_spikes), len(subset_240),2), dtype = 'float')
    
    for spike in range(len(good_spikes)):
        for chan in range(len(subset_240)):
            soma_pk = np.max(npx_sta_array[central_chan,50:70,spike]) - np.min(npx_sta_array[central_chan,50:70,spike])
            soma_t = 50+np.argmin(npx_sta_array[central_chan,50:70,spike])
            
            chan_pk = np.max(npx_sta_array[subset_240[chan],50:70,spike]) - np.min(npx_sta_array[subset_240[chan],:,spike])
            chan_t = 50+np.argmin(npx_sta_array[subset_240[chan],50:70,spike])
            
            backprop_array[spike, chan, 0] = chan_pk/float(soma_pk)
            backprop_array[spike, chan, 1] = (chan_t - soma_t) * 1000/30000.0
    #%%
    for spike in range(len(good_spikes)):
        plt.plot(backprop_array[spike,6:19,1], backprop_prof240[20,6:19,0], color = 'grey', alpha = 0.1, lw = 0.3)
    plt.plot( np.mean(backprop_array[:,6:19,1],axis = 0), backprop_prof240[20,6:19,0], color = 'b' )
    #%%
dot_pr = np.empty((len(patch_spikes), len(subset_240)), dtype = float)

for spike in range(len(patch_spikes)):
    for chan in col_channels[abs(col_distances)<=240]:
        dot_pr[spike, list(col_channels[abs(col_distances)<=240]).index(chan) ] = angle_between( npx_sta_mean[chan], npx_sta_array[chan,:,spike]  )
        #math.degrees(angle_between( npx_sta_mean[list(col_channels[abs(col_distances)<=240]).index(chan),40:80] , npx_sta_array[ list(col_channels[abs(col_distances)<=240]).index(chan),40:80,spike]))
#        math.degrees( angle_between( npx_sta_mean[ list(col_channels[abs(col_distances)<=240).index(chan)] ], npx_sta_array[list(col_channels[abs(col_distances)<=240).index(chan)],:,spike] ))
#        np.dot(npx_sta_mean[list(col_channels[abs(col_distances)<=240].index(chan)),40:80], npx_sta_array[list(col_channels[abs(col_distances)<=240].index(chan)),40:80,spike] ) / np.dot(npx_sta_mean[subset_240[chan],40:80], npx_sta_mean[subset_240[chan],40:80])
    #%%
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
dot_pr_multi = np.empty((len(patch_spikes), 1), dtype = float)

for spike in range(len(patch_spikes)):
    dot_pr_multi[spike] = math.degrees(angle_between(np.concatenate(npx_sta_mean[col_channels[abs(col_distances)<=240]]) , np.concatenate(npx_sta_array[col_channels[abs(col_distances)<=240],:,spike] )))#np.dot(np.concatenate(npx_sta_mean[subset_240,45:60]) , np.concatenate(npx_sta_array[subset_240,45:60,spike] )) / np.dot(np.concatenate(npx_sta_mean[subset_240, 45:60]), np.concatenate(npx_sta_mean[subset_240,45:60]) )
#%%

    #%%
    good_spikes = np.where(dot_pr_multi <= 40)[0]
    
    fig, avg = plt.subplots(25,3, sharey='row')
    #avg = avg.ravel()
    #avg = plt.subplot2grid((25,2),(0,0))
    for chan in range(len(subset_240)):
        #avg.subplot2grid((25,2), (chan,0))
        avg[chan,0].plot(npx_sta_mean[subset_240[chan]])
        avg[chan,0].axis('off')

        for spike in good_spikes:
            avg[chan,1].plot(npx_sta_array[subset_240[chan],:,spike], alpha = 0.5, lw = 0.7)
            avg[chan,1].axis('off')
            
            #
        avg[chan,2].plot(np.mean(npx_sta_array[subset_240[chan],:,good_spikes], axis = 0), color = 'orange' )
        avg[chan,2].axis('off')
        avg[chan,1].text(80, 10,'ch %s' %(subset_240[chan]), alpha = 0.5, fontweight = 'bold', fontsize = 7)
    #%%
    subs = subset_240[12:20]
    fig, avg = plt.subplots(8,3, sharey='row')
    #avg = avg.ravel()
    #avg = plt.subplot2grid((25,2),(0,0))
    for chan in range(len(subs)):
        #avg.subplot2grid((25,2), (chan,0))
        avg[chan,0].plot(npx_sta_mean[subs[chan]])
        avg[chan,0].axis('off')

        for spike in good_spikes:
            avg[chan,1].plot(npx_sta_array[subs[chan],:,spike], alpha = 0.5, lw = 0.7)
            avg[chan,1].axis('off')
            
            #
        avg[chan,2].plot(np.mean(npx_sta_array[subs[chan],:,good_spikes], axis = 0), color = 'orange' )
        avg[chan,2].axis('off')
        avg[chan,1].text(80, 10,'ch %s' %(subs[chan]), alpha = 0.5, fontweight = 'bold', fontsize = 7)

    #%%
    origin_cmap= plt.get_cmap('seismic')
    cm=origin_cmap
    cNorm=colors.Normalize(vmin= 0, vmax= 25 )
    scalarMap= cmx.ScalarMappable(norm=cNorm,cmap=cm)

    for chan in range(25):
        colorVal=scalarMap.to_rgba( chan )
        plt.hist(dot_pr[:,chan] , bins = np.arange(-3.14,3.15,0.03),alpha = 0.7, color = colorVal)
    plt.legend(col_channels[abs(col_distances)<=240])
    #%%
    npx_voltage = np.memmap( paths['npx_v'][0], mode = 'c', dtype=np.int16 )
    npx_voltage = npx_voltage.reshape((expt_meta['npx'][0][0], expt_meta['npx'][0][1]))