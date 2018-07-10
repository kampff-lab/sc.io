#%%Imports

import pandas
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict
import random
import math
import scipy.stats as stats
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn import metrics
#%%
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
    
def rolling_mean(a, w):
    r = np.empty(a.shape)
    r.fill(np.nan)
    for i in range(w - 1, a.shape[0]):
        r[i] = np.mean(a[(i-w+1):i+1])
    return r

def rolling_median(a, w):
    r = np.empty(a.shape)
    r.fill(np.nan)
    for i in range(w - 1, a.shape[0]):
        r[i] = np.median(a[(i-w+1):i+1])
    return r

def rolling_std(a, w):
    r = np.empty(a.shape)
    r.fill(np.nan)
    for i in range(w - 1, a.shape[0]):
        r[i] = np.std(a[(i-w+1):i+1])
    return r
#%%Data frames and structures
listcells = os.listdir('E:/code/for analysis')
suffix = '_aligned_analysis'

data_summary = pandas.read_excel('C:/Users/Andre Marques-Smith/Dropbox/Patch-Clamp Summary 2.0.xlsx')

cells_above_10uV = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 10].tolist()
#excluded_cells = [2, 6, 9, 10, 12, 21, 23]
excluded_cells = [2, 5, 9, 12, 21, 23]
cells_above_50uV = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 50].tolist()
cells_to_analyse = cells_above_50uV
#Excluding c27 - too many spikes confounding GT spike
cells_to_analyse.pop(5)
#%%Data structures
parsed_spikes = np.empty((len(cells_to_analyse), 2), dtype = object)
summary_stats = np.empty((len(cells_to_analyse), 8, 4), dtype = float)      #8 features, mean, std, cv
pcent_spikes = np.empty((len(cells_to_analyse),1), dtype = float)

#%%
output_dir = 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/repro/fig6/'
for cell in cells_to_analyse:
    paths = defaultdict(list)
    cell_idx = listcells.index(data_summary.loc[cell]['Cell']+suffix)

    aligned_directory = 'E:/code/for analysis/'+listcells[cell_idx]
    
    row = cells_to_analyse.index(cell)

#paths = defaultdict(list)
#use cell id or type in string ('c45') below
#cell = listcells.index('c14'+suffix)

    #aligned_directory = 'E:/code/for analysis/'+listcells[cell]
    #
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
    npx_voltage = np.memmap( paths['npx_v'][0], mode = 'c', dtype=np.int16 )
    npx_voltage = npx_voltage.reshape((expt_meta['npx'][0][0], expt_meta['npx'][0][1]))
    #
    npx_sta_array = np.load(paths['npx_sta_array'][0])#, dtype = 'int16')
    npx_sta_mean = np.load(paths['npx_sta_mean'][0])#, dtype = 'float64')
    npx_chan_peaks = np.load(paths['npx_chan_peaks'][0])#, dtype = 'float64')
    
    central_chan = int(npx_chan_peaks[0,0])
#
    m = len(patch_v)/float(len(npx_voltage[0,:]))
    #%
    npx_gt_spikes = [int(i/m) for i in patch_spikes]
    
    gt_pk_lag = np.argmin(npx_sta_mean[central_chan])
    
    samprate = 30000.0
    
    
    #Extract spike features of ground truth spikes 
    gt_spike_features = np.empty((len(npx_gt_spikes), 9), dtype=float)

    for spike in range(len(npx_gt_spikes)):
        t_ia = gt_pk_lag
        t_ib = t_ia+1
        baseline = np.median(npx_sta_array[central_chan,:,spike])
        #t_ia is the first point of threshold crossing for a spike (when it goes 'down')
        #t_ib is the last point (when it goes 'up')
        #To determine t_ia, I use npx_sta_array, which has the extracellular snippets for every patch spike. I go to the time point where, on average, the patch spike peak
        #was, and start treading back one sample at a time, until I cross the threshold (baseline). That procedure is repeated in the opposite direction to get t_ib.
        
        while npx_sta_array[central_chan, t_ia, spike] <= baseline:
            t_ia -= 1
            if t_ia < 30:
                pass
            
        while npx_sta_array[central_chan, t_ib, spike] <= baseline:
            t_ib += 1
            if t_ib > 80:
                pass
        
        neg_pk_i = t_ia + np.argmin( npx_sta_array[central_chan,t_ia:t_ib,spike] )
        neg_pk_v = 2.34375 * npx_sta_array[central_chan, neg_pk_i, spike] - baseline
        
        pos_pk_i = t_ib + np.argmax( npx_sta_array[central_chan,t_ib:t_ib+20,spike] )
        pos_pk_v = 2.34375 * npx_sta_array[central_chan, pos_pk_i, spike] - baseline
        
        pk2pk = pos_pk_v - neg_pk_v
        
        duration = (pos_pk_i - neg_pk_i) * 1000/samprate
        if duration <0:
            if spike >=1:
                duration = gt_spike_features[spike-1, 2]
        
        half_amp = (neg_pk_v/2.34375 + npx_sta_array[central_chan,t_ia,spike])/2
        half_amp_ia = t_ia + np.where(npx_sta_array[central_chan,t_ia:,spike] <=half_amp)[0][0]
        
        half_amp_ib = neg_pk_i + np.where(npx_sta_array[central_chan,neg_pk_i:,spike] >=half_amp)[0][0]
        
        half_width = (half_amp_ib - half_amp_ia) * 1000/samprate
        if half_width <0:
            if spike >=1:
                half_width = gt_spike_features[spike-1, 1]
        
        latency = (neg_pk_i - 60) * 1000/samprate
        
        symmetry = (neg_pk_i - t_ia) / float(0.001+(pos_pk_i - t_ib))
        if symmetry > 10:
            if spike >=1:
                symmetry = gt_spike_features[spike-1, 3]
        
        if pos_pk_v != 0:
            neg_pos_r = abs(neg_pk_v/pos_pk_v)
        else:
            neg_pos_r = abs(neg_pk_v/(pos_pk_v+0.1))
            
        gt_spike_features[spike, 0] = pk2pk
        gt_spike_features[spike, 1] = half_width
        gt_spike_features[spike, 2] = duration
        gt_spike_features[spike, 3] = symmetry
        gt_spike_features[spike, 4] = neg_pos_r
        gt_spike_features[spike, 5] = latency
        gt_spike_features[spike, 6] = abs(neg_pk_v)
        gt_spike_features[spike, 7] = pos_pk_v


#GT spike vector is the average of the extracted spike features of all ground-truth spikes. Here to exclude background activity that may have summated with a GT spike,
#I compute a z vector for each spike, which is how much its features differed from the average spike's features, for that cell.
# Good spikes are below a mean difference threshold, bad spikes are above and get excluded from further reliability analysis.
    gt_spike_vector = np.mean(gt_spike_features[:,0:8], axis = 0)
    gt_spike_vector_std = np.std(gt_spike_features[:,0:8], axis = 0)
    z_vector = np.empty((1,8), dtype = float)
#Get mean z score for each spike
    for spike in range(len(patch_spikes)):
        for feature in range(len(gt_spike_features[0,:])-1):
            z_vector[0,feature] = abs((gt_spike_features[spike, feature] - gt_spike_vector[feature])/gt_spike_vector_std[feature])
            gt_spike_features[spike, 8] = np.mean(z_vector)
#
    good_spikes = np.where(gt_spike_features[:,8]<=1.0)[0]
    bad_spikes = np.where(gt_spike_features[:,8]>1.0)[0]
#
    pcent_spikes[ cells_to_analyse.index(cell) ] = round(100*len(good_spikes)/float(len(patch_spikes)),1)
#
    rolling_plot = np.empty((8, len(good_spikes),2), dtype=float)
#
    feature_labels = ['Peak-peak\namplitude (uV)', 'Half-width (ms)', 'Duration (ms)', 'Duration\nSymmetry', 'Negative-positive\nratio', 'Latency\nto patch (ms)', 'Negative\nPeak (uV)', 'Positive\nPeak (uV)' ]
#
    fig, axs = plt.subplots(4,2, sharex = 'col', sharey=False, figsize=(15,10))
    axs = axs.ravel()
#Rolling plots show the rolling mean over 100 spike events +- the rolling STD of 100 events. Here I calculate rolling mean and stdev for each feature as well.
    for feature in range(8):
        rolling_plot[feature,:,0] = rolling_mean(gt_spike_features[good_spikes, feature], 100 )
        rolling_plot[feature,:,1] = rolling_std(gt_spike_features[good_spikes, feature], 100 )
        
        axs[feature].scatter(range(len(good_spikes)), gt_spike_features[good_spikes, feature ], s = 10, c =  '#011993', alpha = 0.3)
        axs[feature].plot( rolling_plot[feature,:,0], c = 'r' )
        axs[feature].fill_between(range(len(good_spikes)), rolling_plot[feature,:,0] - rolling_plot[feature,:,1], rolling_plot[feature,:,0] + rolling_plot[feature,:,1], color = 'orange'   )
        axs[feature].set_ylabel(feature_labels[feature], fontweight = 'bold', fontsize = 12)
        axs[feature].tick_params(axis='both',labelsize = 10)
    axs[1].set_ylim([0, 1])
    axs[2].set_ylim([0, 2])
    axs[5].set_ylim([-1, 1])
    plt.subplots_adjust(top = 0.95, bottom = 0.05, left = 0.15, right = 0.85)
    fig.suptitle(data_summary.loc[cell]['Cell'], fontweight = 'bold')
    #
    plt.savefig( output_dir + data_summary.loc[cell]['Cell'] + '.png', dpi = 300)
    plt.close()

    X = np.empty(len(good_spikes), dtype = int)
    for i in range(len(good_spikes)):
        X[i] = i
        X = X.reshape(-1, 1)
    linreg = LinearRegression()
#Here I calculate summary statistics for the running mean and standard deviation. I obtain a grand mean of every feature per cell, which is the mean of the running average,
#and a grand standard deviation which is the mean of the running standard deviation. This is done in order to attempt to control for the slow timescale change in features
#that we often observed, which we cannot exclude may be due to simultaneously recording the neuron with a patch electrode. I also calculate the coefficient of variation.
#Moreover, I also calculate the linear regression slope coefficient for every cell and feature, checking if there is an overall tendency for a particular feature
#to move in a particular direction over time. Features are taken from here and plotted in graphpad.    
    for feature in range(8):
        summary_stats[row, feature, 0] = np.nanmean(rolling_plot[feature, :, 0])
        summary_stats[row, feature, 1] = np.nanmean(rolling_plot[feature, :, 1])
        summary_stats[row, feature, 2] = abs(summary_stats[row, feature, 1] / summary_stats[row, feature, 0])
        Y = gt_spike_features[good_spikes, feature]
        linreg.fit(X, Y)
        summary_stats[row, feature, 3] =  linreg.coef_[0]

#
    parsed_spikes[row, 0] = 2.34375*npx_sta_array[central_chan,30:90, good_spikes]
    parsed_spikes[row, 1] = 2.34375*npx_sta_array[central_chan,30:90, bad_spikes]

#%%subplot for "good" spikes - Fig 6A
spike_splits, axspikes = plt.subplots(3, 3, sharex = False, sharey = False, figsize=(7.5,15))
axspikes = axspikes.ravel()
#
for cell in cells_to_analyse:
   # for spike in range(len(parsed_spikes[cells_to_analyse.index(cell), 1])):
   #     axspikes[cells_to_analyse.index(cell)].plot( parsed_spikes[cells_to_analyse.index(cell), 1][spike] - np.median(parsed_spikes[cells_to_analyse.index(cell), 1][spike] ), alpha = 0.1, color = 'orange')
    #
    for spike in range(len(parsed_spikes[cells_to_analyse.index(cell), 0])):
        axspikes[cells_to_analyse.index(cell)].plot(parsed_spikes[cells_to_analyse.index(cell), 0][spike] - np.median( parsed_spikes[cells_to_analyse.index(cell), 0][spike] ), alpha = 0.05, color = 'blue', lw= 0.7)
    axspikes[cells_to_analyse.index(cell)].set_title('%s'%(pcent_spikes[cells_to_analyse.index(cell)][0])+'%', fontweight = 'bold', loc='right', y = 0.05, color = 'blue' )
    axspikes[cells_to_analyse.index(cell)].plot(np.mean( parsed_spikes[cells_to_analyse.index(cell), 0], axis = 0 ), color = 'c' )
    axspikes[cells_to_analyse.index(cell)].set_xticks([0,15, 30, 45, 60])
    axspikes[cells_to_analyse.index(cell)].set_xticklabels([-1, -0.5, 0, 0.5, 1])
    axspikes[cells_to_analyse.index(cell)].set_title(data_summary.loc[cell]['Cell'], fontweight = 'bold', color = 'k' )
plt.subplots_adjust(hspace = 0.3, wspace = 0.3)
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 5_materials/good_spikes.png', dpi = 300)
plt.close()
#%%subplot for "bad" spikes - Supplement to Fig 6-1A
spike_splits, axspikes = plt.subplots(3, 3, sharex = False, sharey = False, figsize=(7.5,15))
axspikes = axspikes.ravel()

for cell in cells_to_analyse:
   # for spike in range(len(parsed_spikes[cells_to_analyse.index(cell), 1])):
   #     axspikes[cells_to_analyse.index(cell)].plot( parsed_spikes[cells_to_analyse.index(cell), 1][spike] - np.median(parsed_spikes[cells_to_analyse.index(cell), 1][spike] ), alpha = 0.1, color = 'orange')
    #
    for spike in range(len(parsed_spikes[cells_to_analyse.index(cell), 1])):
        axspikes[cells_to_analyse.index(cell)].plot(parsed_spikes[cells_to_analyse.index(cell), 1][spike] - np.median( parsed_spikes[cells_to_analyse.index(cell), 1][spike] ), alpha = 0.05, color = 'orange', lw= 0.7)
    axspikes[cells_to_analyse.index(cell)].set_xticks([0,15, 30, 45, 60])
    axspikes[cells_to_analyse.index(cell)].set_xticklabels([-1, -0.5, 0, 0.5, 1])
    axspikes[cells_to_analyse.index(cell)].set_title(data_summary.loc[cell]['Cell'], fontweight = 'bold', color = 'k' )
plt.subplots_adjust(hspace = 0.3, wspace = 0.3)
#
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 5_materials/bad_spikes.png', dpi = 300)
plt.close()