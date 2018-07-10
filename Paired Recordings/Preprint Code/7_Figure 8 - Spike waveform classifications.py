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
#%%Functions
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a
    
def butter_lowpass_filter(data, cutoff, fs, order=5):
        b, a = butter_lowpass(cutoff, fs, order=order)
        y = lfilter(b, a, data)
        return y

def classify_waveform(recording, Tm = 3):
#   
    #Tm = 1.5
    #chan = random.randint(central_chan+30,central_chan+50)
    #recording = npx_sta_mean[chan,:]
    
        
    filtered =  butter_lowpass_filter(recording, 3000, 30000)
    filtered_diff = np.diff(filtered)
    #T = np.std( filtered[0:40])
    T = np.empty((1,120),dtype = float)
    T[0, 0:60] = np.std( filtered[0:40])
    T[0, 60:] = 10*np.std( filtered[0:40])
    Tneg = np.empty((1,120),dtype = float)
    Tneg[0, 0:65] = -np.std( filtered[0:40])*Tm
    Tneg[0, 65:] = -np.std( filtered[0:40])*30
    #%
    filtered_diff_possub = [filtered_diff[i] - T[0,i] for i in range(len(filtered_diff)) ]
    #going_up = np.where(np.diff(np.sign(filtered_diff - T))>0)
    going_up = np.where(np.diff(np.sign(filtered_diff_possub))>0)
    #going_down = np.where(np.diff(np.sign(filtered_diff + 10*T))<0)
    filtered_diff_negsub = [filtered_diff[i] - Tneg[0,i] for i in range(len(filtered_diff)) ]
    going_down = np.where(np.diff(np.sign(filtered_diff_negsub))<0)
    
    #
    allcrosses = sorted(np.concatenate((going_up[0], going_down[0])))
    for i in range(len(allcrosses)):
        if allcrosses[i] in going_up[0]:
            allcrosses[i] = np.sign(allcrosses[i])
        if allcrosses[i] in going_down[0]:
            allcrosses[i] = np.sign(-allcrosses[i])
    
    if len(allcrosses) == 1:
        if allcrosses[0] == 1:
            Tm = 0.5
            Tneg[0, :] = -np.std( filtered[0:40])*Tm
            filtered_diff_negsub = [filtered_diff[i] - Tneg[0,i] for i in range(len(filtered_diff)) ]
            going_down = np.argmin(filtered_diff_negsub)-10 + np.where(np.diff(np.sign( filtered_diff_negsub[ np.argmin(filtered_diff_negsub)-10 :np.argmin(filtered_diff_negsub)+10 ] ))<0)
            allcrosses = None
            allcrosses = sorted(np.concatenate((going_up[0], going_down[0])))
            for i in range(len(allcrosses)):
                if allcrosses[i] in going_up[0]:
                    allcrosses[i] = np.sign(allcrosses[i])
                if allcrosses[i] in going_down[0]:
                    allcrosses[i] = np.sign(-allcrosses[i])
                    
        elif allcrosses[0] == -1:
            Tm = 0.1
            T[0, 60:] = 10*np.std( filtered[0:40]) * Tm
            filtered_diff_possub = [filtered_diff[i] - T[0,i] for i in range(len(filtered_diff)) ]
            going_up = np.argmin(filtered_diff_possub[40:60])  + np.where(np.diff(np.sign( filtered_diff_possub[ np.argmin(filtered_diff_possub[40:60]) :  ] ))>0)
            allcrosses = None
            allcrosses = sorted(np.concatenate((going_up[0], going_down[0])))
            for i in range(len(allcrosses)):
                if allcrosses[i] in going_up[0]:
                    allcrosses[i] = np.sign(allcrosses[i])
                if allcrosses[i] in going_down[0]:
                    allcrosses[i] = np.sign(-allcrosses[i])
            
    
    string = ''.join(str(e) for e in allcrosses)
    return string
#%%Dataframes and structures

listcells = os.listdir('E:/code/for analysis')
suffix = '_aligned_analysis'

data_summary = pandas.read_excel('C:/Users/Andre Marques-Smith/Dropbox/Patch-Clamp Summary 2.0.xlsx')


cells_to_analyse = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 50].tolist()
#%%INs
cells_to_analyse = data_summary.index[data_summary['Cell Type'] == 'IN'].tolist()
#%%example for waveform traces
cells_to_analyse = [19]
#%%
#Data structures for waveform analysis
num_rows = 25
wave_classes = np.empty((num_rows * 4,len(cells_to_analyse)), dtype = object)    
electrodes =  np.empty((num_rows * 4,len(cells_to_analyse)), dtype = object)
all_npx_stas = np.empty( (384,121, len(cells_to_analyse)), dtype = float)

submap = np.full((num_rows * 2,4,len(cells_to_analyse)),-100, dtype = int)
#


for cell in cells_to_analyse:
    paths = defaultdict(list)
    cell_idx = listcells.index(data_summary.loc[cell]['Cell']+suffix)

    aligned_directory = 'E:/code/for analysis/'+listcells[cell_idx]
    
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
    #npx_voltage = np.memmap( paths['npx_v'][0], mode = 'c', dtype=np.int16 )
    #npx_voltage = npx_voltage.reshape((expt_meta['npx'][0][0], expt_meta['npx'][0][1]))
    #
    npx_sta_array = np.load(paths['npx_sta_array'][0])#, dtype = 'int16')
    npx_sta_mean = np.load(paths['npx_sta_mean'][0])#, dtype = 'float64')
    npx_chan_peaks = np.load(paths['npx_chan_peaks'][0])#, dtype = 'float64')
    
    central_chan = int(npx_chan_peaks[0,0])
    
    
    r_cell = np.where(npx_mapcsv == central_chan)[0][0]
#    c_cell = np.where(npx_mapcsv == central_chan)[1][0]
    
    
    submap[:,:,cells_to_analyse.index(cell)] = npx_mapcsv[r_cell-num_rows:r_cell+num_rows]
    
    #electrodes[:, cells_to_analyse.index(cell)] = submap[submap > 0 ]
    electrodes[:, cells_to_analyse.index(cell)] = submap[submap[:,:,cells_to_analyse.index(cell)] > 0 ,cells_to_analyse.index(cell)]
    
    
    all_npx_stas[:,:,cells_to_analyse.index(cell)] = npx_sta_mean / abs(np.max(npx_sta_mean[:,:])-np.min(npx_sta_mean[:,:]))
 #  
    for chan in range(len(electrodes)):
        wave_classes[chan ,cells_to_analyse.index(cell)] = classify_waveform(npx_sta_mean[electrodes[:,cells_to_analyse.index(cell)][chan],:], Tm=1.5)
        #%%Waveform classification function classifies some waveforms redundantly (-111 is the same as -11). This is corrected below
    for cell in range(len(cells_to_analyse)):
        for chan in range(len(electrodes)):
            if wave_classes[chan ,cell] == '':
                wave_classes[chan ,cell] = 'None'
            if wave_classes[chan ,cell] == '-111':
                wave_classes[chan ,cell] = '-11'
            if wave_classes[chan ,cell] == '1-111':
                wave_classes[chan ,cell] = '1-11'
                


#%%Here python plots the frequency of each type of waveform classification it found. Because my function is *very* basic indeed, it will classify small wingnumbers of waveforms
#weirdly, but will get the majority and easy cases right. Here we will simply focus on the most frequent waveform classifications, requiring a minimum proportion of
#waveforms to be classified in order to further analyse it. This is determined by examining the barplot visually and setting the variable min_prop.
wave_types = np.empty((len(np.unique(wave_classes[:,:])),2),dtype = object)
wave_types[:,0] = np.unique(wave_classes[:,:])
#
for i in range(len(wave_types)):
    wave_types[i,1] = len( np.where(wave_classes[:,:] == wave_types[i,0])[0] )
    #wave_types[i,1] /= float(np.sum(wave_types[:,1]))
wave_types[:,1] /= float(np.sum(wave_types[:,1]))
#np.where(wave_classes[:] == wave_types[i,0])[0]
#
plt.bar( range(len(wave_types[:,0])), wave_types[:,1] )
plt.xticks(range(len(wave_types[:,0])), wave_types[:,0], rotation='vertical' )
#%%
min_prop = 0.04
waves_to_count = wave_types[np.where(wave_types[0:-1,1] > min_prop) ,0][0]
#%%
waves_counts_by_cell = np.empty((len(cells_to_analyse), len(waves_to_count)+1), dtype = int)

for cell in cells_to_analyse:
    for wave in range(len(waves_to_count)):
        waves_counts_by_cell[cells_to_analyse.index(cell), wave] = len( np.where(wave_classes[:,cells_to_analyse.index(cell)] == waves_to_count[wave])[0] )
    waves_counts_by_cell[cells_to_analyse.index(cell) ,len(waves_to_count)] = np.sum(waves_counts_by_cell[ cells_to_analyse.index(cell)], axis = 0)

#%%
histos = np.empty((len(electrodes), len(waves_to_count)), dtype = float)
histos[:,0] = range(len(electrodes))

for wave in range(len(waves_to_count)):
    for chan in range(len(electrodes)):
        histos[chan,wave] = len(np.where( wave_classes[chan,:] == waves_to_count[wave] )[0])
#
     
#
total = np.sum(histos)

for row in range(len(histos)):
    for col in range (len(waves_to_count)):
        histos[row, col] /= total
        
print total
#
probe_histos =   np.zeros((len(electrodes),4,len(waves_to_count)), dtype = float)  

for cell in cells_to_analyse:
    for chan in range(len(electrodes)):
        for wave in range(len(waves_to_count)):
            r = np.where(submap[:,:, cells_to_analyse.index(cell)] == electrodes[chan, cells_to_analyse.index(cell)])[0][0]
            c = np.where(submap[:,:, cells_to_analyse.index(cell)] == electrodes[chan, cells_to_analyse.index(cell)])[1][0]
            probe_histos[r, c, wave] = histos[chan, wave]
#

template_map = np.full((len(electrodes)/2,4), -1, dtype = int)

e = len(electrodes)-1
for r in list(reversed(range(len(electrodes)/2))):
    for c in range(4):
        if template_map[r, c] != submap[r,c,-1]:
            template_map[r,c] = e
            e -= 1
template_electrodes = list(reversed(range(len(electrodes))))
probe_histos =   np.zeros((len(electrodes)/2,4,len(waves_to_count)), dtype = float)  

for chan in range(len(electrodes)):
    for wave in range(len(waves_to_count)):
        r = np.where(template_map[:,:] == chan)[0][0]
        c = np.where(template_map[:,:] == chan)[1][0]
        probe_histos[r, c, wave] = histos[chan, wave]
#%%wave_stats
wave_stats = np.empty((len(cells_to_analyse),5), dtype = float)
waves = ['-11', '1-1', '1-11', 'None']

for cell in range(len(cells_to_analyse)):
    for wave in range(len(waves)):
        wave_stats[cell, wave] = len(np.where(wave_classes[:,cell] == waves[wave])[0])
    wave_stats[cell, 4] = len([i for i in wave_classes[:,cell] if i not in waves])
#%%

#%%Plot frequency distribution of waveforms of a particular type per channel of the probe relative to the soma position as a color-code, for the 3 key types of waveform.
histmax = round(np.max(probe_histos),3)
colorbin = round(histmax/4,3)
fig, ax = plt.subplots(1,len(waves_to_count))
#
cm = 'afmhot'

for wave in range(len(waves_to_count)):
    wave_heatmap = ax[wave].imshow(probe_histos[:,:,wave], cmap = cm,interpolation='none', aspect = 1, vmin = 0, vmax = histmax    )
    if waves_to_count[wave] == '-11':
        axtitle = 'Biphasic\nnegative'
    elif waves_to_count[wave] == '1-1':
        axtitle = 'Biphasic\npositive'
    elif waves_to_count[wave] == '1-11':
        axtitle = '\nTriphasic'
    else:
        axtitle = 'Other'
    ax[wave].set_title(axtitle, fontweight = 'bold', y = 1.05)
    ax[wave].set_xticks([])
    ax[wave].set_yticks([0, len(template_map)/2, len(template_map)-1])
    if wave == 0:
        ax[wave].set_yticklabels(['+480 um', 'Soma', '-480 um'  ])
        #ax[wave].set_yticklabels(['-480 um', 'Soma', '+480 um'  ])
    else:
        ax[wave].set_yticklabels(['', '->', ''], fontweight = 'bold')
plt.subplots_adjust(top = 0.86, bottom = 0.11, left = 0.25, right = 0.75, wspace = 0.0)
cbar_ax = fig.add_axes([0.35, 0.05, 0.3, 0.02]) #left bottom width height
cbar = fig.colorbar(wave_heatmap, cax = cbar_ax, orientation = 'horizontal', ticks = [0, histmax/2.0, histmax] )
cbar.ax.set_xticklabels(['0', '%s' %(100*round(histmax/2.0,3))+'%', '%s' %( 100*round(histmax,3) )+'%'])
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 7_materials/' + 'Fig 7C_IN_wavetype_distribution.png' , dpi = 1200)
plt.close()
#%%plt an example cell on probe
#fig, ax = plt.subplots(50,4, sharex=True, sharey = 'row', figsize=(3,11))
#ax = ax.ravel()
for cell in cells_to_analyse:
    central_chan = data_summary.loc[cell]['chan_highest']
    fig, ax = plt.subplots(50,4, sharex=True, sharey = 'row', figsize=(3,11))
    for chan in range(len(electrodes)):
        r = np.where(submap[:,:,cells_to_analyse.index(cell)] == electrodes[chan][cells_to_analyse.index(cell)])[0][0]
        c = np.where(submap[:,:,cells_to_analyse.index(cell)] == electrodes[chan][cells_to_analyse.index(cell)])[1][0]
        ax[r, c].axis('off')
        if electrodes[chan][cells_to_analyse.index(cell)] == central_chan:
            ax[r,c].set_title('%s'%(electrodes[chan][cells_to_analyse.index(cell)]), fontsize = 5, y = 0.35, alpha = 1, loc='right', fontweight='bold')
        else:
            ax[r,c].set_title('%s'%(electrodes[chan][cells_to_analyse.index(cell)]), fontsize = 5, y = 0.35, alpha = 0.4, loc='right')
        if wave_classes[chan, cells_to_analyse.index(cell)] == '1-1':
            ax[r, c].plot(cells_to_analyse.index(cell)*0.1+ all_npx_stas[electrodes[chan][cells_to_analyse.index(cell)],:,cells_to_analyse.index(cell)], c = '#FFBB00', lw = 0.7)
        elif wave_classes[chan, cells_to_analyse.index(cell)] == '-11':
            ax[r, c].plot(cells_to_analyse.index(cell)*0.1+all_npx_stas[electrodes[chan][cells_to_analyse.index(cell)],:,cells_to_analyse.index(cell)], c = '#FB6542', lw = 0.7)
        elif wave_classes[chan, cells_to_analyse.index(cell)] == '1-11':
            ax[r, c].plot(cells_to_analyse.index(cell)*0.1+all_npx_stas[electrodes[chan][cells_to_analyse.index(cell)],:,cells_to_analyse.index(cell)], c = '#3F681C', lw = 0.7)
        else:
            ax[r, c].plot(cells_to_analyse.index(cell)*0.1+all_npx_stas[electrodes[chan][cells_to_analyse.index(cell)],:,cells_to_analyse.index(cell)], c = '#375E97', alpha = 0.2, lw = 0.7)
    for r in range(50):
        for c in range(4):
            ax[r, c].axis('off')
#
    plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 7_materials/' + 'Fig 7_'+ data_summary.loc[cell]['Cell'] +'_PC_example_waves.png' , dpi = 1200)
    plt.close()
#%%Fig 7D: plot a column of spike waveforms aligned by patch spike peak time
for cell in cells_to_analyse:
    central_chan = data_summary.loc[cell]['chan_highest']
    target = np.where(npx_mapcsv == central_chan )[1]
    col = []
    
    for electrode in electrodes[:, cells_to_analyse.index(cell)]:
        r = np.where(npx_mapcsv == electrode)[0]
        c = np.where(npx_mapcsv == electrode)[1]
        if  c ==target:
            col.append( np.where(electrodes[:, cells_to_analyse.index(cell)] == electrode)[0] )
#
    soma_t = np.argmin(npx_sta_mean[central_chan])
    fig, ax = plt.subplots(len(col),1, sharex = True, sharey='row', figsize=( 1.22, 5.32))
    for chan in range(len(col)):
        if wave_classes[col[chan], cells_to_analyse.index(cell)] == '1-11':
            ax[chan].plot( all_npx_stas[ electrodes[col[chan], cells_to_analyse.index(cell)][0],:,cells_to_analyse.index(cell)], c = '#3F681C', lw = 0.7)
            #lag = (np.argmin(npx_sta_mean[ electrodes[col[chan], cells_to_analyse.index(cell)][0] ] ) - soma_t) * 33
            if electrodes[col[chan], cells_to_analyse.index(cell)] == central_chan:
                ax[chan].set_title('ch. %s'%electrodes[col[chan], cells_to_analyse.index(cell)][0], loc = 'left', fontsize = 4, y = 0.7, alpha = 0.3, fontweight = 'bold')
            else:
                ax[chan].set_title('ch. %s'%electrodes[col[chan], cells_to_analyse.index(cell)][0], loc = 'left', fontsize = 4, y = 0.7, alpha = 0.3)
        elif wave_classes[col[chan], cells_to_analyse.index(cell)] == '1-1':
            ax[chan].plot( all_npx_stas[electrodes[col[chan], cells_to_analyse.index(cell)][0] ,:,cells_to_analyse.index(cell)], c = '#FFBB00', lw = 0.7)
            #lag = (np.argmin(npx_sta_mean[ electrodes[col[chan], cells_to_analyse.index(cell)][0] ] ) - soma_t) * 33
            if electrodes[col[chan], cells_to_analyse.index(cell)] == central_chan:
                ax[chan].set_title('ch. %s'%electrodes[col[chan], cells_to_analyse.index(cell)][0], loc = 'left', fontsize = 4, y = 0.7, alpha = 0.3, fontweight = 'bold')
            else:
                ax[chan].set_title('ch. %s'%electrodes[col[chan], cells_to_analyse.index(cell)][0], loc = 'left', fontsize = 4, y = 0.7, alpha = 0.3)
        elif wave_classes[col[chan], cells_to_analyse.index(cell)] == '-11':
            ax[chan].plot( all_npx_stas[electrodes[col[chan], cells_to_analyse.index(cell)][0] ,:,cells_to_analyse.index(cell)], c = '#FB6542', lw = 0.7)
            #lag = (np.argmin(npx_sta_mean[ electrodes[col[chan], cells_to_analyse.index(cell)][0] ] ) - soma_t) * 33
            if electrodes[col[chan], cells_to_analyse.index(cell)] == central_chan:
                ax[chan].set_title('ch. %s'%electrodes[col[chan], cells_to_analyse.index(cell)][0], loc = 'left', fontsize = 4, y = 0.7, alpha = 0.3, fontweight = 'bold')
            else:
                ax[chan].set_title('ch. %s'%electrodes[col[chan], cells_to_analyse.index(cell)][0], loc = 'left', fontsize = 4, y = 0.7, alpha = 0.3)
        else:
            ax[chan].plot( all_npx_stas[electrodes[col[chan], cells_to_analyse.index(cell)][0] ,:,cells_to_analyse.index(cell)], c = '#375E97', alpha = 0.2, lw = 0.7)
            #lag = (np.argmin(npx_sta_mean[ electrodes[col[chan], cells_to_analyse.index(cell)][0] ] ) - soma_t) * 33
            if electrodes[col[chan], cells_to_analyse.index(cell)] == central_chan:
                ax[chan].set_title('ch. %s'%electrodes[col[chan], cells_to_analyse.index(cell)][0], loc = 'left', fontsize = 4, y = 0.7, alpha = 0.3, fontweight = 'bold')
            else:
                ax[chan].set_title('ch. %s'%electrodes[col[chan], cells_to_analyse.index(cell)][0], loc = 'left', fontsize = 4, y = 0.7, alpha = 0.3)
        ax[chan].axis('off')
            
    plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 7_materials/' + 'Fig 7_'+ data_summary.loc[cell]['Cell'] +'_PC_column_waves.png' , dpi = 1200)
    plt.close()
#%%