#Imports
import os
import pandas
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

#%%Functions
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
    
    return [(int(i[:2], 16)/255.0, int(i[2:4], 16)/255.0, int(i[4:], 16)/255.0) for i in colors]

def extracellular_spike_features_adaptive(cell, samprate = 30000.0, delta = 5):
  #
    #cell = all_npx_stas[cells_to_analyse[0]]
    #samprate = 30000.0
    #delta = 5
  #  
    adaptive = 1
    baseline = int(len(cell)/3)
    
    
    t = np.median( cell[:baseline] ) - np.std(cell[:baseline])*adaptive
    
    #while(cell[t:t+5] > t )
    
    t_ia = 40+np.where(cell[40:] <=t)[0][0]
    
    while (cell[t_ia:t_ia+delta] > t).any():

        adaptive += 0.5
        t = np.median( cell[:baseline] ) - np.std(cell[:baseline])*adaptive
        t_ia = np.where(cell <=t)[0][0]
        
    #
    tro_i = np.argmin( cell ) #brought one line up
    t_ib = tro_i + np.where( cell[tro_i:]>=t)[0][0] #replaced t_ia with tro_i here
    
    
    half_amp = (np.min(cell)+t)/2.0
    half_a = np.where(cell<=half_amp)[0][0] - 1
    half_b = half_a + 1 + np.where(cell[half_a + 1:]>=half_amp)[0][0]
    
    symmetry_r = (tro_i - t_ia)/float((t_ib - tro_i))
    
    half_width = (half_b - half_a) * 1000/samprate
    
    tro_pk_t = (np.argmax(cell[t_ia:]) - np.argmin(cell[t_ia:])) * 1000/samprate
    tro_pk_r = abs(np.min(cell[t_ia:])/float(np.max(cell[t_ia:])))
    
    pk2pkamp = cell[np.argmax(cell)] - cell[np.argmin(cell)]
    
    plt.plot(cell, color = '#0433FF', lw=2)
    #plt.scatter(t_ia, cell[t_ia], s = 80, c = 'orange')
    #plt.scatter(t_ib, cell[t_ib], s = 80, c = 'orange')
    plt.hlines(t,0,len(cell),'r', linestyles='--')
    plt.hlines(half_amp,half_a,half_b,'r', lw=2)
    plt.vlines(tro_i, np.min(cell), t, linestyles='--')
    plt.vlines(np.argmax(cell[t_ia:]) + t_ia, t, np.max(cell), linestyles='--', color = 'magenta')
    plt.hlines(np.min(cell)-2, np.argmin(cell), np.argmax(cell), lw = 2, color = 'cyan' )
  #
    return pk2pkamp, half_width, symmetry_r, tro_pk_t, tro_pk_r, t, half_amp, half_a, half_b, tro_i, t_ia
#%%Dataframes and structures
listcells = os.listdir('E:/code/for analysis')
suffix = '_aligned_analysis'

data_summary = pandas.read_excel('C:/Users/Andre Marques-Smith/Dropbox/Patch-Clamp Summary 2.0.xlsx')

all_npx_stas = np.zeros((len(data_summary),121),dtype = float)

for cell in range(len(data_summary['Cell'])):
    aligned_directory = 'E:/code/for analysis/'+data_summary.loc[cell]['Cell']+suffix
    for file in os.listdir(aligned_directory):
        if file.endswith('sta_np_by_channel.npy'):
            npx_sta = np.load(aligned_directory + '/' + file)
            all_npx_stas[cell] = npx_sta[data_summary.loc[cell]['chan_highest']]
        else:
            pass

#%%#Admission criteria for further analysis: pk-pk > 10 uV; # of patch spikes recorded > 200; canonical waveform ie single negative peak leading the positive peak;
#or if positive peak leads negative peak, negative has greater amplitude
cells_above_10uV = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 10].tolist()
excluded_cells = [2, 5, 9, 12, 21, 23] #updated 02.05.2018; c16 highest channel changed after correcting lag error
#excluded_cells = [2, 6, 9, 10, 12, 21, 23]
cells_to_analyse = [ cells_above_10uV[i] for i in range(len(cells_above_10uV)) if i not in excluded_cells]

#%%Supplement to Fig 4-1A: average spike waveforms for excluded cells
fig, ax = plt.subplots(3,2, sharex = 'col', sharey = False)
ax=ax.ravel()

for cell in range(len(excluded_cells)):
    ax[cell].plot(all_npx_stas[cells_above_10uV[excluded_cells[cell]]], color = '#0432FE')
    ax[cell].set_xticks([0, 30, 60, 90, 120]) 
    ax[cell].set_xticklabels(['-2', '-1', '0', '1', '2'])
    ax[cell].set_yticks([ 0])
    ax[cell].set_title(data_summary.loc[cells_above_10uV[excluded_cells[cell]]]['Cell'], y = 1.0, fontweight = 'bold')
    ax[cell].vlines(15,-5,5)
    #plt.ylim([-20,20])
    plt.xlim([0,120])
    plt.subplots_adjust(top = 0.94, bottom = 0.11, left = 0.2, right = 0.6, hspace = 0.3)

#ax[7].axis('off')
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Supp to Fig 3_exclusions.png', dpi = 1200)
#%%Supplement to Fig 4-1B: average spike waveforms for further-analysed cells
fig, ax = plt.subplots(7,3, sharex = 'col', sharey = False)
ax=ax.ravel()

for cell in range(len(cells_to_analyse)):
    ax[cell].plot(all_npx_stas[cells_to_analyse[cell]], color = '#0432FE')
    ax[cell].set_xticks([0, 30, 60, 90, 120]) 
    ax[cell].set_xticklabels(['-2', '-1', '0', '1', '2'])
    ax[cell].set_yticks([ 0])
    ax[cell].set_title(data_summary.loc[cells_to_analyse[cell]]['Cell'], y = 1.0, fontweight = 'bold')
    ax[cell].vlines(15,-5,5)
    #plt.ylim([-20,20])
    plt.xlim([0,120])
    plt.subplots_adjust(top = 0.94, bottom = 0.11, left = 0.2, right = 0.8, hspace = 0.3)
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Supp to Fig 3_inclusions.png', dpi = 1200)

#%%#Extract and organise PSTA spike features
#spike feature column indexes are in order pk2pk, HW, symmetry, duration, ratio and threshold

STA_spike_features = np.empty((len(cells_to_analyse), 12), dtype = float)

for cell in range(len(cells_to_analyse)):
    STA_spike_features[cell,0:11] = extracellular_spike_features_adaptive(all_npx_stas[cells_to_analyse[cell]])
    
for cell in range(len(cells_to_analyse)):
    if STA_spike_features[cell, 1] <0.3:
        if STA_spike_features[cell, 3]<0.4:
            STA_spike_features[cell, 11] = 0
        else:
            STA_spike_features[cell, 11] = 1            
    else:
        STA_spike_features[cell, 11] = 1
#%%Supplement to Fig 4-2: illustration of extraction of spike features for each of 21 further analysed cells.
fig, axs = plt.subplots(7,3, sharex='col', sharey=False)
axs = axs.ravel()

for cell in range(len(cells_to_analyse)):
    if STA_spike_features[cell, 11] ==0:
        axs[cell].plot(all_npx_stas[cells_to_analyse[cell], int(STA_spike_features[cell, 9]-45):int(STA_spike_features[cell, 9]+45)], color = 'cyan', lw = 0.7)
    elif STA_spike_features[cell, 11] ==1:
        axs[cell].plot(all_npx_stas[cells_to_analyse[cell], int(STA_spike_features[cell, 9]-45):int(STA_spike_features[cell, 9]+45)], color = 'magenta', lw = 0.7)    
    axs[cell].set_xticks([0, 30, 60, 90, 120]) 
    axs[cell].set_xticklabels(['-2', '-1', '0', '1', '2'])
    #axs[cell].set_yticks([-20, 0, 20])
    
    axs[cell].set_title(data_summary.loc[cells_to_analyse[cell]]['Cell'], y = 1, fontweight = 'bold')

    axs[cell].hlines(STA_spike_features[cell, 5],0,90, 'r', linestyles='--', lw = 0.8  )           #threshold
    axs[cell].hlines(STA_spike_features[cell, 6], 45-STA_spike_features[cell, 9] + STA_spike_features[cell, 7] , 45 - STA_spike_features[cell, 9]+ STA_spike_features[cell, 8] , 'blue', lw = 0.8  )      #halfwidth
    axs[cell].vlines(45 , np.min(all_npx_stas[cells_to_analyse[cell]]), STA_spike_features[cell, 5], color='k', lw = 0.8) #peak
    pki = 45 - int(STA_spike_features[cell, 9]) + int(STA_spike_features[cell, 10]+np.argmax(all_npx_stas[cells_to_analyse[cell], int(STA_spike_features[cell, 10]):]))
    pkv = np.max(all_npx_stas[cells_to_analyse[cell], int(STA_spike_features[cell, 9]): ])
    axs[cell].vlines(pki, STA_spike_features[cell, 5], pkv , color = 'k', lw = 0.8 )       #positive peak
    axs[cell].hlines(np.min(all_npx_stas[cells_to_analyse[cell]]) -2, 45, pki, color = 'blue'   )  #pk2pk dur
    #axs[cell].set_xlim([STA_spike_features[cell, 9]-30,STA_spike_features[cell, 9]+30])
    axs[cell].vlines(0, STA_spike_features[cell, 5]-5, STA_spike_features[cell, 5]+5, color='grey')
    
    plt.subplots_adjust(hspace = 0.45, wspace = 0.45, top = 0.94, bottom = 0.06, right = 0.9, left = 0.1)
    axs[cell].axis('off')
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'All_spike_feature_measures.png', dpi = 1200)
#%%Figure 4C - scatterplot of spike half-width vs trough to peak duration to cluster putative interneurons vs pyramidal cells
fig = plt.figure(figsize=cm2inch(10,10))
ax1 = fig.add_subplot(111)
cm = plt.get_cmap('cool',2)
cell_types = ax1.scatter(STA_spike_features[:,1],STA_spike_features[:,3], c = STA_spike_features[:,11], vmin = 0, vmax = 1, cmap = cm, s = 45, edgecolor='k', lw=0.5)

ax1.set_xlim([0, 0.5])
ax1.set_ylim([0, 1.2])
plt.xlabel('Spike Half-Width (ms)', labelpad = 5, fontsize = 11, fontweight = 'bold')
plt.ylabel('Trough to Peak Duration (ms)', labelpad = 5, fontsize = 11, fontweight = 'bold')
plt.subplots_adjust(top = 0.90, bottom = 0.2, left = 0.20, right = 0.95)
#plt.tight_layout
plt.grid('on', alpha = 0.4)

#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Fig 3C_scatter_cell_types.png', dpi = 1200)
#%%Example average extracelluar spike traces in Fig 4C

example_cells = [4, 20, 3, 6]

fig, ax = plt.subplots(4,1, figsize=(0.25,6), sharex = True, sharey = False)
#plt.xlim(30,91)
for cell in range(2,4):
    ax[cell].plot( all_npx_stas[cells_to_analyse[example_cells[cell]]], c = '#00FDFD', lw = 0.8)
    ax[cell].axis('off')
    ax[cell].set_title(data_summary.loc[cells_to_analyse[example_cells[cell]]]['Cell'], loc = 'left', fontweight = 'bold', y = 0.70)
    ax[cell].vlines(115, -5, -10)
    ax[cell].set_ylim([-10, 9])
for cell in range(2):
    ax[cell].plot( all_npx_stas[cells_to_analyse[example_cells[cell]]], c = '#F909F9', lw = 0.8)
    ax[cell].axis('off')
    ax[cell].set_title(data_summary.loc[cells_to_analyse[example_cells[cell]]]['Cell'], loc = 'left', fontweight = 'bold', y = 0.70)
    ax[cell].vlines(115, -155, -105)
    ax[cell].set_ylim([-155, 95])
plt.subplots_adjust(hspace = 0.15, top = 0.96, bottom = 0.04, left = 0.25, right = 0.75)
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Fig 3C_inset_traces.png', dpi = 1200)
#%%group non-analysed cells
remaining_cells = [i for i in range(len(data_summary)) if i not in cells_to_analyse]

#%%Plot extracellular peak-peak amplitude vs distance - figure 4A

fig = plt.figure(figsize=cm2inch(10, 10))

ax1 = fig.add_subplot(111)
good_cells_plot = ax1.scatter(data_summary['Distance'][cells_to_analyse] ,data_summary['JTA Peak-Peak Amplitude'][cells_to_analyse], c = '#F99120', edgecolor='k', lw=0.5, s = 45)
bad_cells_plot = ax1.scatter(data_summary['Distance'][remaining_cells] ,data_summary['JTA Peak-Peak Amplitude'][remaining_cells], c = '#807F7F', edgecolor='k', lw=0.5, s = 45, alpha = 0.5)
ax1.grid(b=True, alpha = 0.4)
#box = ax1.get_position()
#ax1.set_position([box.x0*1.6, box.y0*1.2 + box.height * 0.25, box.width*0.95, box.height * 0.8])
#plt.legend((good_cells_plot, bad_cells_plot),('Detectable EAP waveform', 'Non-detectable/non-canonical\nEAP waveform'), fontsize = 9, loc='lower center', bbox_to_anchor=(0.45, -0.50), frameon=False)

#col = plt.colorbar(depth_plot ,orientation = 'horizontal', ticks = np.arange(300,2100,300), shrink = 0.75, pad = 0.17)
#col.set_label('Cortical Depth (um)', fontweight = 'bold')
#col.ax.tick_params(labelsize=9)   

plt.xticks(np.arange(0, np.max(data_summary['Distance'])+1, 20.0  ))
#plt.subplots_adjust(top = 0.98, bottom = 0.2, left = 0.20, right = 0.96)
plt.tick_params(labelsize=10)
plt.xlabel('Distance from Channel (um)', fontsize=11,labelpad = 5, fontweight = 'bold')
plt.ylabel('Peak-Peak Amplitude (uV)', fontsize=11, labelpad = 5, fontweight = 'bold')
plt.subplots_adjust(top = 0.90, bottom = 0.2, left = 0.20, right = 0.95)
plt.show()
#%%
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Fig 3A_scatter_depth_version 3.png', dpi = 1200)

#%%PLot Inset of recording depth vs amplitude scatter for Fig 4B
fig = plt.figure(figsize=cm2inch(3.4, 3.4))

ax1 = fig.add_subplot(111)
good_cells_plot = ax1.scatter(data_summary['Depth'][cells_to_analyse] ,data_summary['JTA Peak-Peak Amplitude'][cells_to_analyse], c = '#F99120', edgecolor='k', lw=0.2, s = 10)
bad_cells_plot = ax1.scatter(data_summary['Depth'][remaining_cells] ,data_summary['JTA Peak-Peak Amplitude'][remaining_cells], c = '#807F7F', edgecolor='k', lw=0.2, s = 10, alpha = 0.4)
#ax1.grid(b=True, alpha = 0.4)

#col = plt.colorbar(depth_plot ,orientation = 'horizontal', ticks = np.arange(300,2100,300), shrink = 0.75, pad = 0.17)
#col.set_label('Cortical Depth (um)', fontweight = 'bold')
#col.ax.tick_params(labelsize=9)   

plt.xticks(np.arange(0, np.max(data_summary['Depth'])+1, 400.0  ))
plt.subplots_adjust(top = 0.98, bottom = 0.25, left = 0.25, right = 0.98)
plt.tick_params(labelsize=6)
#plt.xlabel('Depth (um)', fontsize=7,labelpad = 2, fontweight = 'bold')
plt.ylabel('Amplitude (uV)', fontsize=7, labelpad = 2, fontweight = 'bold')
plt.yticks(rotation=90)
plt.show()

#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Fig 3_scatter_depth_inset.png', dpi = 1200)
#%%Organise data for plotting in graphpad (Fig 4E); cortical depth vs cell type vs peak-peak amplitude
subframe = data_summary.loc[cells_to_analyse]
#%%
depths_detected = data_summary.loc[cells_to_analyse]['Depth']
depths_undetected = data_summary.loc[remaining_cells]['Depth']

IN = subframe['Cell Type'] == 'IN'
PC = subframe['Cell Type'] == 'PC'

IN_depth = subframe[IN]['Depth']
IN_cells = subframe[IN]['Cell']
IN_amp = subframe[IN]['JTA Peak-Peak Amplitude']

PC_depth = subframe[PC]['Depth']
PC_cells = subframe[PC]['Cell']
PC_amp = subframe[PC]['JTA Peak-Peak Amplitude']
#%%color-code PCs, INs, and unclassified for supp figure
cell_type_color = np.empty((len(data_summary),1), dtype = object)

for row in range(len(cell_type_color)):
    if row in cells_to_analyse:
        if data_summary.loc[row]['Cell Type'] == 'PC':
            cell_type_color[row] = '#FD00FF'
        elif data_summary.loc[row]['Cell Type'] == 'IN':
            cell_type_color[row] = '#00FFFF'
    else:
        cell_type_color[row] = '#929191'
#%%Supplementary figure 4-3: depth vs distance and cell type vs distance
cm = plt.get_cmap('cool',3)
#Supplement to Figure 4-3A - Depth vs distance
fig = plt.figure(figsize=cm2inch(10, 10))

ax1 = fig.add_subplot(111)
for row in range(len(cell_type_color)):
    ax1.scatter(data_summary.loc[row]['Depth'], data_summary.loc[row]['Distance'],c = cell_type_color[row][0], edgecolor='k', lw=0.5, s = 30)
ax1.grid(b=True, alpha = 0.4)

#col = plt.colorbar(depth_plot ,orientation = 'horizontal', ticks = [-1, 0, 1], shrink = 0.75, pad = 0.17)
#col.set_label('Putative Cell Type ', fontweight = 'bold')
#col.ax.tick_params(labelsize=9)   

plt.yticks(np.arange(0, np.max(data_summary['Distance'])+1, 20.0  ))
plt.subplots_adjust(top = 0.95, bottom = 0.15, left = 0.22, right = 0.92)
plt.tick_params(labelsize=10)
plt.xticks([0, 400, 800, 1200, 1600, 2000])
plt.ylabel('Distance from Channel (um)', fontsize=11,labelpad = 5, fontweight = 'bold')
plt.xlabel('Cortical Depth (um)', fontsize=11, labelpad = 5, fontweight = 'bold')
plt.show()
plt.subplots_adjust(top = 0.90, bottom = 0.2, left = 0.20, right = 0.95)
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Supp_Fig 3_distance_vs_depth.png', dpi = 1200)
#%%Calculate firing rate over whole session
cell_rows_without_EAP = remaining_cells
cell_rows_with_EAP = cells_to_analyse

rates_dark = data_summary.loc[cell_rows_without_EAP]['# Patch Spikes' ] / (60*data_summary.loc[cell_rows_without_EAP]['Recording Length (mins)' ])
                             
rates_light = data_summary.loc[cell_rows_with_EAP]['# Patch Spikes' ] / (60*data_summary.loc[cell_rows_with_EAP]['Recording Length (mins)' ])
#%%Supplement to Figure 4-3C - amplitude vs firing rate over whole session
plt.scatter (data_summary['# Patch Spikes' ] / (60*data_summary['Recording Length (mins)' ]), data_summary['JTA Peak-Peak Amplitude' ])
                          
fig = plt.figure(figsize=cm2inch(10, 10))

ax1 = fig.add_subplot(111)
good_cells_plot = ax1.scatter(data_summary['# Patch Spikes' ][cells_to_analyse] / (60*data_summary['Recording Length (mins)' ][cells_to_analyse])  ,data_summary['JTA Peak-Peak Amplitude'][cells_to_analyse], c = '#F99120', edgecolor='k', lw=0.2, s = 40)
bad_cells_plot = ax1.scatter(data_summary['# Patch Spikes' ][remaining_cells] / (60*data_summary['Recording Length (mins)' ][remaining_cells]) ,data_summary['JTA Peak-Peak Amplitude'][remaining_cells], c = '#807F7F', edgecolor='k', lw=0.2, s = 40, alpha = 0.4)
#ax1.grid(b=True, alpha = 0.4)

#col = plt.colorbar(depth_plot ,orientation = 'horizontal', ticks = np.arange(300,2100,300), shrink = 0.75, pad = 0.17)
#col.set_label('Cortical Depth (um)', fontweight = 'bold')
#col.ax.tick_params(labelsize=9)   

#plt.xticks(np.arange(0, np.max(data_summary['Depth'])+1, 400.0  ))
plt.subplots_adjust(top = 0.98, bottom = 0.25, left = 0.25, right = 0.98)
plt.tick_params(labelsize=10)
plt.xlabel('Firing rate over whole session (Hz)', fontsize=11,labelpad = 5, fontweight = 'bold')
plt.ylabel('Amplitude (uV)', fontsize=11, labelpad = 5, fontweight = 'bold')
plt.yticks(rotation=90)
plt.show()
#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Supp_Fig 3_fr_vs_pk2pk.png', dpi = 1200)
#%%Compute Pearson correlation between PSTA amplitude and firing rate
import scipy.stats
print scipy.stats.spearmanr(data_summary['# Patch Spikes' ] / (60*data_summary['Recording Length (mins)' ]), data_summary['JTA Peak-Peak Amplitude' ])
                                         
plt.scatter( data_summary['# Patch Spikes' ] / (60*data_summary['Recording Length (mins)' ]), data_summary['JTA Peak-Peak Amplitude' ])
#%%Supplement to Figure 4-3B - Cell type distance vs amplitude
cells_above_10uV = data_summary.index[data_summary['JTA Peak-Peak Amplitude'] >= 10].tolist()
PCs = data_summary.index[data_summary['Cell Type'] == 'PC'].tolist()
INs = data_summary.index[data_summary['Cell Type'] == 'IN'].tolist()
#
excluded_cells = [2, 6, 9, 10, 12, 21, 23]
PCs2 = [ cells_above_10uV[i] for i in range(len(cells_above_10uV)) if i not in excluded_cells]
PCs_plot = [i for i in PCs2 if i not in INs ]
INs_plot = [i for i in INs if i in cells_above_10uV]
#           
fig = plt.figure(figsize=cm2inch(10,10))
ax1 = fig.add_subplot(111)
cm = plt.get_cmap('cool',2)

PCs_scatter = ax1.scatter(data_summary['Distance'][PCs_plot] ,data_summary['JTA Peak-Peak Amplitude'][PCs_plot], c = '#FD00FF', edgecolor='k', lw=0.2, s = 40)
INs_scatter = ax1.scatter(data_summary['Distance'][INs_plot] ,data_summary['JTA Peak-Peak Amplitude'][INs_plot], c = '#00FFFF', edgecolor='k', lw=0.2, s = 40)


plt.xlabel('Distance from channel (um)', labelpad = 5, fontsize = 11, fontweight = 'bold')
plt.ylabel('Peak-peak amplitude (uV)', labelpad = 5, fontsize = 11, fontweight = 'bold')
plt.subplots_adjust(top = 0.90, bottom = 0.2, left = 0.20, right = 0.95)
#plt.tight_layout
plt.grid('on', alpha = 0.4)

#%%
plt.savefig( 'C:/Users/Andre Marques-Smith/Dropbox/Paired Recordings biorxiv/Figure 3_materials/' + 'Supp to Fig 3_distance_cell_types.png', dpi = 1200)