#%%Imports
#comment test
import scipy.signal as signal
import numpy as np
import stfio
import matplotlib.pyplot as plt
import os
from collections import defaultdict
import pandas

#Functions
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a

#b,a = butter_bandpass(200, 2000, 30000)

def patch_spikes(wcrecording, T, mode='vclamp'):
    import numpy as np
    if mode == 'vclamp':
        wcrecording = -1* wcrecording
    else:
        pass

    #samples = int(len(wcrecording)/no_chunks)
    #chunks = range(0, no_chunks*samples+samples, samples)
    #medians = [np.median(wcrecording [ chunks[i]:chunks[i+1] ]) for i in range(len(chunks)-1)] #[ chunks[i]:chunks[i+1] ]) for i in range(len(chunks)-1)]
    #SlidingCatSub = []
    #SlidingCatSub[ chunks[i]:chunks[i+1] ] = np.concatenate([ wcrecording[ chunks[i]:chunks[i+1] ] - medians[i] for i in range(len(chunks)-1)])
    
    juxta_threshold = T*np.std(wcrecording)
    ThreshSub = wcrecording - juxta_threshold
    ThreshCrosses = np.diff(np.sign(ThreshSub))
    SpikeWstart = []
    SpikeWend = []
    SpikeWstart.append(np.where(ThreshCrosses == 2))
    SpikeWend.append(np.where(ThreshCrosses == -2))

    SpikeWDurs = [SpikeWend[0][0][i] - SpikeWstart[0][0][i] for i in range(len(SpikeWstart[0][0]))]
    
    duration_rejects = [index for index in range(len(SpikeWDurs)) if SpikeWDurs[index] < 5]
    
    a = SpikeWstart[0][0]
    b = SpikeWend[0][0]
    
    FinalWstart = [a[item] for item in range(len(a)) if item not in duration_rejects]
    FinalWend = [b[item] for item in range(len(b)) if item not in duration_rejects]
    
    FinalWDurs = [FinalWend[i] - FinalWstart[i] for i in range(len(FinalWstart))]

    
    wc_spike_samples = [list() for x in range(len(FinalWstart))]
    
    for i in range(len(wc_spike_samples)):
        wc_spike_samples[i].append( (FinalWstart[i] + np.argmax(wcrecording [FinalWstart[i] : FinalWend[i]])))
    wc_spike_samples = np.concatenate(wc_spike_samples)
    
        
    return wc_spike_samples, juxta_threshold

#%%
data_summary = pandas.read_excel('C:/Users/Andre Marques-Smith/Dropbox/Patch-Clamp Summary 2.0.xlsx')
input_dir = 'G:/transposed/'
listcells = os.listdir(input_dir)

#%% Select cell to process from listcells
cell = 35
pd_cell = listcells[cell]
paths = defaultdict(list)
#%%Set paths for each file of the chosen cell and load patch and extracellular recordings

cell_dir = input_dir + listcells[cell]
for file in os.listdir(cell_dir):
    if file.endswith('patch_ch1.bin'):
        paths['patch_v'].append(cell_dir + '/' + file)
    elif file.endswith('patch_sync.bin'):
        paths['patch_sync'].append(cell_dir + '/' + file)
    elif file.endswith('npx_raw.bin'):
        paths['npx_v'].append(cell_dir + '/' + file)
    elif file.endswith('npx_sync.bin'):
        paths['npx_sync'].append(cell_dir + '/' + file)
    elif file.endswith('meta.npy'):
        paths['meta'].append(cell_dir + '/' + file)
    else:
        pass
    
expt_meta = np.load(paths['meta'][0]).item()

np_spikes_cropped_test1 = np.memmap( paths['npx_v'][0], mode = 'c', dtype=np.int16, order = 'C')
np_spikes_cropped_test1 = np_spikes_cropped_test1.reshape((expt_meta['npx'][0][0], expt_meta['npx'][0][1]), order = 'F')
npx_sync = np.fromfile( paths['npx_sync'][0],dtype='int16')

wc_ch1_final_test1 = np.fromfile(paths['patch_v'][0] ,dtype='float64')
patch_sync = np.fromfile(paths['patch_sync'][0] , dtype='float64') 

#%%Offset subtraction  (already done) for shared data
#for chan in range(384):
#    np_spikes_cropped_test1[chan] -= np.median(np_spikes_cropped_test1[chan,:], axis = 0 ).astype(int)
#    print chan

#print 'Offset subtracted'
#%%High-pass filter at 200 Hz

Fcutoff = 200
Fsampling = 30000
Wn = np.float32(Fcutoff / (Fsampling / 2.0))
filterOrder = 6
#
b, a = signal.butter(filterOrder, Wn, btype='highpass', analog=0, output='ba')
#filteredData = signal.filtfilt(b, a, data, axis)

for chan in range(len(np_spikes_cropped_test1)):
    np_spikes_cropped_test1[chan,:] = signal.filtfilt(b, a, np_spikes_cropped_test1[chan,:])
    print chan

print 'High-pass filtering finished'
#Common median referencing
comm_avg_ref = np.median(np_spikes_cropped_test1[0:384], 0).astype(int)
#
for chan in range(len(np_spikes_cropped_test1)):
    np_spikes_cropped_test1[chan] -= comm_avg_ref
    print chan

print 'Common average referencing finished'
#%%High-pass filter the patch-clamp recording and run spike detection function
Fcutoff = 100
Fsampling = 30000
Wn = np.float32(Fcutoff / (Fsampling / 2.0))
filterOrder = 6
#
b, a = signal.butter(filterOrder, Wn, btype='highpass', analog=0, output='ba')
    
filter_patch = signal.filtfilt(b,a, wc_ch1_final_test1 )    
pd_idx = np.where(data_summary['Cell'] == pd_cell)
patch_type = (data_summary.loc[pd_idx]['Patch Type']).tolist()[0][-2:]

if patch_type == 'VC':
    spikemode = 'vclamp'
elif patch_type == 'IC':
    spikemode = 'iclamp'
elif patch_type == 'WC':
    spikemode = 'iclamp'
#
wc_spike_samples, juxta_thresh = patch_spikes(filter_patch, 7, mode=spikemode)

m = len(wc_ch1_final_test1)/float(len(np_spikes_cropped_test1[0]))
#m is the conversion factor to go from patch sample to neuropixel sample
np_buffer = 60 #100 for 3ms either side, 60 for 2ms

wc_buffer = int(np_buffer * m)

sta_wc = []
for i in range(len(wc_spike_samples)):
    sta_wc.append(  wc_ch1_final_test1[ wc_spike_samples[i] - wc_buffer : wc_spike_samples[i] + wc_buffer + 1] )

sta_wc_mean = np.mean(sta_wc, axis = 0)
sta_wc_timebase = np.arange(-wc_buffer/(m*30.0) , (wc_buffer+1)/(m*30.0), 1/(m*30.0))

sta_samples = (int(wc_spike_samples[0]/m) +1 + np_buffer) - (int(wc_spike_samples[0]/m) - np_buffer)
sta_windows_chan_array = np.empty((384, sta_samples, len(wc_spike_samples)),dtype='float')
#%%Organise an array with patch spike snippets (+- 2ms from spike peak) for all patch spikes, per np channel

for wc_spike in range(len(wc_spike_samples)):
    for channel in range(384):
        sta_windows_chan_array[channel,:,wc_spike] = (np_spikes_cropped_test1[channel] [int(wc_spike_samples[wc_spike]/m) - np_buffer : int(wc_spike_samples[wc_spike]/m) + np_buffer + 1 ])
        print wc_spike*100.0/len(wc_spike_samples)
#%%Compute mean STA on patch spike per np channel, convert from bit to microvolt, compute and sort chans by pk2pk amp to find channel with highest EAP amplitude
sta_np_by_channel = np.empty((384, len(sta_windows_chan_array[0,:,0])),dtype='float64' )
for chan in range(len(sta_windows_chan_array)):
    sta_np_by_channel[chan,:] = 2.34375*np.mean(sta_windows_chan_array[chan,:,:], axis = 1)


#
pk2pk = []
for i in range(384):
    pk2pk.append( np.max(sta_np_by_channel[i][int(0.75*np_buffer):int(1.25*np_buffer)]) - np.min(sta_np_by_channel[i][int(0.75*np_buffer):int(1.25*np_buffer)]))

channels_by_pk = np.empty((384,2))
channels_by_pk[:,0] = range(384)
channels_by_pk[:,1] = pk2pk
channels_by_pk=-channels_by_pk
channels_by_pk=-channels_by_pk[channels_by_pk[:,1].argsort()]
#
print 'number of patch spikes is', len(wc_spike_samples), 'spikes'
print "# channels over 20 uV is", np.sum(channels_by_pk[:,1] > 20)
print 'recording length is', len(npx_sync) * 1/30000/60.0, 'minutes'
print 'juxta threshold is', juxta_thresh
print 'central channel is', int(channels_by_pk[0,0]), 'with amplitude', channels_by_pk[0,1]
print 'cell just finished is', pd_cell
#plot average waveforms for 8 channels with highest pk2pk
labels = channels_by_pk[0:8,0].astype(int)
for i in range(8):         
    plt.plot( sta_np_by_channel[channels_by_pk[i,0].astype(int)], label = labels[i])    
    plt.legend()

#%%Output files for further analysis
output_dir = 'E:/Code/for analysis/' + listcells[cell] +'_aligned_analysis'
os.makedirs(output_dir)
#
np.save(output_dir + '/%s_expt_meta.npy' %(listcells[cell]), expt_meta)

np_spikes_cropped_test1.tofile(output_dir + '/%s_npx_preprocessed.bin'  %(listcells[cell]), sep="")
#
np.save(output_dir + '/%s_patch_preprocessed'  %(listcells[cell]), filter_patch)
np.save(output_dir + '/%s_patch_original'  %(listcells[cell]), wc_ch1_final_test1)
np.save(output_dir + '/%s_sta_windows_chan_array' %(listcells[cell]), sta_windows_chan_array)
np.save(output_dir + '/%s_sta_np_by_channel' %(listcells[cell]), sta_np_by_channel)

np.save(output_dir + '/%s_wc_spike_samples' %(listcells[cell]), wc_spike_samples)
np.save(output_dir + '/%s_sta_wc' %(listcells[cell]), sta_wc)

np.save(output_dir + '/%s_channels_by_pk' %(listcells[cell]), channels_by_pk)
#%%
