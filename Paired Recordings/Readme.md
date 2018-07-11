# Introduction

This is the companion repository to our preprint where we describe experiments recording from the same cortical neuron *in vivo* using [Neuropixel probes](https://www.nature.com/articles/nature24636) and patch-clamp:

**Marques-Smith, A., Neto, J.P., Lopes, G., Nogueira, J., Calcaterra, L., Frazão, J., Kim, D., Phillips, M., Dimitriadis, G., Kampff, A.R. (xx July 2018)** *Recording from the same neuron with high-density CMOS probes and patch-clamp: a ground-truth dataset and an experiment in collaboration.* biorXiv DOI & URL.

In this repository you will find:
- The code used to preprocess data, run analysis and generate the figures presented in Marques-Smith et al., 2018;
- A list of proposed collaborative follow-up projects based on that same dataset and details on how to contribute;
- Instructions on navigating the dataset

## Conditions of use
All the data we shared is free for you to use for scientific research purposes. We ask only that you cite both:
- The original publication describing the dataset: Marques-Smith, A., Neto, J.P., Lopes, G., Nogueira, J., Calcaterra, L., Frazão, J., Kim, D., Phillips, M., Dimitriadis, G., Kampff, A.R. (xx July 2018). *Recording from the same neuron with high-density CMOS probes and patch-clamp: a ground-truth dataset and an experiment in collaboration.* biorXiv DOI & URL.
- The dataset itself: 

## Dataset and Metadata
The full dataset is freely-available for download [here](https://drive.google.com/open?id=13GCOuWN4QMW6vQmlNIolUrxPy-4Wv1BC).

That directory also contains:
- A spreadsheet with metadata for each paired-recording;
- A document 

### Dataset organisation
Each folder is titled "cxx' (cell xx), corresponds to a paired-recording and contains the following files:

- cxx_expt_meta.csv  
*Dimensions that Neuropixel and patch-clamp recording files should be reshaped to and their data type (int, float), when importing into Python or Matlab.*

- cxx_npx_raw.bin  
*Neuropixel recording (384 channels), 1D binary file.*

- cxx_npx_sync.bin  
*Neuropixel sync channel, binary file.*

- cxx_patch_ch1.bin  
*Patch-clamp current or voltage, depending if cell was recorded in voltage- or current-clamp.*

- cxx_patch_sync.bin  
*Patch-clamp sync channel.*

- cxx_wc_spike_samples.npy  
*Patch-clamp recording samples corresponding to patched cell's spike peaks, numpy 1D array.*

- cxx_extracellular_spikes.npy   
*Only for the 21 cells where an extracellular spike waveform could be detected. Neuropixel recording samples corresponding to the peak of extracellular spikes in the channel closest to patched cell's soma, numpy 1D array.*

### How to load data
*Neuropixel*  
Neuropixel recordings were saved into a 1D binary vector, from a matrix organised as 384 rows (channels) x n columns (samples). The data was written to file from this matrix in column-major (F) order, ie, the first sample in the recording was written to file for every channel, then the second sample was written for every channel, etc, etc. To load the data correctly you have to tell your software how it's organised and how to unpack it. The file titled cxx_expt_meta.csv provides the info you need for this.

In Python, you could do:  
```python 
import numpy as np

npx_path = 'full path to the downloaded cxx_npx_raw.bin file'
npx_channels = 384

npx_recording = np.memmap( npx_path, mode = 'r', dtype=np.int16, order = 'C')

npx_samples = len(npx_recording)/npx_channels

npx_recording = npx_recording.reshape((npx_channels, npx_samples), order = 'F')
```
which will result in a 2D array where channels are rows (384) and columns are samples. Pay attention to the ['mode' parameter](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.memmap.html) in the memmap function, as whis will determine whether the file is read-only or if you're writing changes to disk.

**Neuropixel raw data is provided as an int16 binary. [Neuropixel ADCs](https://github.com/cortex-lab/neuropixels/wiki/Gain_settings) are 10 bits, with a range of -0.6 to 0.6 V, and acquisition was at 500x gain, yielding a resolution of 2.34 µV/bit.  
To obtain readings in µV, you should  multiply the int16 values by 2.34.** 

*Patch-clamp*
Patch-clamp data are already provided in float64 type and current (pA) or voltage (mV) units, depending if the recording was performed in voltage- or current-clamp.

To load a patch-clamp recording in python, do:  
```python
import numpy as np

patch_path = 'full path to the downloaded cxx_patch_ch1.bin file'
patch_recording = np.fromfile(patch_path, dtype='float64')
```
