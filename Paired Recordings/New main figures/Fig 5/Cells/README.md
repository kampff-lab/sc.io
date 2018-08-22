# Spike detectability
Here we are asking what are the limits of spike detectability given that in GT experiments we know the waveform(s) we're looking for.    

### First pass
Following a suggestion by Nick Steinmetz, I looked at a single channel and computed the dot product between the template (average spike waveform) of that channel and each of the individual extracellular GT spikes. I plotted the frequency distribution of this and compared with that obtained from the dot product between the same template and random non-GT spikes detected in that channel. These are the plots shown in this folder.    

See below for two examples. Frequency distribution of template x GT spike dot product is in purple, and versus non-GT random spike is in yelloy.    

This is the plot for c46, a unit with 246 µV peak-peak.
![Easy unit](https://github.com/kampff-lab/sc.io/blob/master/Paired%20Recordings/New%20main%20figures/Fig%205/Cells/dot_prod_distribution_c46.png)

This is the plot for c7, a unit with 12 µV peak-peak.
![Difficult unit](https://github.com/kampff-lab/sc.io/blob/master/Paired%20Recordings/New%20main%20figures/Fig%205/Cells/dot_prod_distribution_c7.png)

Can the two distributions be separated and how small can a spike be for this to remain true? Preliminary analysis suggests even units with small (10-20 µV) spikes can be separated from background (albeit with less success), given a known template to look for.

### Further work
1. Improving this analysis following other suggestions from Nick Steinmetz. 
2. Synthetically increasing or decreasing spike/template amplitude to figure out i) how small a spike can be and still be separable from background or ii) how much bigger would it need to have been for good separation.
3. How about multi-channel waveforms? With high-density probes, spikes from the same neuron (even low amplitude units) have a footprint spanning multiple (4-30) channels. In principle, multi-channel templates should improve the ability to detect and separate low-amplitude spikes from background, if we know the template to look for.
