# Spike detectability
What are the limits of spike detectability in GT conditions where we know the waveform(s) we're looking for.    

### To-do
* For a single channel, compute dot product between the template (average spike waveform) of that channel and each of the individual extracellular GT spikes.
* For same channel, compute dot product between the template and segments of the recording.
* Plot the frequency distribution of dot product results in two conditions. Can the distributions be separated?
* Systematically vary threshold used in each case and plot ROC curve. 
**Outlook:** Potential for using generative models + template matching for spike detection.

# Spike propagation trajectory across multiple channels
Given a probe insertion close to parallel to the somato-dendritic axis of cortex, is it possible to track spike propagation of a unit across channels?


### To-do
* Work out AP spatial spread along D-V axis
* Map propagation trajectory for different units
* Analyse features of propagation in space
* Tie into sub-cellular AP signatures
**Outlook:** Viability of investigating BAPs in vivo; use of multi-channel EAP propagation stereotipy as an added dimension for sorting.


# Dark Neurons
What else can we tell about them from this dataset?

### To-do
* Analyse EAP and patch spike waveforms and compare to good units at similar distance range
**Outlook:** Further characterisation of dark neurons and their identity.


# Synaptic Connectivity
How can we take advantage of a paired-recording dataset to explore analyses for connectivity?

### To-do
* Spike sort recordings where whole-cell patch-clamp was achieved *and* there was a good EAP signature.
* For every unit in the recording, run STA on the whole-cell unit's membrane potential.
* For every putative pre-synaptic unit, run STA and obtain extracellular footprint; investigate for AxTPs
**Outlook:** Finding additional signatures for unit-unit connectivity.
