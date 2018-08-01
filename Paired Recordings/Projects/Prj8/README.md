# Backpropagating Action Potentials and Multi-channel Spike Waveforms  
### High-density CMOS probes provide high spatial sampling of extracellular waveforms, as can be seen from the examples we provide in Figure 7 and its supplements. Recently, Jia et al (2018, biorXiv) from the Allen Brain Institute provided a compelling demonstration of how this feature of Neuropixels probes can be useful. Here, I would like to extend the analysis of spatiotemporal dynamics of EAP waveforms, looking at backpropagation of action potentials, their trajectory motifs and reliability.  

Andre Marques-Smith: I am re-working Figure 7 and splitting it into two. The new Fig 7 will be dedicated exclusively to furthering our previous analysis of EAP backpropagation, inspired by the findings of Jia et al (2018, biorXiv).

New Fig 7A:  
![Fig 7A](https://github.com/kampff-lab/sc.io/blob/Project_8/Paired%20Recordings/Projects/Prj8/new_fig_7A.png)  
Here we are looking at EAP propagation along a single column of Neuropixel channels (corresponding to the somatic channel's column). The scatterplot is distance to soma versus EAP amplitude relative to somatic channel. The black line is the average EAP propagation profile of 21 cells. I have determined also the distance at which EAP amplitude attenuates to 50%, 25% and 12% of the maximum amplitude. The latter was done by interpolating the propagation profile (black line) to 1 µm intervals, then finding the intersection with lines y = 0.5, 0.25 or 0.12.

Next steps:
- Based on the average profile, I decided to focus on the spatial propagation window of -240 µm to +240 µm above the soma;
- I will re-draw the Joy Division plots (old Fig. 7A, 7B and Supplements to Fig 7) to focus only on this window, as propagation outside it is very limited and consists mostly of noise;
- For every cell, I will plot also its propagation trajectory, as done by Jia et al. This plot will show for every channel in the column of -240:240µm distance to the soma the latency between negative EAP peak at the soma and negative EAP peak at that channel. I like these plots :)
- Preliminary plotting suggests we have the same pattern of some cells having uni- and others bi-directional propagation; second, our Joy Division plots already show the 5 putative interneurons have very different propagation profiles;
- I am interested in asking whether for the same neuron, EAP propagation always follows the same motif. To ensure we can examine single spikes, this analysis will have to be restrained to the 10 cells with higher EAP amplitude (>50 µV), perhaps even only to a subset of these with particularly high amplitudes. 
    - If the mechanistic possibility exists for AP backpropagation to be modulated (e.g. by inhibitory interneurons or cell-intrinsic mechanisms), what could it be exploited for? Main hypothesis I'm aware of is plasticity, but contextual modulation/gating of top-down inputs sounds interesting as well.
    - On the other hand, AP backpropagation could be stereotyped for a cell*. Presumably, differences in backprop motifs arise from a combination of i) biophysical properties, ii) morphology and iii) orientation of that neuron's neurites relative to the probe. One might imagine the particular combination of these 3 properties embodied by a cell to be quite specific; with that in mind, can multi-channel backprop motifs be exploited as an additional dimension for spike-sorting? Is variability in motif between cells greater than within cells? This could be interesting to explore.

*_This could be stereotyped in anaesthetised state but not in awake state_

Questions I have (help!):
- These plots are limited, for clarity, to a single column of channels. I'm having trouble working out a trajectory representation that encompasses all 4 columns, as EAP negative peaks are often simultaneous for channels on the same row;
- What would be a good method for determining whether within-cell variability in propagation motif is greater/lower than between-cell variability?
