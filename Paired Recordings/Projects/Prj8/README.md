# Backpropagating Action Potentials and Multi-channel Spike Waveforms  
### High-density CMOS probes provide high spatial sampling of extracellular waveforms, as can be seen from the examples we provide in Figure 7 and its supplements. Recently, Jia et al (2018, biorXiv) from the Allen Brain Institute provided a compelling demonstration of how this feature of Neuropixels probes can be useful. Here, I would like to extend the analysis of spatiotemporal dynamics of EAP waveforms, looking at backpropagation of action potentials, their trajectory motifs and reliability.  

Andre Marques-Smith: I am re-working Figure 7 and splitting it into two. The new Fig 7 will be dedicated exclusively to furthering our previous analysis of EAP backpropagation, inspired by the findings of Jia et al (2018, biorXiv).

New Fig 7A - EAP propagation profile  
![Fig 7A](https://github.com/kampff-lab/sc.io/blob/Project_8/Paired%20Recordings/Projects/Prj8/Fig_7A_EAP_propagation_profile.png)  
Here we are looking at EAP propagation along a single column of Neuropixel channels (corresponding to the somatic channel's column). The scatterplot is EAP amplitude relative to somatic channel vs distance to soma. The black line is the average EAP propagation profile of 21 cells. I have determined also the distance at which EAP amplitude attenuates to 50%, 25% and 12% of the maximum amplitude. The latter was done by interpolating the propagation profile (black line) to 1 µm intervals, then finding the intersection with lines y = 0.5, 0.25 or 0.12.

New Fig 7B - average EAP propagation trajectory  
![Fig 7B](https://github.com/kampff-lab/sc.io/blob/Project_8/Paired%20Recordings/Projects/Prj8/Fig_7B_Trajectories.png)  
This illustrates the propagation trajectory for each cell's average multi-channel EAP waveform. I obtained local minima at each channel and plotted the time difference of that to the somatic spike trough versus the distance to the soma.

Next steps:
- Showcase and describe examples of cells with representative/interesting EAP trajectories (e.g. INs vs PYRs).
- Calculate trajectory features, i.e. velocity above vs below soma.
- Investigate within-cell single-spike variability of EAP trajectory in high peak-peak cells.
- Is within-cell variability > or < than between-cell?

*_This could be stereotyped in anaesthetised state but not in awake state_
