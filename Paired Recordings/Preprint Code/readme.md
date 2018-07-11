This folder hosts the code that was used to preprocess and analyse data and generate figures contained in Marques-Smith et al., 2018:

Marques-Smith, A., Neto, J.P., Lopes, G., Nogueira, J., Calcaterra, L., Frazão, J., Kim, D., Phillips, M., Dimitriadis, G., Kampff, A.R. (xx July 2018) **Recording from the same neuron with high-density CMOS probes and patch-clamp: a ground-truth dataset and an experiment in collaboration.** biorXiv DOI.


Code in files 2-7 relies on data structures that are created by 1. Briefly, 1 will high-pass filter and common-average reference the extracellular recording data, then high-pass filter the patch-clamp data and run spike detection on it. Subsequently it will create key data structures which form the basis for further analysis in files 2-7.

Files 2-7 have no dependencies on each other. You can therefore run them in any order, once you have executed and saved the outputs of 1.

You will also need to download 2 additional files in this folder:
- npx_map.csv: a schematic of Neuropixel channel layout;
- Data Summary.xlsx: spreadsheet containing required metadata for every recording.
