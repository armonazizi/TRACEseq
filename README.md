# TRACEseq Code
### This repository contains the code associated with Sharma et al. Nature Communications 2021
[link to paper](https://www.biorxiv.org/content/10.1101/2020.05.25.115329v1)

The code in the repository generates all bubble plots shown in the publication (Figures 4,5,Extended Figures 3,5).

To run the code: 
  1. Download or clone the repository locally.
  2. Open *generate_bubble_plots.R* in Rstudio.
  3. Navigate to the local directory containing the repository code.
  4. Run code (this will source *bubble_plot_functions.R* and generate all plots)
  
The following libraries must be installed for the code to function.
* openxlsx
* ggpubr
* packcircles
* ggplot2
* grid


The script *align_and_filter_pre_tuba_HBB.py* is utilized in our analyses within the TUBAseq pipeline as described in the manuscript. The script filters a fastq file and bins reads into HR, NHEJ, Unmodified, and Other categories.
