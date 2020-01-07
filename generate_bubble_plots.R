# Armon Azizi (aazizi@stanford.edu)
#
# Bubble plot code associated with Sharma et al. 20XX
#
# The R code in this file contains calls to generate the bubble plots in the publication.

source("bubble_plot_functions.R")

library(openxlsx)

##### LINEAGE PLOTS HBB #####

# Figure 4 Bubble plots
# Read barcode matrix, generate colors and samples of interest, and run plotting function.
barcode_matrix<-read.table("data/hbb_invivo_normalized_counts.tsv", header = TRUE, row.names = "barcode")
colors<-c("#1B9E77","#AE6D1C","#E41A1C","#000099","#9B58A5","#FF9933","#D8367D","#749829","#ABCFE5")
samples<-c("m20_Week_16_HBBbc_19","m20_Week_16_HBBbc_33","m20_Week_16_HBBbc_3410","m20_Week_12_Secondary_HBBbc_19","m20_Week_12_Secondary_HBBbc_33","m20_Week_12_Secondary_HBBbc_HSPC")
make_bubble_set_multiplot_other(barcode_matrix=barcode_matrix,
                                list_of_samples=samples,
                                num_barcodes = 3,
                                num_cells=200,
                                title.size=0.5,
                                numcolumns=6,
                                colors = colors)




# Extended data figure 3: all hbb bubble plots
# read in excel file of samples to plot
comp_matrix<-read.xlsx("data/hbb_comparison_matrix.xlsx")
hbb_labels<-comp_matrix

# generate list of samples and titles (label low input samples)
samples<-as.vector(unlist(hbb_labels))
samples<-samples[!is.na(samples)]
titles<-samples
low_input_samples<-c("m29_Week_16_HBBbc_3410","m31_Week_16_HBBbc_3410")
for(i in low_input_samples){
  titles[titles==i]<-paste0(i, "\nLOW_INPUT")
}

# generate large color set for all clones and shuffle
colors<-get_palette("ucscgb", k = 60)
set.seed(1)
colors<-sample(colors)

# generate bubble plots and plot a 3 column grid of samples.
make_bubble_set_multiplot_other(barcode_matrix=barcode_matrix,
                                list_of_samples=samples,
                                titles=titles,
                                num_barcodes = 3,
                                num_cells=200,
                                cutoff_threshold = NULL,
                                title.size=0.25,
                                numcolumns=3,
                                colors=colors,
                                save=FALSE,
                                returnplots = FALSE)




##### LINEAGE PLOTS BFP #####

# read in normalized count matrix for BFP/AAVS1 samples
BFP_barcode_matrix<-read.table("data/bfp_invivo_normalized_counts.tsv", header = TRUE, row.names = "barcode", check.names = FALSE)

# generate colors and samples
samples<-c("m7_Week_16_BFP_19+","m7_Week_16_BFP_33P","m7_Week_12_Secondary_BFP_19+","m7_Week_12_Secondary_BFP_33P")
colors<-get_palette("ucscgb",k=10)

# Run plotting function
make_bubble_set_multiplot_other(barcode_matrix=BFP_barcode_matrix,
                                list_of_samples=samples,
                                num_barcodes = 5,
                                num_cells=200,
                                title.size=0.5,
                                colors=colors,
                                numcolumns=4)


# read matrix of samples to plot.
# generate samples and titles list and label low input samples.
samples<-read.xlsx("data/bfp_comparison_matrix.xlsx")
samples<-as.vector(unlist(samples))
samples<-samples[!is.na(samples)]
titles<-samples
low_input_samples<-c("m38_Week_16_BFP_33P")
for(i in low_input_samples){
  titles[titles==i]<-paste0(i, "\nLOW_INPUT")
}
colors<-get_palette("ucscgb", k = 60)
set.seed(1)
colors<-sample(colors)

# run plotting function and generate plot grid
make_bubble_set_multiplot_other(barcode_matrix=BFP_barcode_matrix,
                                list_of_samples=samples,
                                titles=titles,
                                num_barcodes = 3,
                                num_cells=200,
                                title.size=0.25,
                                numcolumns=2,
                                colors=colors,
                                save=FALSE,
                                returnplots = FALSE)
