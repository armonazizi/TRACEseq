# Armon Azizi (aazizi@stanford.edu)
#
# Bubble plot code associated with Sharma et al. 20XX
#
# The R code in this file contains functions and calls to generate the bubble plots in the publication.

require(ggpubr)
require(packcircles)
require(ggplot2)

#' Generate a cloud of circles to represent cell population
#' Modified from  clonevol v0.99.11
#' @param colors: matching colors of the cells. Length(colors) will be
#' used as the number of cells to generate
#' @param maxiter: number of iterations for the layout algorithm
#' of the packcircles package
#
generate.cloud.of.cells <- function(colors, maxiter=1000){
  
  n = length(colors)
  limits = c(-50,50)
  inset = diff(limits)/3
  # scale radius to make sure cloud of cells gather approximately like a sphere
  # with 200 cells, limits = c(-50,50), radius ~ 3.3 makes sure circles
  radius = sqrt(11*200/n)
  xyr = data.frame(
    x = runif(n, min(limits) + inset, max(limits) - inset),
    y = runif(n, min(limits) + inset, max(limits) - inset),
    r = rep(radius, n))
  
  # Next, we use the `circleLayout` function to try to find a non-overlapping
  # arrangement, allowing the circles to occupy any part of the bounding square.
  # The returned value is a list with elements for the layout and the number
  # of iterations performed.
  res = circleLayout(xyr, limits, limits, maxiter = maxiter)
  
  ## plot data for the `after` layout returned by circleLayout
  cells = circlePlotData(res$layout)
  
  # color the circles
  cells$color = sample(colors, nrow(cells), replace=TRUE)
  
  return(cells)
}

#' Plot a cloud of cells
#' modified from  clonevol v0.99.11
#' @param cells: cell cloud as returned from generate.cloud.of.cells
#' @param colors: vector of colors to use to color bubbles
#' @param frame: draw a frame surrounding the cloud of cells
#' @param cell.border.color: color of the border of the circles used
#' to draw a cell (default = black)
#' @param cell.border.size: line size of the cell border, the smaller
#' the figure size is, the smaller this value needs to be (default = 0.1)
#' @param clone.grouping: how the cells of the same clone being grouped
#' values are c("random", "horizontal", "vertical"), default="random"
#' @param limits: plot limits for ggplot frame. Modifying these limits allows the user to change the relative size of the plot.
#' @import ggplot2
plot.cloud.of.cells.mod <-function(cells, colors, title='', title.size=1,
                                   alpha=1, frame=FALSE,
                                   cell.border.color='black', cell.border.size=0.1,
                                   clone.grouping='random', limits=c(-50,50)){
  
  diameter<-max(cells[cells$id==1,]$y)-min(cells[cells$id==1,]$y)
  collapsed_data<-aggregate(cells[,c("x","y")],list(cells$id),mean)
  colnames(collapsed_data)<-c("id", "x", "y")
  collapsed_data$color<-colors
  collapsed_data$color[order(collapsed_data$y, collapsed_data$x)]<-colors
  
  cells_colored<-merge(cells[,1:3], collapsed_data[, c("id", "color")], by="id")
  
  p = (ggplot(cells_colored) +
         geom_polygon(aes(x, y, group=id, fill=color), color=cell.border.color,
                      alpha=alpha, size=cell.border.size) +
         coord_equal(xlim=limits, ylim=limits) +
         theme_bw() +
         theme(axis.ticks.length=unit(0.1, 'mm')) +
         theme(axis.text=element_blank(),
               axis.ticks=element_blank(),
               axis.title=element_blank(),
               legend.position='none',
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank()) +
         #theme_void() +
         theme(plot.margin=unit(c(0,0,0,0), 'mm')) +
         #labs(title=title) +
         #scale_fill_manual(values=colors) +
         scale_fill_identity()
  )
  if (title != ""){
    p = p + ggtitle(title) + theme(plot.title=element_text(size=10*title.size))
  }
  
  if (!frame){
    p = p + theme(panel.border=element_blank())
  }else{
    p = p + theme(panel.border=element_rect(linetype='dotted'))
  }
  return(p)
}


# Multiplot-Function
# From: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow=TRUE)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#' generates a set of bubble plots given a counts matrix and list of samples
#' @param barcode_matrix: dataframe/matrix of counts. Rows are barcodes and columns are samples
#' @param list_of_samples: character vector of sample names corresponding to column names in barcode_matrix
#' @param num_barcodes: number of barcodes to plot from each sample. i.e. if n=3, the top 3 barcodes from each sample will be plotted on all bubble plots.
#' @param num_cells: number of cells to plot for each sample
#' @param title.size: title size
#' @param titles: titles to plot for each sample (same order as list_of_samples)
#' @param cutoff_threshold: if provided, fractions below this threshold are not plotted.
#' @param numcolumns: number of columns in final plot
#' @param pallette: can provide ggpubr pallette
#' @param colors: set of colors to use (must be long enough for all barcodes)
#' @param save: save plots?
#' @param saveprefix: file prefix to use for saving pdfs
#' @param returnplots: should the function return all plots as a list of ggplot objects?
#' @param returnvalues: instead of plotting, return a dataframe of counts that would have been plotted.
make_bubble_set_multiplot_other<-function(barcode_matrix, 
                                          list_of_samples, 
                                          num_barcodes, 
                                          num_cells, 
                                          title.size=1, 
                                          titles=NULL,
                                          cutoff_threshold=NULL,
                                          numcolumns=4, 
                                          palette='Dark2',
                                          colors=NULL,
                                          save=FALSE,
                                          saveprefix=NULL,
                                          returnplots=FALSE,
                                          returnvalues=FALSE){
  
  # get list of samples with nonzero barcodes
  nonzero_samples<-list_of_samples[colSums(barcode_matrix[,list_of_samples])>0]
  
  if(is.null(titles)){titles<-samples}
  
  # generate list of barcodes
  list_of_barcodes<-c()
  for(sample in nonzero_samples){
    list_of_barcodes<-c(list_of_barcodes, rownames(barcode_matrix)[rev(order(barcode_matrix[,sample]))][1:num_barcodes])
  }
  
  # subset barcode matrix
  new_barcode_matrix<-barcode_matrix[unique(list_of_barcodes),nonzero_samples]
  
  # Add "other" barcodes
  new_barcode_matrix["other",]<-colSums(barcode_matrix[!(rownames(barcode_matrix)%in%unique(list_of_barcodes)),nonzero_samples])
  new_barcode_matrix<-sweep(new_barcode_matrix,2,colSums(new_barcode_matrix),"/")
  
  # if cutoff is provided, remove fractions below threshold.
  if(!is.null(cutoff_threshold)){
    new_barcode_matrix[new_barcode_matrix<cutoff_threshold]<-0
  }
  
  new_barcode_matrix<-floor(new_barcode_matrix*num_cells)
  new_barcode_matrix<-new_barcode_matrix[rowSums(new_barcode_matrix)>0,]
  
  # set colors
  if("other"%in%rownames(new_barcode_matrix)){
    new_barcode_matrix$color<-c(get_palette(palette, k=length(rownames(new_barcode_matrix))-1),"#CCCCCC")
    if(!is.null(colors)){
      if(length(colors)<length(rownames(new_barcode_matrix))-1){errorCondition("Not enough colors provided")}
      new_barcode_matrix$color<-c(colors[1:length(rownames(new_barcode_matrix))-1],"#CCCCCC")
    } 
  }else{
    new_barcode_matrix$color<-get_palette(palette, k=length(rownames(new_barcode_matrix)))
    if(!is.null(colors)){
      if(length(colors)<length(rownames(new_barcode_matrix))){errorCondition("Not enough colors provided")}
      new_barcode_matrix$color<-colors[1:length(rownames(new_barcode_matrix))]
    }
  }
  
  if(returnvalues){return(new_barcode_matrix)}
  
  plots <- list()  # new empty list
  n=1
  for(sample in list_of_samples){
    
    if(sample %in% nonzero_samples){
      color_vector<-c()
      
      for(bc in rownames(new_barcode_matrix)){
        color_vector<-c(color_vector, rep(new_barcode_matrix[bc,"color"], new_barcode_matrix[bc,sample]))
      }
      
      
      cell_cloud<-generate.cloud.of.cells(color_vector, maxiter = 100)
      
      p1<-plot.cloud.of.cells.mod(cell_cloud, colors = color_vector, clone.grouping="horizontal",title=titles[n], title.size)
      plots[[n]] <- p1 #add plot to plot list
      n <- n+1
    } else {
      message(paste(sample," - Error: Sum of Barcodes in Sample = Zero"))
      cell_cloud<-generate.cloud.of.cells(rep("#f7f7f7", num_cells), maxiter = 100)
      
      p1<-plot.cloud.of.cells.mod(cell_cloud, colors = rep("#f7f7f7", num_cells), clone.grouping="horizontal",title=titles[n], title.size)
      plots[[n]] <- p1 #add plot to plot list
      n <- n+1
    }
    
  }
  
  if(!save){
    multiplot(plotlist = plots, cols=numcolumns)
  }else{
    pdf(paste0(saveprefix,".bubbleplot.pdf"), width=numcolumns, height=ceiling(length(plots)/numcolumns))
    multiplot(plotlist = plots, cols=numcolumns)
    dev.off()
    pdf(paste0(saveprefix,".legend.pdf"), width=7, height=5)
    plot.new()
    legend("topleft",
           legend=rownames(new_barcode_matrix), 
           fill=new_barcode_matrix$color)
    dev.off()
  }
  
  if(returnplots){
    return(plots)
  }
}
