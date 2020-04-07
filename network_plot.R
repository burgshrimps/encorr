# Functions
load_neurons <- function(dataset) {
  ### Loads neurons (nodes of network graph) from CSV file and performs basic data correction and filtering ###
  
  neurons <- read.csv(paste(dataset, '_neurons.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
  
  # For now filter out neurons with area labels 'cortex' and 'x'
  neurons <- neurons[neurons$area != 'cortex' & neurons$area != 'x', ]  
  
  # Rename areas labels 'pCA3x' to 'pCA3'
  neurons[neurons == 'pCA3x'] <- 'pCA3'
  
  return(neurons)
}


load_neurons_effective <- function(filename) {
  ### Loads list of effective neurons from CSV for all datasets ###
  
  neurons_effective <- read.csv(paste(dir, filename, sep=''), header=TRUE, sep=',', check.names=FALSE)
  
  return(neurons_effective)
}


reorder_neurons_for_circ <- function(neurons) {
  ### Change order of neurons in regards to the order they appear in the network plot based on brain area ###
  
  neurons$area <- factor(neurons$area, levels = c('sub', 'dCA1', 'CA1b', 'pCA1b', 'pCA1', 'dCA3', 'pCA3'))
  neurons$area.num <- as.numeric(neurons$area)
  neurons <- neurons[order(neurons$area.num), ]
  
  return(neurons)
}


load_phase_connections <- function(dataset, phase, neurons) {
  ### Loads connections between neurons measered by cross-correlation from CSV ###
  
  phase_connections <- read.csv(paste(dataset, '_', phase, '.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
  phase_connections$type.num <- as.numeric(factor(phase_connections$type))
  
  # Make sure to only include connections between neurons which are actually in the dataset, e.g. not between neurons labeled with 'x'
  phase_connections <- phase_connections[phase_connections$ref %in% neurons$id & phase_connections$tar %in% neurons$id, ]
  
  return(phase_connections)
}


plot_circ <- function(neurons, areas, phase, title, dataset, shift, v_size, v_label_size) {
  ### Plots neurons as vertices in a circle with connections between them as edges ###
  
  # Create network graph from vertex and edge list
  net <- graph_from_data_frame(d=phase, vertices=neurons, directed=FALSE)
  
  # Shift coordinates of each vertex by a manual parameter such that graph layout roughly matches anatomical layout of hippocampus
  l <- layout_in_circle(net)
  end <- length(l[,1])+shift-1
  l <- rbind(l, l)[shift:end, ]
  
  # Set vertex size
  V(net)$size <- v_size
  
  # Define color scheme for vertices and edges
  v_color <- c('yellow', 'red', 'greenyellow', 'greenyellow', 'green', 'cyan', 'orange') 
  e_color <- c('red', 'blue')
  V(net)$color <- v_color[V(net)$area.num]
  E(net)$color <- e_color[E(net)$type.num]
  
  # Use absolute value of connection weight as width for edge if not NA
  if (is.null(E(net)$weight)){
    E(net)$width <- 0
  } else {
    E(net)$width <- abs(E(net)$weight)
  }
  
  # Plot network and save as PDF
  pdf(paste(dataset, '_network_',title, '.pdf', sep=''))
  plot(net, layout=l, vertex.label.cex=v_label_size)
  title(paste(dataset, title, sep=' '))
  legend("topleft", legend=areas, col=v_color[unique(neurons$area.num)], cex=0.7, pch=21, pt.bg=v_color[unique(neurons$area.num)])
  dev.off()
}


# Set dataset parameters
animal <- 'LE87'
date <- '20190520'
dir <- '/Users/burgshrimps/Documents/lin/analysis/XCORR/'
dataset <- paste(animal, date, sep='_')
wd <- paste(dir, animal, '/', date, '/conn_stat', sep='')
setwd(wd)


# Load neurons
library(plyr)
neurons <- load_neurons(dataset)
neurons_ordered_circ <- reorder_neurons_for_circ(neurons)
neurons_effective <- load_neurons_effective('effective_neurons.csv')
areas <- unique(neurons$area)


# Load phase connections
baseline <- load_phase_connections(dataset, 'baseline', neurons)
study <- load_phase_connections(dataset, 'study', neurons)
exp_old <- load_phase_connections(dataset, 'exp_old', neurons)
exp_new <- load_phase_connections(dataset, 'exp_new', neurons)


#
library(igraph)
plot_circ(neurons_ordered_circ, areas, baseline, 'baseline', dataset, 10, 5, 0.4)

