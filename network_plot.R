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
  
  neurons_effective <- read.csv(filename, header=TRUE, sep=',', check.names=FALSE)
  
  return(neurons_effective)
}


reorder_neurons_for_circ <- function(neurons) {
  ### Change order of neurons in regards to the order they appear in the network plot based on brain area ###
  
  neurons$area <- factor(neurons$area, levels = c('sub', 'dCA1', 'CA1b', 'pCA1b', 'pCA1', 'dCA3', 'pCA3'))
  neurons$area.num <- as.numeric(neurons$area)
  neurons <- neurons[order(neurons$area.num), ]
  
  return(neurons)
}


filter_effective_neurons <- function(neurons, effective_neurons_names, animal) {
  ### Filter neuron vertix list and only return those that are marked as effective ###
  
  neurons <- neurons[neurons$id %in% effective_neurons_names[animal][,], ]
  return(neurons)
}

filter_effective_connections <- function(neurons_effective, phase, logop) {
  ### Filters connections such that either both, ref and tar id of connection, have to be in effective neurons list or ###
  ### only one of them. Specified by logical operator parameter 'logop', either & or | ###
  
  if (logop == '&') {
    phase_effective <- phase[phase$ref %in% neurons_effective$id & phase$tar %in% neurons_effective$id, ]
  } else if (logop == '|') {
    phase_effective <- phase[phase$ref %in% neurons_effective$id | phase$tar %in% neurons_effective$id, ]
  }
  
  return(phase_effective)
}


load_phase_connections <- function(dataset, phase, neurons) {
  ### Loads connections between neurons measered by cross-correlation from CSV ###
  
  phase_connections <- read.csv(paste(dataset, '_', phase, '.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
  phase_connections$type.num <- as.numeric(factor(phase_connections$type))
  
  # Make sure to only include connections between neurons which are actually in the dataset, e.g. not between neurons labeled with 'x'
  phase_connections <- phase_connections[phase_connections$ref %in% neurons$id & phase_connections$tar %in% neurons$id, ]
  
  return(phase_connections)
}


compute_coords_around_circle <- function(x0, y0, r, n) {
  ### Computes coordinates of n points evenly distributed around a circle with center (x0, y0) and radius r. ###
  
  t <- 2*pi*1:n/n
  coords <- t(sapply(t, function(t)c(x0 + r * cos(t), y0 + r + sin(t))))
  
  return(coords)
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


compute_anatomical_layout <- function(net, neurons_effective, center, r) {
  ### Computes the coordinates of each vertix based on their corresponding brain area ###
  
  l <- layout_randomly(net)
  l[neurons_effective$area == 'dCA1', ] <- compute_coords_around_circle(center, center, r, length(l[neurons_effective$area == 'dCA1', 1]))
  l[neurons_effective$area == 'dCA3', ] <- compute_coords_around_circle(-center, -center, r, length(l[neurons_effective$area == 'dCA3', 1]))
  l[neurons_effective$area == 'pCA3', ] <- compute_coords_around_circle(center, -center, r, length(l[neurons_effective$area == 'pCA3', 1]))
  l[neurons_effective$area == 'pCA1', ] <- compute_coords_around_circle(-center, center, r, length(l[neurons_effective$area == 'pCA1', 1]))
  l[neurons_effective$area == 'sub', ] <- compute_coords_around_circle(center, 0, r, length(l[neurons_effective$area == 'sub', 1]))
  l[neurons_effective$area == 'CA1b', ] <- compute_coords_around_circle(0, center, r, length(l[neurons_effective$area == 'CA1b', 1]))
  l[neurons_effective$area == 'pCA1b', ] <- compute_coords_around_circle(-center+1, center-1, r, length(l[neurons_effective$area == 'pCA1b', 1]))
  
  return(l)
}


plot_anatomical <- function(neurons_effective, phase_effective, center, radius, dataset, title) {
  ### Plot effective neurons and connections between them in a layout which roughly represents their anatomical layout in hippocampus ###
  
  net <- graph_from_data_frame(d=phase_effective, vertices=neurons_effective, directed=FALSE)
  l <- compute_anatomical_layout(net, neurons_effective, center, radius)
  
  # Set color scheme
  v_color <- c('yellow', 'red', 'greenyellow', 'greenyellow', 'green', 'cyan', 'orange') 
  e_color <- c('red', 'blue')
  V(net)$color <- v_color[V(net)$area.num]
  E(net)$color <- e_color[E(net)$type.num]
  
  # Plot and save as PDF
  pdf(paste(dataset, '_network_', title, '_effective_anatomical', '.pdf', sep=''))
  plot(net, layout=l)
  title(paste(dataset, ' ', title, sep=''))
  legend("bottom", legend=areas_effective, col=v_color[unique(neurons_effective$area.num)], pch=21, pt.bg=v_color[unique(neurons_effective$area.num)])
  dev.off()
}


# Set dataset parameters
animal <- 'LE82'
date <- '20190424'
dir <- '/Users/burgshrimps/Documents/lin/analysis/XCORR/'
dataset <- paste(animal, date, sep='_')
wd <- paste(dir, animal, '/', date, '/conn_stat', sep='')
setwd(wd)


# Set parameters for circular plot
shift <- 10
v_size <- 5
v_label_size <- 0.4


# Set parameters for anatomical plot
center <- 5
radius <- 1


# Load neurons
library(plyr)
neurons <- load_neurons(dataset)
neurons_ordered_circ <- reorder_neurons_for_circ(neurons)
effective_neurons_names <- load_neurons_effective(paste(dir, 'effective_neurons.csv', sep=''))
areas <- unique(neurons$area)


# Load phase connections
baseline <- load_phase_connections(dataset, 'baseline', neurons)
study <- load_phase_connections(dataset, 'study', neurons)
exp_old <- load_phase_connections(dataset, 'exp_old', neurons)
exp_new <- load_phase_connections(dataset, 'exp_new', neurons)


# Plot all neurons and all connections in circular layout
library(igraph)
plot_circ(neurons_ordered_circ, areas, baseline, 'baseline', dataset, shift, v_size, v_label_size)
plot_circ(neurons_ordered_circ, areas, study, 'study', dataset, shift, v_size, v_label_size)
plot_circ(neurons_ordered_circ, areas, exp_old, 'exp_old', dataset, shift, v_size, v_label_size)
plot_circ(neurons_ordered_circ, areas, exp_new, 'exp_new', dataset, shift, v_size, v_label_size)


# Filter only effective neurons and connections
neurons_effective <- filter_effective_neurons(neurons_ordered_circ, effective_neurons_names, animal)
areas_effective <- unique(neurons_effective$area)
baseline_effective <- filter_effective_connections(neurons_effective, baseline, '&')
study_effective <- filter_effective_connections(neurons_effective, study, '&')
exp_old_effective <- filter_effective_connections(neurons_effective, exp_old, '&')
exp_new_effective <- filter_effective_connections(neurons_effective, exp_new, '&')


# Plot effective neurons and connections in anatomical layout
plot_anatomical(neurons_effective, baseline_effective, center, radius, dataset, 'baseline')
plot_anatomical(neurons_effective, study_effective, center, radius, dataset, 'study')
plot_anatomical(neurons_effective, exp_old_effective, center, radius, dataset, 'exp_old')
plot_anatomical(neurons_effective, exp_new_effective, center, radius, dataset, 'exp_new')




