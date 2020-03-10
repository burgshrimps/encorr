install.packages("igraph") 
install.packages("network") 
install.packages("sna")
install.packages("ggraph")
install.packages("visNetwork")
install.packages("threejs")
install.packages("networkD3")
install.packages("ndtv")
install.packages('circlize')


animal <- 'LE84'
date <- '20190712'
shift <- 135
v_size <- 3
v_label_size <- 0.2
dataset_name <- paste(animal, date, sep='_')
wd <- paste('/Users/burgshrimps/Documents/lin/analysis/XCORR/', animal, '/', date, '/conn_stat', sep='')
setwd(wd)

library(plyr)
neurons <- read.csv(paste(dataset_name, '_neurons.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
neurons <- neurons[neurons$area != 'cortex' & neurons$area != 'x', ]
neurons[neurons == 'pCA3x'] <- 'pCA3'
neurons$area <- factor(neurons$area, levels = c('sub', 'dCA1', 'CA1b', 'pCA1b', 'pCA1', 'dCA3', 'pCA3'))
neurons$area.num <- as.numeric(neurons$area)
neurons <- neurons[order(neurons$area.num), ]
areas <- unique(neurons$area)


baseline <- read.csv(paste(dataset_name, '_baseline.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
baseline$type.num <- as.numeric(factor(baseline$type))
baseline <- baseline[baseline$ref %in% neurons$id & baseline$tar %in% neurons$id, ]

study <- read.csv(paste(dataset_name, '_study.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
study$type.num <- as.numeric(factor(study$type))
study <- study[study$ref %in% neurons$id & study$tar %in% neurons$id, ]

exp_old <- read.csv(paste(dataset_name, '_exp_old.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
exp_old$type.num <- as.numeric(factor(exp_old$type))
exp_old <- exp_old[exp_old$ref %in% neurons$id & exp_old$tar %in% neurons$id, ]

exp_new <- read.csv(paste(dataset_name, '_exp_new.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
exp_new$type.num <- as.numeric(factor(exp_new$type))
exp_new <- exp_new[exp_new$ref %in% neurons$id & exp_new$tar %in% neurons$id, ]


library('igraph')
plot_circ <- function(neurons, phase, title, dataset_name, shift, v_size, v_label_size) {
  net <- graph_from_data_frame(d=phase, vertices=neurons, directed=FALSE)
  l <- layout_in_circle(net)
  end <- length(l[,1])+shift-1
  l <- rbind(l, l)[shift:end, ]
  V(net)$size <- v_size
  
  colrs <- c('yellow', 'red', 'greenyellow', 'greenyellow', 'green', 'cyan', 'orange')
  colors <- c('red', 'blue')
  
  if (is.null(E(net)$weight)){
    E(net)$width <- 0
  } else {
    E(net)$width <- abs(E(net)$weight)
  }
  
  E(net)$color <- colors[E(net)$type.num]
  V(net)$color <- colrs[V(net)$area.num]
  pdf(paste(dataset_name, '_network_',title, '.pdf', sep=''))
  plot(net, layout=l, vertex.label.cex=v_label_size)
  title(paste(dataset_name, ' ', title, sep=''))
  
  legend("topleft", legend=areas, col=colrs[unique(neurons$area.num)], cex=0.7, pch=21, pt.bg=colrs[unique(neurons$area.num)])

  dev.off()
}

plot_circ(neurons, baseline, 'baseline', dataset_name, shift, v_size, v_label_size)
plot_circ(neurons, study, 'study', dataset_name, shift, v_size, v_label_size)
plot_circ(neurons, exp_old, 'exp_old', dataset_name, shift, v_size, v_label_size)
plot_circ(neurons, exp_new, 'exp_new', dataset_name, shift, v_size, v_label_size)


