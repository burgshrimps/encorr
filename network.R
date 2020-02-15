install.packages("igraph") 
install.packages("network") 
install.packages("sna")
install.packages("ggraph")
install.packages("visNetwork")
install.packages("threejs")
install.packages("networkD3")
install.packages("ndtv")
install.packages('circlize')

setwd('/Users/burgshrimps/Documents/lin/analysis/XCORR/LE87/20190520/conn_stat')
dataset_name <- 'LE87_20190520'

neurons <- read.csv(paste(dataset_name, '_neurons.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
neurons$area.num <- as.numeric(factor(neurons$area))

baseline <- read.csv(paste(dataset_name, '_baseline.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
baseline$type.num <- as.numeric(factor(baseline$type))
study <- read.csv(paste(dataset_name, '_study.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
study$type.num <- as.numeric(factor(study$type))
exp_old <- read.csv(paste(dataset_name, '_exp_old.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
exp_old$type.num <- as.numeric(factor(exp_old$type))
exp_new <- read.csv(paste(dataset_name, '_exp_new.csv', sep=''), header=TRUE, sep=',', dec='.', check.names=FALSE)
exp_new$type.num <- as.numeric(factor(exp_new$type))

# pCA3 : darkgoldenrod1
# dCA3 : darkgoldenrod (dark)
# pCA1 : cyan
# dCA1 : cyan4 (dark)

library('circlize')
plot_chord <- function(neurons, phase, title, dataset_name) {
  circos.clear()
  phase$ref.area <- neurons[match(phase$ref, neurons$id),4]
  phase$tar.area <- neurons[match(phase$tar, neurons$id),4]
  phase$ref.tet <- neurons[match(phase$ref, neurons$id),2]
  phase$tar.tet <- neurons[match(phase$tar, neurons$id),2]
  pdf(paste(dataset_name, '_chord_',title, '.pdf', sep=''))
  ### LE46 ###
  #grid.col <- c('1'='darkgoldenrod1', '2'='cyan4', '3'='cyan4', '4'='cyan4', '7'='cyan', '8'='cyan')
  
  ### LE82 ###
  #grid.col <- c('1'='cyan4', '2'='cyan4', '3'='cyan4', '4'='darkgoldenrod1', '5'='darkgoldenrod1', '6'='cyan4', '7'='cyan4', '9'='darkgoldenrod', '10'='darkgoldenrod', '11'='cyan', '12'='darkgoldenrod', '13'='cyan4', '14'='darkgoldenrod', '15'='darkgoldenrod1', '16'='darkgoldenrod1')
  
  ### LE83 ###
  #grid.col <- c('2'='cyan4', '3'='cyan4', '4'='darkgoldenrod1', '5'='cyan4','8'='cyan4', '9'='darkgoldenrod', '10'='cyan', '13'='darkgoldenrod1', '15'='darkgoldenrod1', '16'='darkgoldenrod1')
  
  ### LE84 ###
  #grid.col <- c('1'='cyan4', '2'='cyan4', '3'='cyan4', '4'='darkgoldenrod1', '5'='cyan4', '6'='cyan4', '7'='cyan4', '8'='cyan4', '9'='darkgoldenrod', '10'='darkgoldenrod', '11'='darkgoldenrod1', '12'='cyan', '13'='darkgoldenrod1', '14'='darkgoldenrod1', '15'='cyan', '16'='darkgoldenrod1')
  
  ### LE87 ###
  grid.col <- c('3'='cyan4', '4'='darkgoldenrod1', '7'='cyan4', '8'='cyan4', '9'='darkgoldenrod', '10'='darkgoldenrod', '11'='cyan', '12'='darkgoldenrod', '13'='cyan', '14'='darkgoldenrod', '15'='darkgoldenrod1', '16'='darkgoldenrod1')
  
  chordDiagram(table(phase[,8:9]), grid.col = grid.col, col = 'grey', link.border = 'darkgrey', preAllocateTracks = list(track.height = uh(4, "mm"),track.margin = c(uh(4, "mm"), 0)), annotationTrack = c("grid", "axis"))
  title(paste(dataset_name, ' ', title, sep=''))
  circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, niceFacing = TRUE)
  }, bg.border = NA)
  
  ### LE46 ###
  #highlight.sector(c('2', '3', '4'), track.index = 1, col = "cyan4", text = "dCA1", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('1'), track.index = 1, col = "darkgoldenrod1", text = "pCA3", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('7', '8'), track.index = 1, col = "cyan", text = "pCA1", cex = 0.6, text.col = "black", niceFacing = TRUE)
  
  ### LE82 ###
  #highlight.sector(c('2', '3', '13'), track.index = 1, col = "cyan4", text = "dCA1", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('4', '5'), track.index = 1, col = "darkgoldenrod1", text = "pCA3", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('9', '10'), track.index = 1, col = "darkgoldenrod", text = "dCA3", cex = 0.6, text.col = "black", niceFacing = TRUE)
  
  ### LE83 ###
  #highlight.sector(c('2', '3', '5', '8'), track.index = 1, col = "cyan4", text = "dCA1", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('4', '13', '15', '16'), track.index = 1, col = "darkgoldenrod1", text = "pCA3", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('9'), track.index = 1, col = "darkgoldenrod", text = "dCA3", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('10'), track.index = 1, col = "cyan", text = "pCA1", cex = 0.6, text.col = "black", niceFacing = TRUE)
  
  ### LE84 ###
  #highlight.sector(c('1', '2', '7', '8'), track.index = 1, col = "cyan4", text = "dCA1", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('4', '11', '13', '16'), track.index = 1, col = "darkgoldenrod1", text = "pCA3", cex = 0.6, text.col = "black", niceFacing = TRUE)
  #highlight.sector(c('12', '15'), track.index = 1, col = "cyan", text = "pCA1", cex = 0.6, text.col = "black", niceFacing = TRUE)
  
  ### LE87 ###
  highlight.sector(c('7'), track.index = 1, col = "cyan4", text = "dCA1", cex = 0.6, text.col = "black", niceFacing = TRUE)
  highlight.sector(c('4', '16'), track.index = 1, col = "darkgoldenrod1", text = "pCA3", cex = 0.6, text.col = "black", niceFacing = TRUE)
  highlight.sector(c('9', '10'), track.index = 1, col = "darkgoldenrod", text = "dCA3", cex = 0.6, text.col = "black", niceFacing = TRUE)
  
  dev.off()
}

plot_chord(neurons, baseline, 'baseline', dataset_name)
plot_chord(neurons, study, 'study', dataset_name)
plot_chord(neurons, exp_old, 'exp_old', dataset_name)
plot_chord(neurons, exp_new, 'exp_new', dataset_name)

library('igraph')
plot_circ <- function(neurons, phase, title, dataset_name) {
  net <- graph_from_data_frame(d=phase, vertices=neurons, directed=FALSE)
  l <- layout_in_circle(net)
  V(net)$size <- 3
  ### LE46 ###
  #colrs <- c("cyan4", "cyan", "darkgoldenrod1")
  
  ### LE82, LE83 ###
  colrs <- c("cyan4", "darkgoldenrod", "cyan", "darkgoldenrod1")
  
  colors <- c('red', 'blue')
  E(net)$width <- abs(E(net)$weight)
  E(net)$color <- colors[E(net)$type.num]
  V(net)$color <- colrs[V(net)$area.num]
  pdf(paste(dataset_name, '_network_',title, '.pdf', sep=''))
  plot(net, layout=l, vertex.label.cex=0.2)
  title(paste(dataset_name, ' ', title, sep=''))
  dev.off()
}

plot_circ(neurons, baseline, 'baseline', dataset_name)
plot_circ(neurons, study, 'study', dataset_name)
plot_circ(neurons, exp_old, 'exp_old', dataset_name)
plot_circ(neurons, exp_new, 'exp_new', dataset_name)


