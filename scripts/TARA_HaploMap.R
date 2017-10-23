#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(maps)
library(plyr)
library(gridExtra)

##########################################################
# Set this directory to where your TARA-*.txt files are :)
##########################################################
setwd("/Users/u1375038/repos/Ebame3/2017_DESMAN")

gen_blank_world_map <- function(df) {
  world_map <- map_data("world")
  
  min_lat <- -60
  max_lat <- 70
  min_lon <- -155
  max_lon <- 75
  
  p <- ggplot(data=world_map,aes(x=long, y=lat, group=group))
  p <- p + geom_polygon(fill = '#777777') + coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))  
  p <- p + scale_colour_manual(values = c("Gr01" = "brown1", "Gr02" = "black", "Gr01_A" = "darkorchid", "Gr01_B" = "TAN", "Gr01_C" = "seagreen3", "Gr02_A" = "brown1", "Gr02_B" = "gold2", "Gr02_C" = "cadetblue3", "Gr09" = "plum", "Gr10" = "moccasin", "Other" = "gray88"))
  return(p)
}

add_TARA_sampling_locations <- function(plot_object, tara_samples, size = 5, alpha=0.2, labels = TRUE){
  plot_object <- plot_object + geom_jitter(data = tara_samples,
                                           position=position_jitter(width=0, height=0),
                                           aes(x=Mean_Long, y=Mean_Lat, group=cluster_id, size=TARA_MED_MAG_00110, color=cluster_id), alpha=alpha, fill="gray")
  plot_object <- plot_object + scale_size(range=c(2, 15))

  if (labels){
    plot_object <- plot_object + geom_text(data = tara_samples,
                                           aes(x=Mean_Long, y=Mean_Lat, group=ocean, label=label),
                                           size=0, hjust = 0, nudge_x = 1, fontface = "bold")
  }
  
  return(plot_object)
}

clean_map <- function(plot_object){
  plot_object <- plot_object +
    xlab(NULL) +
    ylab(NULL) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(colour = "black"))
  return(plot_object)
}


# GET TARA samples data frame
tara_samples <- read.table(file = 'TARA-samples.txt', header = TRUE , sep = "\t" , quote = "")
#tara_samples[tara_samples$Mean_Depth > 0, ]$Mean_Lat <- tara_samples[tara_samples$Mean_Depth > 10, ]$Mean_Lat - 0

# GET TARA CLUSTERS
tara_clusters <- read.table(file = 'TARA-clusters.txt', header = TRUE , sep = "\t" , quote = "")

# ADD LABELS and CLUSTER IDs TO TARA SAMPLES DATA FRAME
tara_samples$label <- "UNKNOWN LABEL"
tara_samples$cluster_id <- "UNKNOWN CLUSTER"
for (sample in levels(tara_samples$samples)){
  single_sample <- tara_samples[tara_samples$samples == sample, ]
  sample_label <- paste(sapply(strsplit(sample, '_'), '[', 1), " ", sapply(strsplit(sample, '_'), '[', 2), ', ', single_sample$Mean_Temperature, "C", sep="")
  tara_samples[tara_samples$samples == sample, ]$label <- sample_label
  cluster_id <- tara_clusters[tara_clusters$samples == sample, ]$Groups
  if (length(cluster_id)){
    tara_samples[tara_samples$samples == sample, ]$cluster_id <- as.character(cluster_id)
  }
}

# SUBSAMPLE TO max 10m DEPTH
tara_samples <- tara_samples[tara_samples$Mean_Depth < 17, ]
tara_samples$samples <- factor(tara_samples$samples)

# GENERATE THE PLOT OBJECT
world_map <- gen_blank_world_map(tara_samples)
world_map <- add_TARA_sampling_locations(world_map, tara_samples, alpha=1)
world_map <- clean_map(world_map)

# SAVE IT
pdf('tara_world.pdf', width=16, height=16)
print(world_map)
dev.off()
