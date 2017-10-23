#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(maps)
library(plyr)
library(gridExtra)
library(getopt)

spec = matrix(c('verbose','v',0,"logical",'help','h',0,"logical",'gammafile','g',1,"character",'samplefile','s',1,"character",'clusterfile','c',1,"character"),byrow=TRUE,ncol=4)

opt=getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE)); 
  q(status=1);
}


gen_blank_world_map <- function(df) {
  world_map <- map_data("world")
  
  min_lat <- -60
  max_lat <- 45
  min_lon <- -155
  max_lon <- 75
  
  p <- ggplot(data=world_map,aes(x=long, y=lat, group=group))
  p <- p + geom_polygon(fill = '#777777') + coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))  
  return(p)
}

add_TARA_sampling_locations <- function(plot_object, tara_samples, selectHaplo, size = 5, alpha=0.2, labels = TRUE){
  print(selectHaplo)
  plot_object <- plot_object + geom_jitter(data = tara_samples,
                                           position=position_jitter(width=0, height=0),
                                           aes_string(x="Mean_Long", y="Mean_Lat", group="cluster_id", size=selectHaplo,alpha=selectHaplo), colour="red")
                                           #aes(x=Mean_Long, y=Mean_Lat, group=cluster_id, size=selectHaplo, alpha=selectHaplo), colour="red")
  plot_object <- plot_object + scale_size(range=c(0.5, 16), limits=c(0.0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1))

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



# GET TARA CLUSTERS
tara_clusters <- read.table(file = opt$clusterfile, header = TRUE , sep = "\t" , quote = "")
tara_samples <- read.table(file = opt$samplefile, header = TRUE , sep = "\t" , quote = "")

gamma <- read.csv(opt$gammafile,header=TRUE,row.names=1)

colnames(gamma) <- gsub("X","H",colnames(gamma))
rownames(gamma) <- gsub("m","M",rownames(gamma))
rownames(gamma) <- gsub("_5M","_05M",rownames(gamma))
rownames(tara_samples) <- tara_samples$samples 
gamma$samples <- rownames(gamma)

tara_samples <- merge(tara_samples, gamma, by="samples",all.x = TRUE)

rownames(tara_samples) <- tara_samples$samples
gamma$samples <- NULL

G <- ncol(gamma)
T <- ncol(tara_samples)
S <- T-G
tara_samples[, S:T][is.na(tara_samples[, S:T])] <- 0


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

for (H in 0:G-1){
  hid <- sprintf("H%d", H)
  filename <- sprintf("WorldMap_H%d.pdf", H)
  print(filename)

  # GENERATE THE PLOT OBJECT
  world_map <- gen_blank_world_map(tara_samples)
  world_map <- add_TARA_sampling_locations(world_map, tara_samples, selectHaplo=hid,alpha=1)
  world_map <- clean_map(world_map)
  
  # SAVE IT
  pdf(filename, width=12, height=5)
  print(world_map + guides(alpha=FALSE))
  dev.off()  
}

