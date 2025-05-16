
library(rgdal)
library(raster)
library(sp)
library(sf)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggmap)
library(ggplot2) 
library(viridis) # color palettes
library(scales)
library(rasterVis)
library(lattice)
library(leaflet) # for interactive maps
library(mapview) # for interactive maps
library(shiny)   # for web applications
library(foreach)
library(geosphere)
library(colorspace)
library(readr)
library(terra)
library(readxl)
library(dendsort)
library(ClusterR)
library(NbClust)
library(factoextra)
library(dendextend)
library(RColorBrewer)

theme_set(theme_bw())


list_ras <- list.files("./", pattern = "pred.tif$", full.names = TRUE)
list_ras

ras <- stack(list_ras)
ras

df1 <-na.omit(as.data.frame(ras,xy=TRUE))

### Find the optimal cluster 
### Initial centroids
set.seed(1)

k1000 <- KMeans_rcpp(df1[,3:5], clusters = 1000, num_init = 20, max_iters = 100,
                    fuzzy = TRUE, initializer = 'kmeans++', verbose = T)

cen1000 <- data.frame(k1000$centroids)

#### final centroids

NClust1 <- NbClust::NbClust(data = cen1000, 
                              distance = "euclidean",
                              min.nc = 8, max.nc = 50, 
                              method="ward.D2")

table(NClust1$Best.nc[1,])

##get the optimal Ncluter

Ncluter <- as.numeric(names(sort(table(NClust1$Best.nc[1,]),decreasing = T)[1]))

###plot

b <- data.frame(NClust1$Best.nc[1,])
colnames(b) <- "Cluster"

ggplot( data = b,aes(x=b$Cluster)) +
  geom_histogram( binwidth=1, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Fungi") +
  xlab("Cluster")+ylab("Frequency")+
  theme_bw() +
  theme(plot.title = element_text(size=15))

### the whole map

set.seed(1)

k.N <- KMeans_rcpp(df1[,3:5], clusters = Ncluter, num_init = 20, max_iters = 1000,
                             fuzzy = TRUE, initializer = 'kmeans++', verbose = T)

cen1 <- data.frame(k.N$centroids)

### split and assign cluster
sp <- sample(1:20, nrow(df1), replace = T)
df1$sp <- sp

for(i in 1:20){
  print(i)
  df1.i <- df1[df1$sp == i,]
  pr <- ClusterR::predict_MBatchKMeans(data=as.matrix(df1.i[, 3:5]), CENTROIDS = as.matrix(cen1), fuzzy = FALSE)
  k <- as.vector(pr)
  df1[df1$sp==i,]$Cluster <- k
  gc()
}


# Hierarchical clustering of pedogenons and color legend ------------
pal <- c("#ab51e3","#ce4257","#9c6644","#ff9770","#80b918","#7678ed",
         "#6BAED6","#33A02C","#2196f3","#04471c","#00cc66","#97a97c",
         "#e9f5db","#f4c1e8","#979dac")


hc.k1 <- hclust(dist(cen1), method="ward.D2")

viz.k1 <- hc.k1 %>%  as.dendrogram(.) %>% 
set(.,"branches_lwd", 1.5) #  %>% 

labels_colors(viz.k11) <- pal

plot(viz.k11,,horiz = FALSE)

branch_labels <- as.vector(labels(as.dendrogram(hc.k11)))
branch_labels

########## plot the clustered map
Fg3 <- df1[,c(1,2,6)]
Fg3 <- rasterFromXYZ(Fg3) 
crs(Fg3)<- crs(ras)
Fg3 


p2 <- ggplot()+
  geom_spatraster(data = rast(Fg3),alpha = 0.8) +
  scale_fill_manual(values = pal,na.value = "white") +
  ggtitle("Fungal Distribution Pattern")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "right") 
p2

writeRaster(Fg3,"Fg2040.tif",format="GTiff",datatype ="FLT4S",overwrite=TRUE)

