
library(raster)
library(sp)
library(sf)
library(dplyr)
library(tidyverse)
library(tidyr)
library(colorspace)
library(readr)
library(terra)
library(soilspec)

theme_set(theme_bw())

###Load data 
biodiversity_data <- read_csv("./biodiversity.data.csv")

###load input layers
list_ras <- list.files(".", pattern=".tif$", full.names = TRUE)
list_ras
covariates.stack <- raster::stack(list_ras)

#### extraction
df1 <- biodiversity_data

lonlat <- df1[,c("Longitude","Latitude")]

extr1 <- raster::extract(covariates.stack,lonlat)
soildat <- cbind(df1, extr1[,-1])


######### model assessing


tab.res <- data.frame(matrix(data=NA, nrow = 3, ncol = 5))

rownames(tab.res) <- c("MDS1","MDS2","MDS3")
colnames(tab.res) <- c('ME', 'RMSE', 'r2', 'R2', 'rhoC')

################### Training

source("D:/OneDrive - The University of Sydney (Staff)/0_Projects/1_SoilBetaDiversity_Au/0_codes/make_cv_RF.R")
setwd("D:/OneDrive - The University of Sydney (Staff)/0_Projects/1_SoilBetaDiversity_Au/Predictions/New1/Bac2")


soildat1 <- na.omit(soildat)
sum(is.na(soildat2))

covar.names <- colnames(extr1)
covar.names

var <-  colnames(df1)[-1]
var

##### training models

for (i in 1:6) {
  res_MDS <- make_cv_RF(
    data = soildat2[, c(var[i], covar.names)],
    property = var[i],
    covars = covar.names,
    raster = ras.all,
    spatial.cv = FALSE,
    nfolds = 10,
    tunning = FALSE
  )
  
  # Evaluation and plotting
  eval_metrics <- eval(res_MDS[[2]][,1], res_MDS[[2]]$pred, obj = 'quant')
  plot(res_MDS[[2]][,1], res_MDS[[2]]$pred, main = paste0('MDS', i)); abline(0, 1, col = 'red')
  
  # Store evaluation results
  tab.res[paste0('MDS', i), ] <- eval_metrics[1:5]
  
  # Print variable importance
  print(res_MDS[[1]]$variable.importance[order(res_MDS[[1]]$variable.importance)])
  
  # Save the result model
  res.model <- res_MDS
  save(res.model, file = paste0('./', var[i], '_model_evaluationCV.Rdata'))
}

write.csv(tab.res,file = "tab.res.csv")
