
library(raster)
library(rgdal)
library(ranger)
library(parallel)
library(sp)

####load input layers
list_ras <- list.files(".", pattern=".tif$", full.names = TRUE)
list_ras
covariates.stack <- raster::stack(list_ras)

# make blocks

bs <- blockSize(covariates.stack,chunksize = 10000000 , minblocks = 50)
bs

#################### Prediction
###Load data 
df1 <- read_csv("./biodiversity.data.csv")

### Map generation

components <-  colnames(df1)[-1]

# Helper function to predict and save tiles
process_tile <- function(i, model, compNumb) {
  tile_temp <- readAll(crop(covariates.stack, extent(covariates.stack, bs$row[[i]], bs$row[[i]] + bs$nrows[[i]], 1, ncol(covariates.stack))))
  
  pred <- raster::predict(tile_temp, model, fun = function(model, ...) predict(model, ...)$predictions)
  predvar05 <- raster::predict(tile_temp, model, fun = function(model, ...) predict(model, type = 'quantiles', quantiles = 0.05, ...)$predictions[,1])
  predvar95 <- raster::predict(tile_temp, model, fun = function(model, ...) predict(model, type = 'quantiles', quantiles = 0.95, ...)$predictions[,1])
  
  saveRDS(pred, file = paste0('./pred', compNumb, '/pred/', i, '.rds'))
  saveRDS(predvar05, file = paste0('./pred', compNumb, '/predvar05/', i, '.rds'))
  saveRDS(predvar95, file = paste0('./pred', compNumb, '/predvar95/', i, '.rds'))
}

# Helper function to merge tiles and write raster
merge_tiles <- function(folder, compNumb, suffix) {
  files <- list.files(folder, full.names = TRUE)
  rast.list <- lapply(files, readRDS)
  
  names(rast.list)[1:2] <- c('x', 'y')
  rast.list$fun <- mean
  rast.list$na.rm <- TRUE
  
  y_pred <- do.call(mosaic, rast.list)
  plot(y_pred)
  
  writeRaster(y_pred, filename = paste0('./', compNumb, suffix, '.tif'), overwrite = TRUE)
}

# Main loop for components
for (compNumb in components) {
  # Load model
  x <- load(paste0('./', compNumb, '_model_evaluationCV.Rdata'))
  model <- get(x)[[1]]
  rm(x)
  
  # Create output directories
  dir.create(paste0('./pred', compNumb), showWarnings = FALSE)
  dir.create(paste0('./pred', compNumb, '/pred'), showWarnings = FALSE)
  dir.create(paste0('./pred', compNumb, '/predvar05'), showWarnings = FALSE)
  dir.create(paste0('./pred', compNumb, '/predvar95'), showWarnings = FALSE)
  
  # Process tiles
  llply(seq_len(bs$n), function(i) process_tile(i, model, compNumb))
  
  # Merge and save predicted maps
  merge_tiles(paste0('./pred', compNumb, '/pred'), compNumb, '_pred')
  merge_tiles(paste0('./pred', compNumb, '/predvar05'), compNumb, '_predvar05p')
  merge_tiles(paste0('./pred', compNumb, '/predvar95'), compNumb, '_predvar95p')
  
  gc()
}