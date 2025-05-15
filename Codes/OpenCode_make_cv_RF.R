
# function to fit a model using ranger
# Arguments:  
# # data: spatialpointdataframe with columns where there is property and raster covariates values at points
# # property: name of the soil property
# # covars: names of the raster layers to be used for model fitting
# # raster: one raster to be used as background if spatial.cv is TRUE
# # spatial.cv: TRUE or FALSE. If TRUE, a spatial cross-validationis made using raster argument as background. 
# # # if otherwise (spatial.cv = F) a random split is made (each point is used only once for validation)
# # nfolds: number of folds for the data split
## tunning: if TRUE, a parameter tuning with grid search is made

make_cv_RF <- function(data, property, covars, raster, spatial.cv, nfolds, tunning){
  
  #k <- nfolds
  library(sf)
  # remove rows with NA
  #data <- na.omit(st_as_sf(data), margin = 1)
  
  # create a data.frame to store the prediction of each fold (record)
  testTable <- data.frame(obs = as.data.frame(data)[property], pred = NA)
  testTable <- cbind(testTable, as.data.frame(data))
  # make the formula
  form <- as.formula(paste(property, paste(covars, collapse=" + "), sep=" ~ "))
  
  if(spatial.cv == T){
    
    # make spatial blocks
    sb <- spatialBlock(speciesData = data,
                       species = property,
                       rasterLayer = raster,
                       theRange = 200000, # size of the blocks
                       k = nfolds,
                       selection = "random",
                       iteration = 10, # find evenly dispersed folds
                       biomod2Format = TRUE,
                       xOffset = 0, # shift the blocks horizontally
                       yOffset = 0, 
                       verbose = FALSE, 
                       progress = FALSE) 
    
    # extract the fold indices from buffering  
    # created in the previous section
    # the folds (list) works for all three blocking strategies
    folds <- sb$folds
    
    for(i in 1:nfolds){
      
      # extracting the training and testing indices
      # this way works with folds list (but not foldID)
      trainSet <- unlist(folds[[i]][1]) # training set indices
      testSet <- unlist(folds[[i]][2]) # testing set indices
      
      dat.train <- as.data.frame(data[trainSet, ][c(property, covars)])
      dat.train <- na.omit(dat.train)
      
      #################### 
      ## model parameter optimization
      library(caret)
      library(ranger)
      
      if(tunning == TRUE){
        grid <-  expand.grid(mtry = c(2, 4, 6), 
                             splitrule = 'variance', 
                             min.node.size = c(2 ,5, 10))
        
        fitControl <- trainControl(method = "CV",
                                   number = 5,
                                   verboseIter = TRUE)
        
        fit = train(
          x = dat.train[ , names(dat.train) != property],
          y = dat.train[ , names(dat.train) == property],
          method = 'ranger',
          tuneGrid = grid,
          num.trees = 250,
          trControl = fitControl, 
          quantreg = TRUE)
        
        # prediction using best parameters
        rf <- ranger(form, data = dat.train, 
                     mtry = fit$bestTune$mtry,
                     min.node.size = fit$bestTune$min.node.size,
                     splitrule = 'variance',
                     num.trees = 250, 
                     quantreg = TRUE) # model fitting on training set
        
        dat.test <- as.data.frame(data[testSet, ][covars])
        testTable$pred[testSet] <- predict(rf, data = dat.test)$predictions # predict the test set
        
      }else{
        
        # prediction using best parameters
        rf <- ranger(form, data = dat.train, 
                     num.trees = 250, 
                     quantreg = TRUE) # model fitting on training set
        
        dat.test <- as.data.frame(data[testSet, ][covars])
        testTable$pred[testSet] <- predict(rf, data = dat.test)$predictions # predict the test set
      }
    }
    
    # best fit using all data
    rf_all_data <- ranger(form, data = as.data.frame(data), 
                          num.trees = 250, 
                          importance = 'impurity')
    
    return(list(rf_all_data, testTable))
  }
  else{
    
    
    # samples that will be used for the validation
    sp <- sample(1:nfolds, nrow(data), replace = T)
    sp <- sample(sp)
    data$sp <- sp
    
    
    for(i in 1:nfolds){
      # extracting the training and testing indices
      # this way works with folds list (but not foldID)
      
      dat.train <- as.data.frame(data[data$sp != i,][c(property, covars)])
      
      
      
      #################### 
      ## model parameter optimization
      library(caret)
      library(ranger)
      
      if(tunning == TRUE){
        grid <-  expand.grid(mtry = c(2, 4, 6),     
                             splitrule = 'variance', 
                             min.node.size = c(2 ,5, 10))
        
        fitControl <- trainControl(method = "CV",
                                   number = 5,
                                   verboseIter = TRUE)
        
        fit = caret::train(
          x = dat.train[ , names(dat.train) != property],
          y = dat.train[ , names(dat.train) == property],
          method = 'ranger',
          tuneGrid = grid,
          num.trees = 250,
          trControl = fitControl, 
          quantreg = TRUE)
        
        # prediction using best parameters
        rf <- ranger(form, data = dat.train, 
                     mtry = fit$bestTune$mtry,
                     min.node.size = fit$bestTune$min.node.size,
                     num.trees = 250, 
                     quantreg = TRUE) # model fitting on training set
        
        testSet <- as.data.frame(data[data$sp == i,][covars])
        testTable$pred[data$sp == i] <- predict(rf, data = testSet)$predictions # predict the test set
      }else{
        
        # prediction using best parameters
        rf <- ranger(form, data = dat.train, 
                     num.trees = 250, 
                     quantreg = TRUE) # model fitting on training set
        
        testSet <- as.data.frame(data[data$sp == i,][covars])
        testTable$pred[data$sp == i] <- predict(rf, data = testSet)$predictions # predict the test set
      }
    }
    
    # best fit using all data
    rf_all_data <- ranger(form, data = as.data.frame(data), 
                          num.trees = 250, 
                          importance = 'impurity', 
                          quantreg = TRUE)
    testTable$pred.oob <- predict(rf_all_data, as.data.frame(data))$predictions
    testTable2 <- cbind(testTable, as.data.frame(data))
    
    return(list(rf_all_data, testTable2))
  }
  
}
