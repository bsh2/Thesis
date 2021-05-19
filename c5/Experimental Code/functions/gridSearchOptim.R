gridStandard <- function(freqIndex){
  
  outdir <- paste0("output/standard/", freqIndex)
  
  lower <- standardConstraints$lower
  upper <- standardConstraints$upper
  
  par1 <- seq(lower[1]+1, upper[1]-1, length.out = gridDimStandard)
  par2<- seq(lower[2]+0.05, upper[2]-0.05, length.out = gridDimStandard)
  par3 <- seq(lower[3]+0.5, upper[3]-0.5, length.out = gridDimStandard)
  par4 <- seq(lower[4]+0.05, upper[4]-0.05, length.out = gridDimStandard)
  par5 <- seq(lower[5]+0.5, upper[5]-0.5, length.out = gridDimStandard)
  
  grid <- expand.grid(par1, par2, par3, par4, par5)
  rm(par1, par2, par3, par4, par5)
  gc()
  
  print("Grid constructed. Dividing into segments")
  division <- length(grid[,1])/gridDivStandard
  
  # CONSTRUCT THE GRIDS
  for(i in 1:gridDivStandard){
    if (i == 1){
      gridTemp <- grid[1:division,]
    } else{
      gridTemp <- grid[((division*(i-1))+1):(division*(i)),]
    }
    saveRDS(gridTemp, file = paste0(outdir,"/grid",i,".Rda"))
    print(paste0("GRID ",i," GENERATED AND SAVED"))
    rm(gridTemp)
    gc()
  }
  
  rm(grid)
  gc()
  
  
  # EVALUATE EACH GRID
  for(i in 1:gridDivStandard){
    print(paste0("Evaluating Grid ",i))
    gridTemp <- readRDS(file = paste0(outdir,"/grid",i,".Rda"))
    cl = makeCluster(c(rep("localhost", detectCores(logical = TRUE))),
                     type = "SOCK")
    clusterExport(cl, c("standardCriterion", "standardIndex", "loads"))
    gridTemp[,6] <- parRapply(cl, gridTemp, standardLoss)
    stopCluster(cl)
    
    # Optimise over the grid
    numCores <- detectCores(logical = TRUE)
    cluster <- makeCluster(numCores)
    registerDoSNOW(cluster)
    standardFitted <- foreach(
      i = 1:dim(gridTemp)[1],
      .combine = "rbind", .verbose = TRUE, .packages = "stats",
      .export = c("standardCriterion", "standardIndex", "loads",
                  "standardPars", "standardCompute", "MAPEfunc",
                  "RMSEfunc", "RSQfunc", "standardConstraints",
                  "standardLoss", "standardFit")) %dopar% {
                    vals <- as.numeric(gridTemp[i,])
                    standardFit(vals)
                  }
    stopCluster(cluster)
    saveRDS(standardFitted, file = paste0(outdir,"/grid",i,".Rda"))
    rm(gridTemp)
    gc()
    print(paste0("GRID: ",i,"/", gridDivStandard," evaluated"))
  }
  
  gc()
  
  print("All grid segments evaluated, compiling results")
  
  for(i in 1:gridDivStandard){
    if (i == 1){
      grid <- readRDS(file = paste0(outdir,"/grid",i,".Rda"))
    } else{
      gridTemp <- readRDS(file = paste0(outdir,"/grid",i,".Rda"))
      grid <- rbind(grid, gridTemp)
      rm(gridTemp)
      gc()
    }
  }

  colnames(grid) <- 
    c("p0Start", "kgStart", "TgStart", "khStart",
      "ThStart", "startRSS", "p0Fitted", "kgFitted", "TgFitted", "khFitted",
      "ThFitted", "optimRSS", "fnCount", "gnCount", "convCode",
      "EigenVal1", "EigenVal2", "EigenVal3", "EigenVal4", "EigenVal5", 
      "testRsq", "testRMSE", "testMAPE", "p0TrueDiff", "kgTrueDiff",
      "TgTrueDiff", "khTrueDiff", "ThTrueDiff")
  
  saveRDS(grid, file = paste0(outdir,"/gridFull.Rda"))
  
  return(grid)
}

gridCalvert <- function(freqIndex){
  
  outdir <- paste0("output/calvert/", freqIndex)
  
  lower <- calvertConstraints$lower
  upper <- calvertConstraints$upper
  
  par1 <- seq(lower[1]+1, upper[1]-1, length.out = gridDimCalvert)
  par2<- seq(lower[2]+0.05, upper[2]-0.05, length.out = gridDimCalvert)
  par3 <- seq(lower[3]+0.5, upper[3]-0.5, length.out = gridDimCalvert)
  par4 <- seq(lower[4]+0.5, upper[4]-0.5, length.out = gridDimCalvert)
  par5 <- seq(lower[5]+0.05, upper[5]-0.05, length.out = gridDimCalvert)
  par6 <- seq(lower[6]+0.5, upper[6]-0.5, length.out = gridDimCalvert)
  
  grid <- expand.grid(par1, par2, par3, par4, par5, par6)
  rm(par1, par2, par3, par4, par5, par6)
  gc()
  
  print("Grid constructed. Dividing into segments")
  division <- length(grid[,1])/gridDivCalvert
  
  # CONSTRUCT THE GRIDS
  for(i in 1:gridDivCalvert){
    if (i == 1){
      gridTemp <- grid[1:division,]
    } else{
      gridTemp <- grid[((division*(i-1))+1):(division*(i)),]
    }
    saveRDS(gridTemp, file = paste0(outdir,"/grid",i,".Rda"))
    print(paste0("GRID ",i," GENERATED AND SAVED"))
    rm(gridTemp)
    gc()
  }
  
  rm(grid)
  gc()
  
  
  # EVALUATE EACH GRID
  for(i in 1:gridDivCalvert){
    print(paste0("Evaluating Grid ",i))
    gridTemp <- readRDS(file = paste0(outdir,"/grid",i,".Rda"))
    cl = makeCluster(c(rep("localhost", detectCores(logical = TRUE))),
                     type = "SOCK")
    clusterExport(cl, c("calvertCriterion", "calvertIndex", "loads"))
    gridTemp[,7] <- parRapply(cl, gridTemp, calvertLoss)
    stopCluster(cl)
    
    # Optimise over the grid
    numCores <- detectCores(logical = TRUE)
    cluster <- makeCluster(numCores)
    registerDoSNOW(cluster)
    calvertFitted <- foreach(
      i = 1:dim(gridTemp)[1],
      .combine = "rbind", .verbose = TRUE, .packages = "stats",
      .export = c("calvertCriterion", "calvertIndex", "loads",
                  "calvertPars", "calvertCompute", "MAPEfunc",
                  "RMSEfunc", "RSQfunc", "calvertConstraints",
                  "calvertLoss", "calvertFit")) %dopar% {
                    vals <- as.numeric(gridTemp[i,])
                    calvertFit(vals)
                  }
    stopCluster(cluster)
    saveRDS(calvertFitted, file = paste0(outdir,"/grid",i,".Rda"))
    rm(gridTemp)
    gc()
    print(paste0("GRID: ",i,"/", gridDivCalvert," evaluated"))
  }
  
  gc()
  
  print("All grid segments evaluated, compiling results")
  
  for(i in 1:gridDivCalvert){
    if (i == 1){
      grid <- readRDS(file = paste0(outdir,"/grid",i,".Rda"))
    } else{
      gridTemp <- readRDS(file = paste0(outdir,"/grid",i,".Rda"))
      grid <- rbind(grid, gridTemp)
      rm(gridTemp)
      gc()
    }
  }
  
  colnames(grid) <- 
    c("p0Start", "kgStart", "Tg1Start", "Tg2Start", "khStart",
      "ThStart", "startRSS", "p0Fitted", "kgFitted", "Tg1Fitted", "Tg2Fitted", "khFitted",
      "ThFitted", "optimRSS", "fnCount", "gnCount", "convCode",
      "EigenVal1", "EigenVal2", "EigenVal3", "EigenVal4", "EigenVal5", 
      "EigenVal6", "testRsq", "testRMSE", "testMAPE", "p0TrueDiff", 
      "kgTrueDiff", "Tg1TrueDiff", "Tg2TrueDiff", "khTrueDiff", "ThTrueDiff")
  
  saveRDS(grid, file = paste0(outdir,"/gridFull.Rda"))
  
  return(grid)
}