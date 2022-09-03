
modified_availableSeqnames <- function (ArrowFiles = NULL, subGroup = "Fragments", threads = getArchRThreads())
{
  threads <- min(threads, length(ArrowFiles))
  o <- h5closeAll()
  seqList <- .safelapply(seq_along(ArrowFiles), function(x) {
    seqnames <- h5ls(ArrowFiles[x]) %>% {
      .[.$group == paste0("/", subGroup), ]$name
    }
    seqnames <- seqnames[!grepl("Info", seqnames)]
    seqnames
  }, threads = threads)
  # if (!all(unlist(lapply(seq_along(seqList), function(x) identical(seqList[[x]],
  #                                                                  seqList[[1]]))))) {
  #   stop("Not All Seqnames Identical!")
  # }
  o <- h5closeAll()
  return(paste0(seqList[[1]]))
}


modified_add_iterative_LSI <- function (ArchRProj = NULL, useMatrix = "TileMatrix", name = "IterativeLSI", 
                                    iterations = 2, clusterParams = list(resolution = c(0.2), 
                                                                         sampleCells = 10000, n.start = 10), varFeatures = 25000, 
                                    dimsToUse = 1:30, LSIMethod = 2, scaleDims = TRUE, corCutOff = 0.75, 
                                    binarize = TRUE, outlierQuantiles = c(0.02, 0.98), filterBias = TRUE, 
                                    sampleCellsPre = 10000, projectCellsPre = FALSE, sampleCellsFinal = NULL, 
                                    selectionMethod = "var", scaleTo = 10000, totalFeatures = 5e+05, 
                                    filterQuantile = 0.9, excludeChr = c(), saveIterations = TRUE, 
                                    UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", 
                                                      verbose = FALSE, fast_sgd = TRUE), nPlot = 10000, outDir = getOutputDirectory(ArchRProj), 
                                    threads = getArchRThreads(), seed = 1, verbose = TRUE, force = FALSE, 
                                    logFile = createLogFile("modified_addIterativeLSI")) 
{
  if (verbose) 
    message("Checking Inputs...")
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = iterations, name = "iterations", valid = c("integer"))
  .validInput(input = clusterParams, name = "clusterParams", 
              valid = c("list"))
  .validInput(input = varFeatures, name = "varFeatures", valid = c("integer"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer"))
  .validInput(input = LSIMethod, name = "LSIMethod", valid = c("integer", 
                                                               "character"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", 
                                                               "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = outlierQuantiles, name = "outlierQuantiles", 
              valid = c("numeric", "null"))
  .validInput(input = filterBias, name = "filterBias", valid = c("boolean"))
  .validInput(input = sampleCellsPre, name = "sampleCellsPre", 
              valid = c("integer", "null"))
  .validInput(input = sampleCellsFinal, name = "sampleCellsFinal", 
              valid = c("integer", "null"))
  .validInput(input = selectionMethod, name = "selectionMethod", 
              valid = c("character"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = totalFeatures, name = "totalFeatures", 
              valid = c("integer"))
  .validInput(input = filterQuantile, name = "filterQuantile", 
              valid = c("numeric"))
  .validInput(input = excludeChr, name = "excludeChr", valid = c("character", 
                                                                 "null"))
  .validInput(input = saveIterations, name = "saveIterations", 
              valid = c("boolean"))
  .validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))
  .validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  .validInput(input = outDir, name = "outDir", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  if (varFeatures < 1000) {
    stop("Please provide more than 1000 varFeatures!")
  }
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()), sys.frame(sys.nframe())), 
           "IterativeLSI Input-Parameters", logFile = logFile)
  .requirePackage("Matrix", source = "cran")
  tstart <- Sys.time()
  if (!is.null(ArchRProj@reducedDims[[name]])) {
    if (!force) {
      stop("Error name in reducedDims Already Exists! Set force = TRUE or pick a different name!")
    }
  }
  set.seed(seed)
  outDir <- file.path(outDir, name)
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  cellNames <- rownames(getCellColData(ArchRProj))
  if (!is.null(sampleCellsPre)) {
    if (length(cellNames) < sampleCellsPre) {
      sampleCellsPre <- NULL
    }
  }
  if (!is.null(sampleCellsFinal)) {
    if (length(cellNames) < sampleCellsFinal) {
      sampleCellsFinal <- NULL
    }
  }
  ArrowFiles <- getSampleColData(ArchRProj)[, "ArrowFiles"]
  # stopifnot(any(tolower(useMatrix) %in% c("tilematrix", "peakmatrix")))
  if (tolower(useMatrix) == "tilematrix") {
    useMatrix <- "TileMatrix"
    tileSizes <- lapply(ArrowFiles, function(x) {
      h5read(x, "TileMatrix/Info/Params/")$tileSize[1]
    }) %>% unlist
    if (length(unique(tileSizes)) != 1) {
      stop("Error not all TileMatrices are the same tileSize!")
    }
    tileSize <- unique(tileSizes)
  }
  if (tolower(useMatrix) == "peakmatrix") {
    useMatrix <- "PeakMatrix"
    tileSize <- NA
  }
  if (tolower(useMatrix) %ni% c("tilematrix","peakmatrix")) {
    useMatrix <- useMatrix
    tileSize <- NA
  }
  tstart <- Sys.time()
  chrToRun <- modified_availableSeqnames(ArrowFiles, subGroup = useMatrix)
  .logDiffTime("Computing Total Accessibility Across All Features", 
               tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  if (useMatrix == "TileMatrix") {
    totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, 
                            seqnames = chrToRun, addInfo = FALSE)
    totalAcc$start <- (totalAcc$idx - 1) * tileSize
  }
  else {
    totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, 
                            seqnames = chrToRun, addInfo = TRUE)
  }
  gc()
  cellDepth <- tryCatch({
    df <- getCellColData(ArchRProj = ArchRProj, select = "nFrags")
    v <- df[, 1]
    names(v) <- rownames(df)
    v
  }, error = function(e) {
    tryCatch({
      .getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, 
                  seqnames = chrToRun)
    }, error = function(y) {
      stop("Could not determine depth from nFrags or colSums!")
    })
  })
  cellDepth <- log10(cellDepth + 1)
  if (length(excludeChr) > 0) {
    totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% 
                                               excludeChr), , drop = FALSE]
  }
  .logDiffTime("Computing Top Features", tstart, addHeader = FALSE, 
               verbose = verbose, logFile = logFile)
  nFeature <- varFeatures[1]
  rmTop <- floor((1 - filterQuantile) * totalFeatures)
  topIdx <- head(order(totalAcc$rowSums, decreasing = TRUE), 
                 nFeature + rmTop)[-seq_len(rmTop)]
  topFeatures <- totalAcc[sort(topIdx), ]
  .logDiffTime(paste0("Running LSI (1 of ", iterations, ") on Top Features"), 
               tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)
  j <- 1
  if (!is.null(clusterParams$sampleCells)) {
    if (!is.na(clusterParams$sampleCells[j])) {
      sampleJ <- clusterParams$sampleCells[j]
    }
    else if (!is.na(clusterParams$sampleCells[1])) {
      sampleJ <- clusterParams$sampleCells[1]
    }
    else {
      sampleJ <- sampleCellsPre
    }
  }
  else {
    sampleJ <- sampleCellsPre
  }
  outLSI <- .LSIPartialMatrix(ArrowFiles = ArrowFiles, featureDF = topFeatures, 
                              cellNames = cellNames, cellDepth = cellDepth, useMatrix = useMatrix, 
                              sampleNames = getCellColData(ArchRProj)$Sample, LSIMethod = LSIMethod, 
                              scaleTo = scaleTo, dimsToUse = dimsToUse, binarize = binarize, 
                              outlierQuantiles = outlierQuantiles, sampleCells = if (j != 
                                                                                     iterations) 
                                sampleCellsPre
                              else sampleCellsFinal, projectAll = j == iterations | 
                                projectCellsPre | sampleJ > sampleCellsPre, threads = threads, 
                              useIndex = FALSE, tstart = tstart, verbose = verbose, 
                              logFile = logFile)
  outLSI$scaleDims <- scaleDims
  outLSI$useMatrix <- useMatrix
  outLSI$tileSize <- tileSize
  gc()
  .logThis(outLSI, paste0("outLSI-", j), logFile = logFile)
  if (iterations == 1) {
    .logDiffTime("Finished Running IterativeLSI", tstart, 
                 addHeader = FALSE, verbose = verbose, logFile = logFile)
    ArchRProj@reducedDims[[name]] <- outLSI
    return(ArchRProj)
  }
  clusterDF <- .LSICluster(outLSI = outLSI, filterBias = filterBias, 
                           cellNames = cellNames, cellDepth = cellDepth, dimsToUse = dimsToUse, 
                           scaleDims = scaleDims, corCutOff = corCutOff, clusterParams = clusterParams, 
                           j = j, verbose = verbose, tstart = tstart, logFile = logFile)
  clusters <- clusterDF$clusters
  nClust <- length(unique(clusters))
  .logDiffTime(sprintf("Identified %s Clusters", nClust), tstart, 
               addHeader = FALSE, verbose = verbose, logFile = logFile)
  .logThis(clusterDF, paste0("clusterDF-", j), logFile = logFile)
  if (saveIterations) {
    .logDiffTime("Saving LSI Iteration", tstart, addHeader = FALSE, 
                 verbose = verbose, logFile = logFile)
    .saveIteration(outLSI = outLSI, clusters = clusters, 
                   scaleDims = scaleDims, dimsToUse = dimsToUse, corCutOff = corCutOff, 
                   outDir = outDir, nPlot = nPlot, UMAPParams = UMAPParams, 
                   ArchRProj = ArchRProj, j = j, threads = threads, 
                   logFile = logFile)
  }
  variableFeatures <- topFeatures
  while (j < iterations) {
    j <- j + 1
    variableFeatures <- .identifyVarFeatures(outLSI = outLSI, 
                                             clusterDF = clusterDF, ArrowFiles = ArrowFiles, useMatrix = useMatrix, 
                                             prevFeatures = variableFeatures, scaleTo = scaleTo, 
                                             totalAcc = totalAcc, totalFeatures = totalFeatures, 
                                             selectionMethod = selectionMethod, varFeatures = varFeatures, 
                                             tstart = tstart, threads = threads, verbose = verbose, 
                                             logFile = logFile)
    .logDiffTime(sprintf("Running LSI (%s of %s) on Variable Features", 
                         j, iterations), tstart, addHeader = TRUE, verbose = verbose, 
                 logFile = logFile)
    if (!is.null(clusterParams$sampleCells)) {
      if (!is.na(clusterParams$sampleCells[j])) {
        sampleJ <- clusterParams$sampleCells[j]
      }
      else if (!is.na(clusterParams$sampleCells[1])) {
        sampleJ <- clusterParams$sampleCells[1]
      }
      else {
        sampleJ <- sampleCellsPre
      }
    }
    else {
      sampleJ <- sampleCellsPre
    }
    outLSI <- .LSIPartialMatrix(ArrowFiles = ArrowFiles, 
                                featureDF = variableFeatures, useMatrix = useMatrix, 
                                cellNames = cellNames, cellDepth = cellDepth, sampleNames = getCellColData(ArchRProj)$Sample, 
                                LSIMethod = LSIMethod, scaleTo = scaleTo, dimsToUse = dimsToUse, 
                                binarize = binarize, outlierQuantiles = outlierQuantiles, 
                                sampleCells = if (j != iterations) 
                                  sampleCellsPre
                                else sampleCellsFinal, projectAll = j == iterations | 
                                  projectCellsPre | sampleJ > sampleCellsPre, threads = threads, 
                                useIndex = FALSE, tstart = tstart, verbose = verbose, 
                                logFile = logFile)
    outLSI$scaleDims <- scaleDims
    outLSI$useMatrix <- useMatrix
    outLSI$tileSize <- tileSize
    .logThis(outLSI, paste0("outLSI-", j), logFile = logFile)
    if (j != iterations) {
      clusterDF <- .LSICluster(outLSI = outLSI, dimsToUse = dimsToUse, 
                               scaleDims = scaleDims, corCutOff = corCutOff, 
                               filterBias = filterBias, cellNames = cellNames, 
                               cellDepth = cellDepth, j = j, clusterParams = clusterParams, 
                               verbose = verbose, tstart = tstart, logFile = logFile)
      clusters <- clusterDF$clusters
      nClust <- length(unique(clusters))
      .logDiffTime(sprintf("Identified %s Clusters", nClust), 
                   tstart, addHeader = FALSE, verbose = verbose, 
                   logFile = logFile)
      .logThis(clusterDF, paste0("clusterDF-", j), logFile = logFile)
      if (saveIterations) {
        .logDiffTime("Saving LSI Iteration", tstart, 
                     addHeader = FALSE, verbose = verbose, logFile = logFile)
        .saveIteration(outLSI = outLSI, clusters = clusters, 
                       scaleDims = scaleDims, dimsToUse = dimsToUse, 
                       corCutOff = corCutOff, outDir = outDir, nPlot = nPlot, 
                       UMAPParams = UMAPParams, ArchRProj = ArchRProj, 
                       j = j, threads = threads, logFile = logFile)
      }
    }
    gc()
  }
  .logDiffTime("Finished Running IterativeLSI", tstart, addHeader = FALSE, 
               verbose = verbose, logFile = logFile)
  ArchRProj@reducedDims[[name]] <- outLSI
  return(ArchRProj)
}
