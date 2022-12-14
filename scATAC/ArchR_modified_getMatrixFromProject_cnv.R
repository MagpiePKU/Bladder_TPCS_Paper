getMatrixFromProject_mod_cnv <- function (ArchRProj = NULL, useMatrix = "CNVMatrix", useSeqnames = NULL, 
                                      verbose = TRUE, binarize = FALSE, threads = getArchRThreads(), 
                                      logFile = createLogFile("getMatrixFromProject")) 
{
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character", 
                                                                   "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()), sys.frame(sys.nframe())), 
           "getMatrixFromProject Input-Parameters", logFile = logFile)
  ArrowFiles <- getArrowFiles(ArchRProj)
  cellNames <- ArchRProj$cellNames
  seL <- .safelapply(seq_along(ArrowFiles), function(x) {
    .logDiffTime(paste0("Reading ", useMatrix, " : ", names(ArrowFiles)[x], 
                        "(", x, " of ", length(ArrowFiles), ")"), t1 = tstart, 
                 verbose = FALSE, logFile = logFile)
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], 
                                subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    if (length(allCells) != 0) {
      o <- getMatrixFromArrow(ArrowFile = ArrowFiles[x], 
                              useMatrix = useMatrix, useSeqnames = useSeqnames, 
                              cellNames = allCells, ArchRProj = ArchRProj, 
                              verbose = FALSE, binarize = binarize, logFile = logFile)
      .logDiffTime(paste0("Completed ", useMatrix, " : ", 
                          names(ArrowFiles)[x], "(", x, " of ", length(ArrowFiles), 
                          ")"), t1 = tstart, verbose = FALSE, logFile = logFile)
      o
    }
    else {
      NULL
    }
  }, threads = threads)
  .logDiffTime("Organizing colData", t1 = tstart, verbose = verbose, 
               logFile = logFile)
  cD <- lapply(seq_along(seL), function(x) {
    colData(seL[[x]]) -> y
    colnames(y)[length(colnames(y))] <- 'seurat_clusters'
    y
  }) %>% Reduce("rbind", .)
  .logDiffTime("Organizing rowData", t1 = tstart, verbose = verbose, 
               logFile = logFile)
  rD1 <- (seL[[1]]@rowRanges)
  rD <- lapply(seq_along(seL), function(x) {
    identical((seL[[x]]@rowRanges), rD1)
  }) %>% unlist %>% all
  if (!rD) {
    stop("Error with rowData being equal for every sample!")
  }
  nAssays <- names(assays(seL[[1]]))
  asy <- lapply(seq_along(nAssays), function(i) {
    .logDiffTime(sprintf("Organizing Assays (%s of %s)", 
                         i, length(nAssays)), t1 = tstart, verbose = verbose, 
                 logFile = logFile)
    m <- lapply(seq_along(seL), function(j) {
      assays(seL[[j]])[[nAssays[i]]]
    }) %>% Reduce("cbind", .)
    m
  }) %>% SimpleList()
  names(asy) <- nAssays
  .logDiffTime("Constructing SummarizedExperiment", t1 = tstart, 
               verbose = verbose, logFile = logFile)
  se <- SummarizedExperiment(assays = asy, colData = cD, rowData = rD1)
  rm(seL)
  gc()
  .logDiffTime("Finished Matrix Creation", t1 = tstart, verbose = verbose, 
               logFile = logFile)
  se
}
