ggFootprint_modified <- function (seFoot = NULL, name = NULL, pal = NULL, smoothWindow = NULL,
          flank = NULL, flankNorm = NULL, baseSize = 6, normMethod = NULL,
          logFile = NULL)
{
  errorList <- list()
  rowDF <- SummarizedExperiment::rowData(seFoot)
  footMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[, 2] ==
                                                    "footprint"), ], name)
  biasMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[, 2] ==
                                                    "bias"), ], name)
  footDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "footprint"),
                  ]
  biasDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "bias"),
                  ]
  errorList$footMat <- footMat
  errorList$biasMat <- biasMat
  errorList$footDF <- footDF
  errorList$biasDF <- biasDF
  if (!is.null(smoothWindow)) {
    .logMessage("Applying smoothing window to footprint",
                logFile = logFile)
    footMat <- apply(footMat, 2, function(x) .centerRollMean(x,
                                                             smoothWindow))
    biasMat <- apply(biasMat, 2, function(x) .centerRollMean(x,
                                                             smoothWindow))
  }
  .logMessage("Normalizing by flanking regions", logFile = logFile)
  idx <- which(abs(footDF$x) >= flank - flankNorm)
  footMat <- t(t(footMat)/colMeans(footMat[idx, , drop = FALSE]))
  biasMat <- t(t(biasMat)/colMeans(biasMat[idx, , drop = FALSE]))
  errorList$footMatNorm <- footMat
  errorList$biasMatNorm <- footMat
  if (tolower(normMethod) == "none") {
    title <- ""
  }
  else if (tolower(normMethod) == "subtract") {
    title <- "Tn5 Bias Subtracted\n"
    footMat <- footMat - biasMat
  }
  else if (tolower(normMethod) == "divide") {
    title <- "Tn5 Bias Divided\n"
    footMat <- footMat/biasMat
  }
  else {
    stop("normMethod not recognized!")
  }
  .logMessage(paste0("NormMethod = ", normMethod), logFile = logFile)
  footMatMean <- .groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
  footMatSd <- .groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatMean <- .groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatSd <- .groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) .centerRollMean(x,
                                                                          11)))
  errorList$footMatMean <- footMatMean
  errorList$footMatSd <- footMatSd
  errorList$biasMatMean <- biasMatMean
  errorList$biasMatSd <- biasMatSd
  errorList$smoothFoot <- smoothFoot
  plotIdx <- seq_len(nrow(footMatMean))
  plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x) {
    data.frame(x = footDF$x, mean = footMatMean[, x], sd = footMatSd[,
                                                                     x], group = colnames(footMatMean)[x])[plotIdx, ,
                                                                                                           drop = FALSE]
  }) %>% Reduce("rbind", .)
  plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))
  plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x) {
    data.frame(x = biasDF$x, mean = biasMatMean[, x], sd = biasMatSd[,
                                                                     x], group = colnames(biasMatMean)[x])[plotIdx, ,
                                                                                                           drop = FALSE]
  }) %>% Reduce("rbind", .)
  plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))
  errorList$plotFootDF <- plotFootDF
  errorList$plotBiasDF <- plotBiasDF
  out <- tryCatch({
    if (is.null(pal)) {
      pal <- paletteDiscrete(values = gtools::mixedsort(SummarizedExperiment::colData(seFoot)$Group))
    }
    plotMax <- plotFootDF[order(plotFootDF$mean, decreasing = TRUE),
                          ]
    plotMax <- plotMax[abs(plotMax$x) > 20 & abs(plotMax$x) <
                         50, ]
    plotMax <- plotMax[!duplicated(plotMax$group), ]
    plotMax <- plotMax[seq_len(ceiling(nrow(plotMax)/4)),
                       ]
    plotMax$x <- 25
    ggFoot <- ggplot(plotFootDF, aes(x = x, y = mean, color = group)) +
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd,
                      linetype = NA, fill = group), alpha = 0.4) +
      geom_line() + scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) + xlab("Distance to motif center (bp)") +
      coord_cartesian(expand = FALSE, ylim = c(quantile(plotFootDF$mean,
                                                        1e-04), 1.15 * quantile(smoothFoot, 0.999)),
                      xlim = c(min(plotFootDF$x), max(plotFootDF$x))) +
      theme_ArchR(baseSize = baseSize) + ggtitle(name) +
      guides(fill = FALSE) + guides(color = FALSE) + ylab(paste0(title,
                                                                 "Normalized Insertions")) + ggrepel::geom_label_repel(data = plotMax,
                                                                                                                       aes(label = group), size = 3, xlim = c(75, NA))
  }, error = function(e) {
    cat('fail\n')
  })
  ggFoot
}

