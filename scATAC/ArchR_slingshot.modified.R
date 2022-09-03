#Run This to accesss hidden functions
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
    tryCatch({
        eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
}




return_SlingShotTrajectories <- function(
  ArchRProj = NULL,
  name = "SlingShot",
  useGroups = NULL,
  principalGroup = NULL,
  groupBy = NULL,
  embedding = NULL,
  reducedDims = NULL,
  force = FALSE,
  seed = 1,
  stretch = 0,approx_points = 100
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = useGroups, name = "useGroups", valid = c("character"))
  .validInput(input = principalGroup, name = "principalGroup", valid = c("character"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = embedding, name = "embedding", valid = c("character", "null"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character", "null"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("numeric"))

  .requirePackage("slingshot")

  set.seed(seed)

  if(!is.null(embedding)){
    rD <- getEmbedding(ArchRProj, embedding = embedding)
  }else{
    rD <- getReducedDims(ArchRProj, reducedDims = reducedDims)
  }

  groups <- getCellColData(ArchRProj, groupBy)

  if(!is.null(useGroups)){
    idx <- which(groups[,1] %in% useGroups)
    rD <- rD[idx, , drop = FALSE]
    groups <- groups[idx, , drop = FALSE]
  }

  sds <- slingshot(
    data = rD, 
    clusterLabels = groups[rownames(rD), ], 
      start.clus = principalGroup,
	stretch = stretch,
		approx_points = approx_points
  )

  sds

}
