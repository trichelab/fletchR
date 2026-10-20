#' add a uwot::umap projection to reducedDims(x) AND SAVE THE MODEL TO METADATA
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param reducedDims   name of reducedDim used (default: NMF if found else LSI)
#' @param name          name of reducedDim to store the UMAP embedding ("UMAP") 
#' @param nNeighbors    integer nearest neighbors to use for UMAP (40)
#' @param minDist       passed to uwot::umap(min_dist); controls packing (0.4)
#' @param metric        distance metric for uwot::umap to use ("cosine") 
#' @param model         precomputed UMAP model to use for projection (NULL)
#' @param ...           more arguments to pass to uwot::umap()/_transform()
#'
#' @return              SingleCellExperiment with reducedDim(x, "NMF")
#'
#' @details             uses the same default neighbors (40) as ArchR by default
#'
#' @seealso             addNMF
#' @seealso             iSEEarchR
#' @seealso             uwot::umap 
#'
#' @import              SingleCellExperiment
#' @import              uwot
#'
#' @export
#'
addUMAP <- function(x, reducedDims=NULL, name="UMAP", nNeighbors=40, minDist=0.4, metric="cosine", model=NULL, ...) { 
  
  if (is(x, 'ArchRProject')) return(ArchR::addUMAP(x, what, ...))

  if (is.null(reducedDims)) {
    if ("NMF" %in% reducedDimNames(x)) {
      reducedDims <- "NMF"
    } else if ("LSI" %in% reducedDimNames(x)) {
      reducedDims <- "NMF"
    } else {
      stop("Neither NMF nor LSI in reducedDimNames(x), need to supply one")
    }
  }
  rd = reducedDim(x, reducedDims)
  if(inherits(rd,"LinearEmbeddingMatrix")){
    rd = SingleCellExperiment::sampleFactors(rd)
  }
  if (!is.null(model)) {
    message("Using pre-computed model for UMAP projection...")
    res <- model
    message("Computing umap_transform() with ", reducedDims, " as input...")
    
    res$embedding <- uwot::umap_transform(rd, 
                                          model=model, ...)
  } else {
    message("Computing UMAP model using ", reducedDims, " as input...")
    res <- uwot::umap(rd, n_neighbors=nNeighbors, 
                      metric=metric, min_dist=minDist, ret_model=TRUE, ...)
  }
  message("Saving model to metadata(x)$UMAP...")
  metadata(x)$UMAP <- res 
  message("Adding embeddings to reducedDim(x, '", name, "')...")
  reducedDim(x, name) <- res$embedding
  message("Done.")
  return(x)

}
