#' Project cells onto NMF factors
#'
#' @param mat Feature-by-cell nonnegative matrix.
#' @param nmf Fitted RcppML NMF model.
#'
#' @return Cell-by-dimension matrix of NMF factors.
#'
projectNMF = function(mat, nmf, scale_w = T) {
  if(scale_w){
    W = sweep(nmf@w, 2, nmf@d, `*`)
  }else {
    W = nmf@w
  }
  H = RcppML::nnls(
    w = W,
    A = mat
  )
  
  matNMF = t(H)
  
  rownames(matNMF) = colnames(mat)
  colnames(matNMF) = paste0(
    "NMF", seq_len(ncol(matNMF))
  )
  
  return(matNMF)
}