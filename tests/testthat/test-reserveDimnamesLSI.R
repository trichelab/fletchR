describe("projectSVD", {
  
  it("preserves row and column names", {
    # set.seed(1)
    
    mat = matrix(
      rpois(100 * 20, lambda = 2),
      nrow = 100,
      ncol = 20,
      dimnames = list(
        paste0("peak", seq_len(100)),
        paste0("cell", seq_len(20))
      )
    )
    
    svd = irlba::irlba(mat, nv = 5, nu = 5)
    projected = projectSVD(mat, svd$u, svd$d, 5)
    
    expect_identical(rownames(projected), colnames(mat))
    expect_identical(colnames(projected), paste0("LSI", seq_len(5)))
  })
  
})

describe("calcLSI", {
  # set.seed(1)
  mat = matrix(
    rpois(100 * 20, lambda = 2),
    nrow = 100,
    ncol = 20,
    dimnames = list(
      paste0("peak", seq_len(100)),
      paste0("cell", seq_len(20))
    )
  )
  ndim = 5
  result = calcLSI(mat, nDimensions = ndim)
  
  it("preserves row and column names in LinearEmbeddingMatrix and sampleFactors", {

    # Check LinearEmbeddingMatrix
    expect_s4_class(result, "LinearEmbeddingMatrix")
    expect_identical(rownames(result), colnames(mat))
    expect_identical(colnames(result), paste0("LSI", seq_len(ndim)))
    
    # Check sampleFactors
    embeddings = SingleCellExperiment::sampleFactors(result)
    
    expect_identical(rownames(embeddings), colnames(mat))
    expect_identical(colnames(embeddings), paste0("LSI", seq_len(ndim)))
  })
  
  it("add dimnames in svd result",{
    expect_identical(rownames(result@metadata$svd$u), rownames(mat))
    expect_identical(rownames(result@metadata$svd$v), colnames(mat))
    })
})
