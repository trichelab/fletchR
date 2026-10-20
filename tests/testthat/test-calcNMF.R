# makeNMFTestMatrix = function() {
#   set.seed(42)
#   
#   W = matrix(
#     runif(40 * 3, min = 0.1, max = 1),
#     nrow = 40,
#     ncol = 3
#   )
#   
#   H = matrix(
#     runif(3 * 20, min = 0.1, max = 1),
#     nrow = 3,
#     ncol = 20
#   )
#   
#   mat = W %*% H
#   
#   rownames(mat) = paste0("feature", seq_len(nrow(mat)))
#   colnames(mat) = paste0("cell", seq_len(ncol(mat)))
#   
#   return(mat)
# }

mat = readRDS(testthat::test_path("..", "testdata", "logTFIDF_iter2.rds"))#makeNMFTestMatrix()

describe("RcppML::nmf seed interactions",{
  it("seed yields same result",{
    nmfFit = RcppML::nmf(
      data = mat,
      k = 3,
      seed = 1,
      maxit = 500
    )
    
    nmfFit2 = RcppML::nmf(
      data = mat,
      k = 3,
      seed = 1,
      maxit = 500
    )
    testthat::expect_all_true(c(
        identical(nmfFit@w,nmfFit2@w),
        identical(nmfFit@d,nmfFit2@d),
        identical(nmfFit@h,nmfFit2@h))
    )
  })
  it("different seed yield different result",{
    nmfFit = RcppML::nmf(
      data = mat,
      k = 3,
      seed = 1,
      maxit = 500
    )
    
    nmfFit2 = RcppML::nmf(
      data = mat,
      k = 3,
      seed = 2,
      maxit = 500
    )
    testthat::expect_all_false(c(
      identical(nmfFit@w,nmfFit2@w),
      identical(nmfFit@d,nmfFit2@d),
      identical(nmfFit@h,nmfFit2@h))
    )
    })
  it("nmf seed independent",{
    set.seed(1)
    nmfFit = RcppML::nmf(
      data = mat,
      k = 3,
      seed = 1,
      maxit = 500
    )
    
    set.seed(2)
    nmfFit2 = RcppML::nmf(
      data = mat,
      k = 3,
      seed = 1,
      maxit = 500
    )
    testthat::expect_all_true(c(
      identical(nmfFit@w,nmfFit2@w),
      identical(nmfFit@d,nmfFit2@d),
      identical(nmfFit@h,nmfFit2@h))
    )
    })
  it("lem yield and independent to outer seed",{
    set.seed(1)
    lem = calcNMF(
      mat = mat,
      nDimensions = 3,
      seed = 1
    )
    
    set.seed(2)
    lem2 = calcNMF(
      mat = mat,
      nDimensions = 3,
      seed = 1
    )
    
    testthat::expect_all_true(c(
      identical(lem@metadata$nmf@w,lem2@metadata$nmf@w),
      identical(lem@metadata$nmf@d,lem2@metadata$nmf@d),
      identical(lem@metadata$nmf@h,lem2@metadata$nmf@h))
    )
    })
  it("different seed yield different nmf",{
    lem = calcNMF(
      mat = mat,
      nDimensions = 3,
      seed = 1
    )
    
    lem2 = calcNMF(
      mat = mat,
      nDimensions = 3,
      seed = 2
    )
    
    testthat::expect_all_false(c(
      identical(lem@metadata$nmf@w,lem2@metadata$nmf@w),
      identical(lem@metadata$nmf@d,lem2@metadata$nmf@d),
      identical(lem@metadata$nmf@h,lem2@metadata$nmf@h))
    )
  })
  it("nmf rresult is same in and out of calcNMF",{
    nmfFit = RcppML::nmf(
      data = mat,
      k = 3,
      seed = 1,
    )
    
    lem = calcNMF(
      mat = mat,
      nDimensions = 3,
      seed = 1
    )
    
    testthat::expect_all_true(c(
      identical(lem@metadata$nmf@w,nmfFit@w),
      identical(lem@metadata$nmf@d,nmfFit@d),
      identical(lem@metadata$nmf@h,nmfFit@h))
    )
    })
})


describe("Projection scaled vs unscaled w",{
  lem = calcNMF(
    mat = mat,
    nDimensions = 5,
    seed = 3
  )
  lem2 = calcNMF(
    mat = mat,
    nDimensions = 5,
    seed = 1
  )
  
  it("projection is similar to original H",{
    testthat::expect_all_true(
      (cor(SingleCellExperiment::sampleFactors(lem),t(lem@metadata$nmf@h)) |> Matrix::diag()) > 0.99
    )
    testthat::expect_equal(
      unname(SingleCellExperiment::sampleFactors(lem)),
      unname(t(lem@metadata$nmf@h)),
      tolerance = 1e-2
    )
  })
  it("unscaled w yields scaled h",{
    projected_unscaled = projectNMF(
      mat = mat,
      nmf = lem@metadata$nmf,
      scale_w = F
    )
    
    testthat::expect_all_true(
      (cor(projected_unscaled,t(lem@metadata$nmf@h)) |> Matrix::diag()) > 0.99
    )
    
    testthat::expect_equal(
      unname(sweep(projected_unscaled,2,lem@metadata$nmf@d,`/`)),
      unname(t(lem@metadata$nmf@h)),
      tolerance = 1e-2
    )
  
  })
  
  
 
  # it("factorizes similarly across seed",{
  #   testthat::expect_all_true(
  #     (cor(
  #       SingleCellExperiment::sampleFactors(lem),
  #       SingleCellExperiment::sampleFactors(lem2)
  #       ) |> Matrix::diag()
  #      ) > 0.99
  #   )
  #   
  #   testthat::expect_equal(
  #     unname(SingleCellExperiment::sampleFactors(lem)),
  #     unname(SingleCellExperiment::sampleFactors(lem2)),
  #     tolerance = 1e-2
  #   )
  #   })
})


# 
# describe("projectNMF", {
#   
#   it("projects all cells into the expected dimensions", {
#     
#     result = projectNMF(
#       mat = mat,
#       nmf = nmfFit
#     )
#     
#     expect_equal(dim(result), c(20L, 3L))
#     expect_identical(rownames(result), colnames(mat))
#     expect_identical(colnames(result), paste0("NMF", 1:3))
#   })
#   
#   it("returns finite and nonnegative cell factors", {
#     
#     result = projectNMF(
#       mat = mat,
#       nmf = nmfFit
#     )
#     
#     expect_true(all(is.finite(result)))
#     expect_true(all(result >= -1e-10))
#   })
#   
#   it("reconstructs a low-rank nonnegative matrix", {
#     
#     result = projectNMF(
#       mat = mat,
#       nmf = nmfFit
#     )
#     
#     W = sweep(nmfFit@w, 2, nmfFit@d, `*`)
#     reconstructed = W %*% t(result)
#     
#     relativeError = sqrt(
#       sum((mat - reconstructed)^2) / sum(mat^2)
#     )
#     
#     expect_lt(relativeError, 0.1)
#   })
#   
#   it("preserves cell order after projection", {
#     
#     reorderedMat = mat[, rev(seq_len(ncol(mat))), drop = FALSE]
#     
#     result = projectNMF(
#       mat = reorderedMat,
#       nmf = nmfFit
#     )
#     
#     expect_identical(
#       rownames(result),
#       colnames(reorderedMat)
#     )
#   })
# })
# 
# 
# describe("calcNMF", {
#   
#   it("returns a LinearEmbeddingMatrix", {
#     
#     expect_s4_class(lem, "LinearEmbeddingMatrix")
#     expect_true("nmf" %in% names(S4Vectors::metadata(lem)))
#   })
#   
#   it("preserves cell names and dimension names", {
#     
#     expect_equal(dim(projected), c(20L, 3L))
#     expect_identical(rownames(projected), colnames(mat))
#     expect_identical(colnames(projected), paste0("NMF", 1:3))
#     expect_identical(dimnames(lem), dimnames(projected))
#   })
#   
#   it("returns finite and nonnegative embeddings", {
#     
#     expect_true(all(is.finite(projected)))
#     expect_true(all(projected >= -1e-10))
#   })
#   
#   it("projects outlier cells without excluding them from the result", {
#     
#     outliers = rep(FALSE, ncol(mat))
#     outliers[c(2, 5, 10)] = TRUE
#     
#     outlierLem = calcNMF(
#       mat = mat,
#       nDimensions = 3,
#       outliers = outliers
#     )
#     
#     outlierProjected = SingleCellExperiment::sampleFactors(outlierLem)
#     
#     expect_equal(nrow(outlierProjected), ncol(mat))
#     expect_identical(rownames(outlierProjected), colnames(mat))
#     
#     expect_true(all(
#       is.finite(outlierProjected[outliers, , drop = FALSE])
#     ))
#     
#     expect_true(all(
#       outlierProjected[outliers, , drop = FALSE] >= -1e-10
#     ))
#   })
#   
#   it("does not use outlier values to learn feature factors", {
#     
#     outliers = rep(FALSE, ncol(mat))
#     outliers[c(2, 5, 10)] = TRUE
#     
#     modifiedMat = mat
#     modifiedMat[, outliers] = modifiedMat[, outliers] * 100
#     
#     set.seed(123)
#     lem1 = calcNMF(
#       mat = mat,
#       nDimensions = 3,
#       outliers = outliers
#     )
#     
#     set.seed(123)
#     lem2 = calcNMF(
#       mat = modifiedMat,
#       nDimensions = 3,
#       outliers = outliers
#     )
#     
#     nmf1 = S4Vectors::metadata(lem1)$nmf
#     nmf2 = S4Vectors::metadata(lem2)$nmf
#     
#     expect_equal(nmf1@w, nmf2@w, tolerance = 1e-6)
#     expect_equal(nmf1@d, nmf2@d, tolerance = 1e-6)
#     
#     factors1 = SingleCellExperiment::sampleFactors(lem1)
#     factors2 = SingleCellExperiment::sampleFactors(lem2)
#     
#     expect_equal(
#       factors1[!outliers, , drop = FALSE],
#       factors2[!outliers, , drop = FALSE],
#       tolerance = 1e-6
#     )
#   })
#   
#   it("removes zero-sum features before fitting", {
#     
#     modifiedMat = mat
#     modifiedMat[1, ] = 0
#     
#     modifiedLem = calcNMF(
#       mat = modifiedMat,
#       nDimensions = 3
#     )
#     
#     nmf = S4Vectors::metadata(modifiedLem)$nmf
#     
#     expect_equal(nrow(nmf@w), nrow(mat) - 1L)
#     
#     expect_equal(
#       nrow(SingleCellExperiment::sampleFactors(modifiedLem)),
#       ncol(mat)
#     )
#   })
#   
#   it("removes features that are zero in non-outlier cells", {
#     
#     outliers = rep(FALSE, ncol(mat))
#     outliers[20] = TRUE
#     
#     modifiedMat = mat
#     modifiedMat[1, !outliers] = 0
#     modifiedMat[1, outliers] = 100
#     
#     modifiedLem = calcNMF(
#       mat = modifiedMat,
#       nDimensions = 3,
#       outliers = outliers
#     )
#     
#     nmf = S4Vectors::metadata(modifiedLem)$nmf
#     
#     expect_equal(nrow(nmf@w), nrow(mat) - 1L)
#   })
# })
