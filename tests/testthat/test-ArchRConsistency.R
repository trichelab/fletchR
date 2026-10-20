describe("FletchR", {
  library(ArchR)
  proj <- ArchR::getTestProject()
  proj <- ArchR::addIterativeLSI(
    proj,
    dimsToUse = 1:5,
    varFeatures = 1000,
    force = TRUE
  )
  
  SCE <- fletchR::archRtoSCE(proj)
  # SCE = addTfIdf(SCE,useMatrix = "counts",subsetLSI = T,assayName = "TfIdf",binarize=T,outlierQuantiles = c(0,1))
  SCE = fletchR::addLSI(
    SCE,
    useMatrix = "ArchRTfIdf",
    subsetLSI = T,
    scaleDims = FALSE,
    nDimensions = 5,
    corCutOff = 0.75
  )
  ArchRLSI = SCE@metadata$LSI
  
  # c = cor(reducedDim(SCE,"LSI")[rownames(ArchRLSI$matSVD),],ArchRLSI$matSVD) |> diag()
  
  logtfidftest = fletchR:::filterAndGetMat(SCE, useMatrix = "ArchRTfIdf", replaceZeros = F)
  
  outliers = rownames(ArchRLSI$matSVD) %in% ArchRLSI$outliers
  tt = logtfidftest[, !outliers]
  tt = tt[rowSums(tt) > 0, ]
  
  z = projectSVD(
    tt,
    u = ArchRLSI$svd$u,
    d = ArchRLSI$svd$d,
    v = ArchRLSI$svd$v,
    nDim = 5
  )
  d = projectSVD(tt,
                 u = ArchRLSI$svd$u,
                 d = ArchRLSI$svd$d,
                 nDim = 5)
  
  #feature selection
  archRTopFeats = readRDS(testthat::test_path("..", "testdata", "first_pass_topFeatures.rds"))
  archRTopFeats$end = archRTopFeats$start + 10000 - 1 #tileSize = 10000 for this set
  archRTopFeats = archRTopFeats |> as("GRanges") |> as.character()
  
  first_pass_clust = readRDS(testthat::test_path("..", "testdata", "first_pass_clust.rds"))
  first_pass_LSI = readRDS(testthat::test_path("..", "testdata", "first_pass_LSI.rds"))
  
  second_pass_LSI = readRDS(testthat::test_path("..", "testdata", "second_pass_LSI.rds"))
  
  second_pass_variableFeatures = readRDS(testthat::test_path("..", "testdata", "second_pass_variableFeatures.rds"))
  second_pass_variableFeatures$end = second_pass_variableFeatures$start + 10000 - 1 #tileSize = 10000 for this set
  ArchRVarFeat = as.character(as(second_pass_variableFeatures, "GRanges"))
  
  set.seed(1)
  
  LSIres1 = .firstPass(
    x = SCE,
    useMatrix = "counts",
    preTfIdf = NULL,
    binarize = T,
    # excludeChr = excludeChr,
    nDimensions = 5,
    # filterQuantile = filterQuantile,
    outlierQuantile = c(0, 1),
    varFeatures = 1000,
    seed = 1
    # totalFeatures = totalFeatures,
    # depthCol = depthCol,
    # scaleTo = scaleTo
  )
  nFrags = colData(SCE)[["nFrags"]]
  
  clustF = clusterSeurat(mat = removeDepthCorrelatedDims(ArchRRowZscores(SingleCellExperiment::sampleFactors(LSIres1)), depth = nFrags))
  
  set.seed(1)
  LSIres2 = .iterativePass(
    x = SCE,
    useMatrix = "counts",
    preTfIdf = NULL,
    clusters = clustF, # passing true results to make sure the second test works
    binarize = T,
    # excludeChr = excludeChr,
    nDimensions = 5,
    # filterQuantile = filterQuantile,
    outlierQuantile = c(0, 1),
    varFeatures = 1000,
    seed = 1
    # totalFeatures = totalFeatures,
    # depthCol = depthCol,
    # scaleTo = scaleTo
  )
  
  LSIres_all = featureSelectionLSI(
    x = SCE, 
    useMatrix = "counts",
    preTfIdf = NULL, 
    binarize = T, 
    nDimensions = 5,
    outlierQuantile = c(0,1),
    varFeatures = 1000,
    seed = 1)
  
  it("Projects RSV identical to ArchR", {
    testthat::expect_identical(z, ArchRLSI$matSVD[rownames(z), ])
  })
  it("Projects RSV properly with calculated V", {
    testthat::expect_equal(z, d, tolerance = 1e-5)
  })
  it("Project RSV with outliers properly", {
    testthat::expect_equal(SingleCellExperiment::sampleFactors(reducedDim(SCE, "LSI"))[rownames(ArchRLSI$matSVD), ], ArchRLSI$matSVD, tolerance = 1e-4)
  })
  it("Select the same features on the first pass", {
    expect_all_true(all(LSIres1@metadata$selectedFeatures %in% archRTopFeats))
  })
  it("Calculates First Pass LSI similar to ArchR", {
    expect_equal(SingleCellExperiment::sampleFactors(LSIres1), first_pass_LSI$matSVD, tolerance = 1e-4)
  })
  it("clusters the same on the first pass", {
    expect_identical(clustF, first_pass_clust)
  })
  it("select the same features on the second pass", {
    expect_all_true(LSIres2@metadata$selectedFeatures %in% ArchRVarFeat)
  })
  it(".iterativePass() calculate same LSI as archR the second pass", {
    expect_equal(SingleCellExperiment::sampleFactors(LSIres2), second_pass_LSI$matSVD, tolerance = 1e-4)
  })
  it("iterativeLSI yields same result as ArchR",{
    expect_equal(SingleCellExperiment::sampleFactors(LSIres_all), second_pass_LSI$matSVD, tolerance = 1e-4)
    })
  it("select similar features to ArchR", {
    expect_all_true(LSIres_all@metadata$selectedFeatures %in% ArchRVarFeat)
  })
  it("adds iterativeLSI as reducDim properly",{
    SCE = addIterativeLSI(
      x = SCE, 
      name = "test",
      useMatrix = "counts",
      preTfIdf = NULL, 
      binarize = T, 
      nDimensions = 5,
      outlierQuantile = c(0,1),
      varFeatures = 1000)
    
    expect_true("test" %in% reducedDimNames(SCE))
    expect_true(inherits(reducedDim(SCE, "test"),"LinearEmbeddingMatrix"))
    expect_equal(nrow(reducedDim(SCE, "test")), ncol(SCE))
    expect_equal(ncol(reducedDim(SCE, "test")), 5)
    expect_equal(rownames(reducedDim(SCE, "test")), colnames(SCE))
    expect_equal(
      SingleCellExperiment::sampleFactors(reducedDim(SCE, "test")),
      second_pass_LSI$matSVD,
      tolerance = 1e-4
    )
    })
  
  describe("NMF method runs",{
    it("runs NMF",{
      SCE = addIterativeLSI(
        x = SCE, 
        name = "testNMF",
        useMatrix = "counts",
        preTfIdf = NULL, 
        binarize = T, 
        nDimensions = 5,
        outlierQuantile = c(0,1),
        varFeatures = 1000,
        method = "NMF")
      
      testthat::expect_true("testNMF" %in% SingleCellExperiment::reducedDimNames(SCE))
      })
    # it(,{})


    })
})
