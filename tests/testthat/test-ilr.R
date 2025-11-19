test_that("ilr_transform is isometric to legacy pivotCoord (no zeros)", {
  suppressWarnings(skip_if_not_installed("robCompositions"))
  set.seed(123)

  # simulate strictly positive counts (no zeros)
  X <- matrix(
    sample(1:100, size = 5 * 10, replace = TRUE),
    nrow = 10, ncol = 5
  )
  colnames(X) <- paste0("P", 1:5)

  # --- legacy version ---
  # robCompositions::pivotCoord uses its own ILR basis
  ILR_legacy <- robCompositions::pivotCoord(X)

  # --- new function ---
  ILR_new <- ilr_transform(X, zero.method = "none")

  # basic dimension check
  expect_equal(dim(ILR_legacy), dim(ILR_new))

  # --- geometry: pairwise distances must match (up to rotation) ---
  d_legacy <- as.matrix(dist(ILR_legacy))
  d_new    <- as.matrix(dist(ILR_new))

  expect_equal(d_legacy, d_new, tolerance = 1e-10, ignore_attr = TRUE)

  # --- geometry: PCA spectrum (eigenvalues of covariance) must match ---
  Sigma_legacy <- cov(scale(ILR_legacy, center = TRUE, scale = FALSE))
  Sigma_new    <- cov(scale(ILR_new,    center = TRUE, scale = FALSE))

  ev_legacy <- eigen(Sigma_legacy, symmetric = TRUE, only.values = TRUE)$values
  ev_new    <- eigen(Sigma_new,    symmetric = TRUE, only.values = TRUE)$values

  expect_equal(ev_legacy, ev_new, tolerance = 1e-10)
})


test_that("ilr_transform is isometric to legacy pivotCoord with cmultRepl", {
  suppressWarnings(skip_if_not_installed("robCompositions"))

  # allow zeros, use real data
  X <- as.matrix(ECU.MF[, 2:6])

  # --- legacy version ---
  X_cmult   <- zCompositions::cmultRepl(X)
  ILR_legacy <- robCompositions::pivotCoord(X_cmult)

  # --- new function ---
  ILR_new <- ilr_transform(X, zero.method = "cmultRepl")

  # basic dimension check
  expect_equal(dim(ILR_legacy), dim(ILR_new))

  # --- geometry: pairwise distances must match (up to rotation) ---
  d_legacy <- as.matrix(dist(ILR_legacy))
  d_new    <- as.matrix(dist(ILR_new))

  expect_equal(d_legacy, d_new, tolerance = 1e-10, ignore_attr = TRUE)

  # --- geometry: PCA spectrum must match ---
  Sigma_legacy <- cov(scale(ILR_legacy, center = TRUE, scale = FALSE))
  Sigma_new    <- cov(scale(ILR_new,    center = TRUE, scale = FALSE))

  ev_legacy <- eigen(Sigma_legacy, symmetric = TRUE, only.values = TRUE)$values
  ev_new    <- eigen(Sigma_new,    symmetric = TRUE, only.values = TRUE)$values

  expect_equal(ev_legacy, ev_new, tolerance = 1e-10)
})
