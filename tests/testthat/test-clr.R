test_that("clr_transform matches legacy CLR code (no zeros)", {
  set.seed(123)

  # simulate strictly positive counts (no zeros)
  X <- matrix(
    sample(1:100, size = 5 * 10, replace = TRUE),
    nrow = 10, ncol = 5
  )
  colnames(X) <- paste0("P", 1:5)

  # --- legacy version ---

  Co_df <- as.data.frame(X)
  Co_df$sum <- rowSums(Co_df)
  for (i in 1:5) {
    Co_df[, i] <- Co_df[, i] / Co_df[, 6]
  }

  # CLR on first 5 parts, using explicit geometric mean column
  clr_df <- Co_df[, 1:5]
  clr_df$geo_mean <- apply(clr_df, 1, function(x) exp(mean(log(x))))
  for (i in 1:5) {
    clr_df[, i] <- log(clr_df[, i] / clr_df[, 6])
  }
  clr_legacy <- as.matrix(clr_df[, 1:5])

  # --- new function ---
  clr_new <- clr_transform(X, zero.method = "none")

  # compare (only first 5 cols for clarity — function returns all 5 anyway)
  expect_equal(clr_new, clr_legacy, tolerance = 1e-10)
})

test_that("clr_transform matches legacy CLR code with cmultRepl", {
  set.seed(123)

  # allow zeros, use real data

  # --- legacy function ---
  X <- as.matrix(ECU.MF[, 2:6])
  Co_df <- zCompositions::cmultRepl(X)

  Co_df$sum <- rowSums(Co_df)
  for (i in 1:5) {
    Co_df[,i] <- Co_df[,i] / Co_df[,6]
  }

  clr_df <- Co_df[,1:5]
  clr_df$geo_mean <- apply(clr_df, 1, function(x){exp(mean(log(x)))})
  for (i in 1:5) {
    clr_df[,i] <- log(clr_df[,i] / clr_df[,6])
  }
  clr_legacy <- as.matrix(clr_df[,1:5])

  # --- new function ---
  clr_new <- clr_transform(X, zero.method = "cmultRepl")

  # both should match for the 5 parts
  expect_equal(clr_new, clr_legacy, tolerance = 1e-10)
})

