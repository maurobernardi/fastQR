## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = TRUE, message = FALSE,
                     error = FALSE, fig.width = 8, fig.height = 4.8,
                     fig.align = "center", out.width = "100%", dpi = 160, cache = FALSE)
figure_helpers <- if (file.exists("../scripts/figure-references.R"))
  "../scripts/figure-references.R" else "scripts/figure-references.R"
source(figure_helpers, local = knitr::knit_global())
options(width = 90)
stopifnot(requireNamespace("fastQR", quietly = TRUE),
          requireNamespace("ggplot2", quietly = TRUE))


## ----decomp-data------------------------------------------------------------------------
set.seed(2026)
n <- 120L
p <- 6L
X <- matrix(rnorm(n * p), n, p)
thin <- fastQR::qr(X, type = "householder", complete = FALSE)
full <- fastQR::qr(X, type = "householder", complete = TRUE)
knitr::kable(data.frame(
  representation = c("Thin", "Full"),
  Q_rows = c(nrow(thin$Q), nrow(full$Q)),
  Q_columns = c(ncol(thin$Q), ncol(full$Q)),
  R_rows = c(nrow(thin$R), nrow(full$R)),
  R_columns = c(ncol(thin$R), ncol(full$R))))


## ----decomp-validation------------------------------------------------------------------
relative_error <- function(actual, reference) {
  norm(as.matrix(actual - reference), "F") /
    max(norm(as.matrix(reference), "F"), .Machine$double.eps)
}
checks <- c(
  thin_reconstruction = relative_error(thin$Q %*% thin$R, X),
  full_reconstruction = relative_error(full$Q %*% full$R, X),
  thin_orthogonality = relative_error(crossprod(thin$Q), diag(p)),
  full_orthogonality = relative_error(crossprod(full$Q), diag(n)),
  Gram_identity = relative_error(crossprod(thin$R), crossprod(X)))
knitr::kable(data.frame(check = names(checks), relative_error = unname(checks)),
             digits = 16)
stopifnot(max(checks) < 1e-10)


## ----direct-qr-estimation---------------------------------------------------------------
demo_y <- drop(X %*% seq_len(p)) + sin(seq_len(n))
demo_projection <- crossprod(thin$Q, demo_y)
demo_coefficients <- backsolve(thin$R, demo_projection)
demo_fitted <- drop(X %*% demo_coefficients)
demo_residuals <- demo_y - demo_fitted
qr_derivation_checks <- c(
  coefficient_agreement = max(abs(demo_coefficients - lm.fit(X, demo_y)$coefficients)),
  fitted_projection = relative_error(demo_fitted, thin$Q %*% demo_projection),
  residual_orthogonality = max(abs(crossprod(thin$Q, demo_residuals))) /
    max(1, sqrt(sum(demo_y^2))))
stopifnot(max(qr_derivation_checks) < 1e-10)
knitr::kable(data.frame(check = names(qr_derivation_checks),
                        error = unname(qr_derivation_checks)), digits = 12)


## ----pivoting---------------------------------------------------------------------------
packed <- fastQR::qr_fast(X, pivot = TRUE)
Q <- fastQR::qr_Q(packed$qr, tau = packed$qraux,
                  rank = packed$rank, complete = FALSE)
R <- fastQR::qr_R(packed$qr, rank = NULL, complete = FALSE)
P <- fastQR::qr_pivot2perm(packed$pivot)
pivot_checks <- c(
  reconstruction = relative_error(Q %*% R %*% P, X),
  permuted_Gram = relative_error(crossprod(R), crossprod(X %*% t(P))))
knitr::kable(data.frame(check = names(pivot_checks), error = pivot_checks),
             digits = 16, row.names = FALSE)
stopifnot(max(pivot_checks) < 1e-10)
pivoted_coefficients <- backsolve(R, crossprod(Q, demo_y))
original_order_coefficients <- drop(t(P) %*% pivoted_coefficients)
pivot_derivation_error <- max(abs(original_order_coefficients - demo_coefficients))
stopifnot(pivot_derivation_error < 1e-10)


## ----triangular-and-solve---------------------------------------------------------------
R_qr <- fastQR::qrchol(X)
R_only <- fastQR::rchol(X)
stopifnot(relative_error(crossprod(R_qr), crossprod(X)) < 1e-10,
          relative_error(crossprod(R_only), crossprod(X)) < 1e-10)

A <- crossprod(X) + diag(p)
b <- seq_len(p)
solution <- fastQR::qrsolve(A, b, type = "QR")
relative_error(A %*% solution, b)
stopifnot(relative_error(A %*% solution, b) < 1e-10)


## ----update-rows------------------------------------------------------------------------
set.seed(2027)
n <- 120L
p <- 6L
X <- matrix(rnorm(n * p), n, p)
initial <- fastQR::qr(X, complete = TRUE)
U <- matrix(rnorm(p * 3L), p, 3L)
X_added <- rbind(X, t(U))
added <- fastQR::qrupdate(
  Q = initial$Q, R = initial$R, k = n + 1L, U = U,
  type = "row", fast = FALSE, complete = TRUE)
removed <- fastQR::qrdowndate(
  Q = added$Q, R = added$R, k = n + 1L, m = ncol(U),
  type = "row", fast = FALSE, complete = TRUE)
relative_error <- function(a, b) {
  norm(as.matrix(a - b), "F") / max(norm(as.matrix(b), "F"), .Machine$double.eps)
}
row_checks <- c(
  after_insertion = relative_error(added$Q %*% added$R, X_added),
  after_deletion = relative_error(removed$Q %*% removed$R, X))
knitr::kable(data.frame(operation = names(row_checks), error = row_checks),
             digits = 16, row.names = FALSE)
stopifnot(max(row_checks) < 1e-10)


## ----update-columns---------------------------------------------------------------------
V <- matrix(rnorm(n * 2L), n, 2L)
k <- 3L
X_wider <- cbind(X[, 1:2, drop = FALSE], V, X[, 3:p, drop = FALSE])
wider <- fastQR::qrupdate(
  initial$Q, initial$R, k = k, U = V, type = "column",
  fast = FALSE, complete = TRUE)
narrower <- fastQR::qrdowndate(
  wider$Q, wider$R, k = k, m = ncol(V), type = "column",
  fast = FALSE, complete = TRUE)
stopifnot(relative_error(wider$Q %*% wider$R, X_wider) < 1e-10,
          relative_error(narrower$Q %*% narrower$R, X) < 1e-10)


## ----update-r-only----------------------------------------------------------------------
R_initial <- initial$R[seq_len(p), , drop = FALSE]
R_added <- fastQR::rupdate(X, R_initial, U = U, type = "row", fast = FALSE)
R_removed <- fastQR::rdowndate(R_added, U = U, type = "row", fast = FALSE)

# R-only column updates append new columns to the right.
R_wider <- fastQR::rupdate(X, R_initial, U = V, type = "column", fast = FALSE)
R_narrower <- fastQR::rdowndate(R_wider, k = p + 1L, m = ncol(V),
                               type = "column", fast = FALSE)
r_checks <- c(
  add_rows = relative_error(crossprod(R_added), crossprod(X_added)),
  remove_rows = relative_error(crossprod(R_removed), crossprod(X)),
  add_columns = relative_error(crossprod(R_wider), crossprod(cbind(X, V))),
  remove_columns = relative_error(crossprod(R_narrower), crossprod(X)))
knitr::kable(data.frame(operation = names(r_checks), error = r_checks),
             digits = 16, row.names = FALSE)
stopifnot(max(r_checks) < 1e-10)


## ----sliding-window---------------------------------------------------------------------
set.seed(2028)
window <- 100L
steps <- 30L
p <- 5L
X_stream <- cbind(1, matrix(rnorm((window + steps) * (p - 1L)),
                            window + steps, p - 1L))
beta_true <- c(1, 2, -1, 0, .5)
y_stream <- drop(X_stream %*% beta_true + rnorm(window + steps, sd = .7))
X_window <- X_stream[seq_len(window), , drop = FALSE]
y_window <- y_stream[seq_len(window)]
R_window <- fastQR::qr(X_window, complete = FALSE)$R
h <- crossprod(X_window, y_window)
coef_update <- coef_fresh <- matrix(NA_real_, steps, p)
gram_error <- numeric(steps)
for (s in seq_len(steps)) {
  old_x <- X_window[1, , drop = FALSE]
  old_y <- y_window[1]
  new_x <- X_stream[window + s, , drop = FALSE]
  new_y <- y_stream[window + s]
  R_minus <- fastQR::rdowndate(R_window, U = t(old_x), type = "row")
  X_minus <- X_window[-1, , drop = FALSE]
  R_window <- fastQR::rupdate(X_minus, R_minus, U = t(new_x), type = "row")
  h <- h - crossprod(old_x, old_y) + crossprod(new_x, new_y)
  X_window <- rbind(X_minus, new_x)
  y_window <- c(y_window[-1], new_y)
  coef_update[s, ] <- backsolve(R_window, forwardsolve(t(R_window), h))
  coef_fresh[s, ] <- lm.fit(X_window, y_window)$coefficients
  gram_error[s] <- relative_error(crossprod(R_window), crossprod(X_window))
}
max_coefficient_error <- max(abs(coef_update - coef_fresh))
stopifnot(max_coefficient_error < 1e-9, max(gram_error) < 1e-10)
knitr::kable(data.frame(
  check = c("Maximum coefficient difference", "Maximum relative Gram error"),
  value = c(max_coefficient_error, max(gram_error))), digits = 16)


## ----sliding-path, fig.cap=guide_caption("sliding-path", "Coefficient estimates as the window advances. Solid lines use updated factors; dashed lines are fresh least-squares fits. The two computations agree to numerical precision.")----
path_data <- rbind(
  data.frame(step = rep(seq_len(steps), p), coefficient = factor(rep(seq_len(p), each = steps)),
             estimate = as.vector(coef_update), method = "Updated factors"),
  data.frame(step = rep(seq_len(steps), p), coefficient = factor(rep(seq_len(p), each = steps)),
             estimate = as.vector(coef_fresh), method = "Fresh fit"))
ggplot2::ggplot(path_data, ggplot2::aes(step, estimate, colour = coefficient,
                                      linetype = method)) +
  ggplot2::geom_line(linewidth = .8) +
  ggplot2::scale_linetype_manual(values = c("Fresh fit" = "dashed", "Updated factors" = "solid")) +
  ggplot2::labs(x = "Window replacement", y = "Coefficient estimate", colour = "Coefficient",
                linetype = NULL) + ggplot2::theme_minimal(base_size = 12)


## ----regression-data--------------------------------------------------------------------
set.seed(2029)
n <- 160L
p <- 6L
n_test <- 400L
X <- matrix(rnorm(n * p), n, p)
X_test <- matrix(rnorm(n_test * p), n_test, p)
beta_true <- c(2, -1.5, .75, 0, 0, 0)
y <- drop(1 + X %*% beta_true + rnorm(n, sd = 1.5))
y_test <- drop(1 + X_test %*% beta_true + rnorm(n_test, sd = 1.5))
design <- cbind(intercept = 1, X)
design_test <- cbind(intercept = 1, X_test)
ols <- fastQR::qrls(y, design, X_test = design_test, type = "QR")
ols_factors <- fastQR::qr(design, complete = FALSE)
ols_projection <- crossprod(ols_factors$Q, y)
ols_explicit <- backsolve(ols_factors$R, ols_projection)
ols_derivation_error <- max(abs(ols$coeff - ols_explicit))
stopifnot(ols_derivation_error < 1e-10,
          max(abs(ols$coeff - lm.fit(design, y)$coefficients)) < 1e-10)
knitr::kable(data.frame(
  term = c("Intercept", paste0("x", seq_len(p))),
  truth = c(1, beta_true), estimate = as.numeric(ols$coeff)), digits = 4)
ols_test_mse <- mean((y_test - drop(design_test %*% ols$coeff))^2)
ols_test_mse


## ----regression-predictions, fig.cap=guide_caption("regression-predictions", "Independent test responses versus least-squares predictions. The diagonal is a reference line, not a fitted trend.")----
ggplot2::ggplot(data.frame(observed = y_test,
                          predicted = drop(design_test %*% ols$coeff)),
                ggplot2::aes(observed, predicted)) +
  ggplot2::geom_point(alpha = .4, colour = "#b5121b") +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#526070") +
  ggplot2::coord_equal() + ggplot2::theme_minimal(base_size = 12)


## ----ridge-path-------------------------------------------------------------------------
x_mean <- colMeans(X)
x_scale <- apply(X, 2, sd)
y_mean <- mean(y)
Z <- sweep(sweep(X, 2, x_mean, "-"), 2, x_scale, "/")
y_centered <- y - y_mean
lambda <- exp(seq(log(1e-4), log(10), length.out = 50L))
path_standardized <- vapply(lambda, function(lam) {
  as.numeric(fastQR::qrridge(y_centered, Z, lambda = lam, type = "QR")$coeff)
}, numeric(p))
path_triangular <- vapply(lambda, function(lam) {
  as.numeric(fastQR::qrridge(y_centered, Z, lambda = lam, type = "R")$coeff)
}, numeric(p))
stopifnot(max(abs(path_standardized - path_triangular)) < 1e-10)
path_slopes <- sweep(path_standardized, 1, x_scale, "/")
path_intercept <- y_mean - drop(crossprod(x_mean, path_slopes))

# Independent augmented-design QR calculation at one grid value.
j <- 20L
reference <- qr.coef(base::qr(rbind(Z, sqrt(n * lambda[j]) * diag(p))),
                     c(y_centered, rep(0, p)))
stopifnot(max(abs(path_standardized[, j] - reference)) < 1e-10)


## ----ridge-triangular-steps-------------------------------------------------------------
ridge_R_fit <- fastQR::qrridge(y_centered, Z, lambda[j], type = "R")
ridge_rhs <- crossprod(Z, y_centered)
ridge_coordinates <- forwardsolve(t(ridge_R_fit$R), ridge_rhs)
ridge_explicit <- backsolve(ridge_R_fit$R, ridge_coordinates)
ridge_manual_error <- max(abs(ridge_explicit - path_triangular[, j]))
stopifnot(ridge_manual_error < 1e-10)


## ----ridge-coefficients, fig.cap=guide_caption("ridge-coefficients", "Ridge slopes on their original predictor scale, plotted against a logarithmic penalty grid. The coefficients change continuously with the penalty; ridge does not generally create an exactly sparse fit.")----
path_table <- data.frame(lambda = rep(lambda, each = p),
                         coefficient = factor(rep(seq_len(p), length(lambda))),
                         estimate = as.vector(path_slopes))
ggplot2::ggplot(path_table, ggplot2::aes(lambda, estimate, colour = coefficient)) +
  ggplot2::geom_line(linewidth = .8) + ggplot2::scale_x_log10() +
  ggplot2::labs(x = expression(lambda), y = "Coefficient estimate", colour = "Predictor") +
  ggplot2::theme_minimal(base_size = 12)


## ----ridge-cv---------------------------------------------------------------------------
set.seed(2030)
K <- 5L
fold <- sample(rep(seq_len(K), length.out = n))
squared_error <- matrix(NA_real_, n, length(lambda))
for (k in seq_len(K)) {
  train <- fold != k
  valid <- !train
  mu <- colMeans(X[train, , drop = FALSE])
  sc <- apply(X[train, , drop = FALSE], 2, sd)
  ybar <- mean(y[train])
  Z_train <- sweep(sweep(X[train, , drop = FALSE], 2, mu, "-"), 2, sc, "/")
  Z_valid <- sweep(sweep(X[valid, , drop = FALSE], 2, mu, "-"), 2, sc, "/")
  for (j in seq_along(lambda)) {
    b <- fastQR::qrridge(y[train] - ybar, Z_train,
                         lambda = lambda[j], type = "QR")$coeff
    prediction <- ybar + drop(Z_valid %*% b)
    squared_error[valid, j] <- (y[valid] - prediction)^2
  }
}
cv_mse <- colMeans(squared_error)
best <- which.min(cv_mse)
lambda_selected <- lambda[best]
selected_prediction <- path_intercept[best] + drop(X_test %*% path_slopes[, best])
knitr::kable(data.frame(
  method = c("Least squares", "CV-selected ridge"),
  test_MSE = c(ols_test_mse, mean((y_test - selected_prediction)^2))), digits = 4)
lambda_selected


## ----ridge-cv-plot, fig.cap=guide_caption("ridge-cv-plot", "Validation mean squared error along the penalty grid. The dashed line marks the minimum on this grid.")----
ggplot2::ggplot(data.frame(lambda, MSE = cv_mse), ggplot2::aes(lambda, MSE)) +
  ggplot2::geom_line(colour = "#b5121b", linewidth = .9) +
  ggplot2::geom_vline(xintercept = lambda_selected, linetype = "dashed", colour = "#526070") +
  ggplot2::scale_x_log10() + ggplot2::labs(x = expression(lambda), y = "Cross-validated MSE") +
  ggplot2::theme_minimal(base_size = 12)


## ----ridge-df, fig.cap=guide_caption("ridge-df", "Effective degrees of freedom for each fixed ridge penalty, including the unpenalized intercept. The dashed line marks the cross-validated penalty; the curve describes fixed-penalty fits.")----
d2 <- svd(Z, nu = 0, nv = 0)$d^2
df_path <- vapply(lambda, function(lam) 1 + sum(d2 / (d2 + n * lam)), numeric(1))
stopifnot(all(diff(df_path) <= 0), all(df_path >= 1), all(df_path <= p + 1))
ggplot2::ggplot(data.frame(lambda, df = df_path), ggplot2::aes(lambda, df)) +
  ggplot2::geom_line(colour = "#b5121b", linewidth = .9) + ggplot2::scale_x_log10() +
  ggplot2::geom_vline(xintercept = lambda_selected, linetype = "dashed") +
  ggplot2::labs(x = expression(lambda), y = "Effective degrees of freedom") +
  ggplot2::theme_minimal(base_size = 12)


## ----multivariate-regression------------------------------------------------------------
set.seed(2031)
Y <- cbind(y, drop(-.5 + X %*% c(0, 1, 0, -1, 0, .5) + rnorm(n)))
multivariate <- fastQR::qrmls(Y, design, type = "QR")
reference <- qr.coef(base::qr(design), Y)
multi_direct <- backsolve(ols_factors$R, crossprod(ols_factors$Q, Y))
multi_direct_error <- max(abs(multivariate$coeff - multi_direct))
stopifnot(multi_direct_error < 1e-10,
          max(abs(multivariate$coeff - reference)) < 1e-10)
knitr::kable(as.matrix(multivariate$coeff), digits = 4,
             col.names = c("Response 1", "Response 2"))


## ----native-ridge-cv--------------------------------------------------------------------
native_qr <- fastQR::qrridge_cv(y, design, lambda, k = 5L, seed = 2032L,
                                X_test = design_test, type = "QR")
native_r <- fastQR::qrridge_cv(y, design, lambda, k = 5L, seed = 2032L,
                               X_test = design_test, type = "R")
stopifnot(identical(native_qr$foldid, native_r$foldid),
          max(abs(native_qr$mse - native_r$mse)) < 1e-10,
          max(abs(native_qr$coeff - native_r$coeff)) < 1e-10)

# Reconstruct the entire validation curve independently.
native_reference_mse <- vapply(lambda, function(lam) {
  mean(vapply(sort(unique(native_qr$foldid)), function(k) {
    train <- native_qr$foldid != k
    valid <- !train
    nt <- sum(train)
    A <- rbind(design[train, , drop = FALSE],
               sqrt(nt * lam) * diag(ncol(design)))
    b <- qr.coef(base::qr(A), c(y[train], rep(0, ncol(design))))
    mean((y[valid] - drop(design[valid, , drop = FALSE] %*% b))^2)
  }, numeric(1)))
}, numeric(1))
stopifnot(max(abs(as.numeric(native_qr$mse) - native_reference_mse)) < 1e-10)
knitr::kable(data.frame(
  backend = c("QR", "R"),
  selected_lambda = c(native_qr$lambda_min, native_r$lambda_min),
  validation_MSE = c(native_qr$mse_min, native_r$mse_min),
  test_MSE = c(mean((y_test - native_qr$predicted)^2),
               mean((y_test - native_r$predicted)^2))), digits = 5,
  caption = "Native cross-validation with a penalized intercept column.")


## ----native-multivariate-ridge----------------------------------------------------------
multi_qr <- fastQR::qrmridge_cv(Y, design, lambda, k = 5L, seed = 2032L,
                                X_test = design_test, type = "QR")
multi_r <- fastQR::qrmridge_cv(Y, design, lambda, k = 5L, seed = 2032L,
                               X_test = design_test, type = "R")
A <- rbind(design, sqrt(n * multi_qr$lambda_min) * diag(ncol(design)))
B_reference <- qr.coef(base::qr(A), rbind(Y, matrix(0, ncol(design), ncol(Y))))
stopifnot(max(abs(multi_qr$coeff - B_reference)) < 1e-10,
          max(abs(multi_r$coeff - B_reference)) < 1e-10,
          max(abs(multi_qr$mse - multi_r$mse)) < 1e-10,
          max(abs(multi_qr$predicted - design_test %*% B_reference)) < 1e-10,
          !anyDuplicated(names(multi_qr)), !anyDuplicated(names(multi_r)))
knitr::kable(data.frame(
  backend = c("QR", "R"),
  selected_lambda = c(multi_qr$lambda_min, multi_r$lambda_min),
  validation_MSE = c(multi_qr$mse_min, multi_r$mse_min),
  maximum_coefficient_error = c(max(abs(multi_qr$coeff - B_reference)),
                                max(abs(multi_r$coeff - B_reference)))),
  digits = 6, caption = "Multiple-response ridge checked against an independent QR fit.")


## ----ridge-evidence-setup, include=FALSE------------------------------------------------
evidence_script_root <- if (dir.exists("../scripts")) "../scripts" else "scripts"
source(file.path(evidence_script_root, "ridge-evidence.R"), local = knitr::knit_global())
source(file.path(evidence_script_root, "check-ridge-evidence.R"),
       local = new.env(parent = knitr::knit_global()))
source(file.path(evidence_script_root, "ridge-evidence-experiments.R"), local = knitr::knit_global())


## ----ridge-evidence-fit-----------------------------------------------------------------
evidence_X <- evidence_results$demo$X
evidence_y <- evidence_results$demo$y
evidence_fit <- evidence_results$demo$fit
knitr::kable(data.frame(lambda = evidence_fit$lambda,
                        alpha = evidence_fit$alpha,
                        sigma2 = evidence_fit$sigma2,
                        status = evidence_fit$status), digits = 5,
             caption = "Empirical Bayes penalty and marginal-model variance estimate for the reproducible example.")
# In a new analysis, after sourcing scripts/ridge-evidence.R:
# evidence_fit <- ridge_evidence_fit(X, y)


## ----ridge-evidence-native-check--------------------------------------------------------
if (is.finite(evidence_fit$lambda)) {
  evidence_native <- fastQR::qrridge(evidence_y, evidence_X,
                                    evidence_fit$lambda, type = "R")
  evidence_native_error <- max(abs(evidence_native$coeff - evidence_fit$coeff))
  stopifnot(evidence_native_error < 1e-10)
  evidence_native_error
}


## ----ridge-evidence-objective, fig.cap=guide_caption("ridge-evidence-objective", "Profile negative marginal log likelihood as a function of the normalized ridge penalty, after profiling the noise variance. A vertical line marks the selected penalty when it is finite. The criterion is shifted by a constant for display.")----
evidence_curve <- evidence_fit$curve
evidence_curve$relative_score <- evidence_curve$score - min(evidence_curve$score)
evidence_plot <- ggplot2::ggplot(evidence_curve, ggplot2::aes(lambda, relative_score)) +
  ggplot2::geom_line(colour = "#b5121b") + ggplot2::scale_x_log10() +
  ggplot2::labs(x = expression(lambda), y = "Profile criterion above displayed minimum") +
  ggplot2::theme_minimal(base_size = 11)
if (is.finite(evidence_fit$lambda)) evidence_plot <- evidence_plot +
  ggplot2::geom_vline(xintercept = evidence_fit$lambda, linetype = 2)
evidence_plot


## ----ridge-evidence-validation----------------------------------------------------------
evidence_check_table <- evidence_results$checks
evidence_check_table[, -1] <- lapply(evidence_check_table[, -1],
  function(x) formatC(x, format = "e", digits = 2))
knitr::kable(evidence_check_table, row.names = FALSE,
 caption = "Agreement of the spectral marginal-likelihood calculation with dense-covariance and augmented-QR references.")


## ----ridge-evidence-comparison----------------------------------------------------------
knitr::kable(evidence_summary, digits = 4, row.names = FALSE,
 caption = "Median test-response MSE, mean-prediction MSE and selection time across five simulated training samples per setting. Timings include the necessary decompositions and coefficient recovery.")


## ----ridge-evidence-time-ratio, include=FALSE-------------------------------------------
evidence_time_pairs <- merge(
 subset(evidence_summary, method == "CV", select = c(n, p, signal, selection_ms)),
 subset(evidence_summary, method == "evidence", select = c(n, p, signal, selection_ms)),
 by = c("n", "p", "signal"), suffixes = c("_CV", "_marginal"))
evidence_time_ratio <- evidence_time_pairs$selection_ms_CV / evidence_time_pairs$selection_ms_marginal


## ----ridge-evidence-prediction, fig.cap=guide_caption("ridge-evidence-prediction", "Independent-test prediction error for empirical Bayes marginal likelihood and five-fold cross-validation. Each point corresponds to one simulated training sample; lines pair the selectors on the same sample. Panels distinguish design dimensions and coefficient-generating scenarios.")----
evidence_sim_plot <- evidence_results$simulations
evidence_sim_plot$design <- paste0("n = ", evidence_sim_plot$n, ", p = ", evidence_sim_plot$p)
ggplot2::ggplot(evidence_sim_plot,
 ggplot2::aes(method, test_mse, group = repetition)) +
 ggplot2::geom_line(colour = "grey70") + ggplot2::geom_point(colour = "#b5121b") +
 ggplot2::facet_grid(signal ~ design, scales = "free_y") +
 ggplot2::labs(x = NULL, y = "Test-response mean squared error") +
 ggplot2::theme_minimal(base_size = 11)


## ----ridge-evidence-boundaries----------------------------------------------------------
evidence_status <- aggregate(list(count = rep(1L, nrow(evidence_results$simulations))),
 evidence_results$simulations[, c("n", "p", "signal", "method", "boundary")], sum)
knitr::kable(evidence_status, row.names = FALSE,
 caption = "Observed selection statuses, including finite search boundaries and the null-model limit.")


## ----ridge-wide-example-----------------------------------------------------------------
set.seed(312)
n_wide <- 40L
p_wide <- 80L
X_wide <- matrix(rnorm(n_wide * p_wide), n_wide, p_wide)
y_wide <- rnorm(n_wide)
alpha_wide <- 8
lambda_wide <- alpha_wide / n_wide
wide_QR <- fastQR::qrridge(y_wide, X_wide, lambda_wide, type = "QR")
wide_R <- fastQR::qrridge(y_wide, X_wide, lambda_wide, type = "R")
A_wide <- rbind(X_wide, sqrt(alpha_wide) * diag(p_wide))
f_wide <- fastQR::qr_lapack(A_wide, output = "thin")
beta_reference <- backsolve(
  f_wide$R, crossprod(f_wide$Q, c(y_wide, rep(0, p_wide))))
c(QR_error = max(abs(wide_QR$coeff - beta_reference)),
  R_error = max(abs(wide_R$coeff - beta_reference)))
stopifnot(max(abs(wide_QR$coeff - beta_reference)) < 1e-10,
          max(abs(wide_R$coeff - beta_reference)) < 1e-10)


## ----ridge-wide-add-rows----------------------------------------------------------------
set.seed(313)
new_X <- matrix(rnorm(5L * p_wide), 5L, p_wide)
new_y <- rnorm(5L)
R_before <- fastQR::qr_lapack(A_wide, output = "R")
H_before <- crossprod(X_wide, y_wide)
R_after <- fastQR::rupdate(A_wide, R_before, U = t(new_X),
                          type = "row", fast = FALSE)
H_after <- H_before + crossprod(new_X, new_y)
beta_after <- backsolve(R_after, forwardsolve(t(R_after), H_after))
fit_after <- fastQR::qrridge(c(y_wide, new_y), rbind(X_wide, new_X),
                            lambda = alpha_wide / (n_wide + 5L), type = "R")
stopifnot(max(abs(beta_after - fit_after$coeff)) < 1e-10)


## ----ridge-wide-add-columns-------------------------------------------------------------
set.seed(314)
k_new <- 3L
W_new <- matrix(rnorm(n_wide * k_new), n_wide, k_new)
A_padded <- rbind(A_wide, matrix(0, k_new, p_wide))
W_augmented <- rbind(W_new, matrix(0, p_wide, k_new),
                     sqrt(alpha_wide) * diag(k_new))
R_enlarged <- fastQR::rupdate(A_padded, R_before, U = W_augmented,
                             type = "column", fast = FALSE)
H_enlarged <- rbind(H_before, crossprod(W_new, y_wide))
beta_enlarged <- backsolve(R_enlarged,
                           forwardsolve(t(R_enlarged), H_enlarged))
fit_enlarged <- fastQR::qrridge(y_wide, cbind(X_wide, W_new),
                               lambda_wide, type = "R")
stopifnot(max(abs(beta_enlarged - fit_enlarged$coeff)) < 1e-10)


## ----ridge-update-experiments, include=FALSE--------------------------------------------
ridge_script <- if (file.exists("../scripts/ridge-update-experiments.R"))
  "../scripts/ridge-update-experiments.R" else "scripts/ridge-update-experiments.R"
source(ridge_script, local = knitr::knit_global())


## ----ridge-update-tables----------------------------------------------------------------
knitr::kable(ridge_comparison, digits = 3, row.names = FALSE,
 caption = "Median milliseconds per fit after a single design change. Speedup is recomputation time divided by updating time; values above one favor updating.")
ridge_validation_display <- ridge_validation
ridge_validation_display[, 4:7] <- lapply(ridge_validation_display[, 4:7],
  function(x) formatC(x, format = "e", digits = 2))
knitr::kable(ridge_validation_display, row.names = FALSE,
 caption = "Relative coefficient and prediction errors against QR of the augmented design with an orthogonally transformed response. Errors are scaled by the larger of one and the reference norm.")


## ----ridge-growing-time, fig.cap=guide_caption("ridge-growing-time", "Cumulative elapsed time while adding observations or predictors at fixed total penalty. Both workflows include their initial factorization. Each panel shows one growing-design trajectory, rather than repeated timing estimates.")----
ridge_time_plot <- rbind(
 transform(ridge_path, method = "Updating", milliseconds = 1000 * update_cumulative),
 transform(ridge_path, method = "Recomputing", milliseconds = 1000 * recompute_cumulative))
ggplot2::ggplot(ridge_time_plot,
 ggplot2::aes(step, milliseconds, colour = method)) +
 ggplot2::geom_line() + ggplot2::geom_point(size = 1.5) +
 ggplot2::facet_wrap(~operation, scales = "free_y") +
 ggplot2::labs(x = "Update step", y = "Cumulative milliseconds", colour = NULL) +
 ggplot2::theme_minimal(base_size = 11) +
 ggplot2::theme(legend.position = "bottom")


## ----ridge-growing-error, fig.cap=guide_caption("ridge-growing-error", "Coefficient and prediction errors along the growing-design sequences. Each reference is recomputed from the current augmented design by QR and orthogonal transformation of the response. Values below machine precision are clipped only for logarithmic display.")----
ridge_error_plot <- rbind(
 transform(ridge_path, quantity = "Coefficients", error = coefficient_error),
 transform(ridge_path, quantity = "Predictions", error = prediction_error))
ggplot2::ggplot(ridge_error_plot,
 ggplot2::aes(step, pmax(error, .Machine$double.eps), colour = quantity)) +
 ggplot2::geom_line() + ggplot2::geom_point(size = 1.5) +
 ggplot2::scale_y_log10() + ggplot2::facet_wrap(~operation) +
 ggplot2::labs(x = "Update step", y = "Scaled error (log scale)", colour = NULL) +
 ggplot2::theme_minimal(base_size = 11) +
 ggplot2::theme(legend.position = "bottom")


## ----conditioning-----------------------------------------------------------------------
set.seed(2032)
n <- 160L
p <- 12L
U <- qr.Q(base::qr(matrix(rnorm(n * p), n, p)))
V <- qr.Q(base::qr(matrix(rnorm(p * p), p, p)))
beta <- seq(-1, 1, length.out = p)
condition_numbers <- 10^(0:8)
accuracy <- do.call(rbind, lapply(condition_numbers, function(kappa_target) {
  singular_values <- exp(seq(0, -log(kappa_target), length.out = p))
  X <- sweep(U, 2, singular_values, "*") %*% t(V)
  y <- drop(X %*% beta)
  f <- fastQR::qr(X, complete = FALSE)
  projected_response <- crossprod(f$Q, y)  # p coordinates in the Q basis
  b_qr <- backsolve(f$R, projected_response) # solve R beta = Q' y
  b_normal <- tryCatch(solve(crossprod(X), crossprod(X, y)),
                        error = function(e) rep(NA_real_, p))
  data.frame(condition = kappa_target,
    reconstruction = norm(X - f$Q %*% f$R, "F") / norm(X, "F"),
    QR_coefficients = sqrt(sum((b_qr - beta)^2) / sum(beta^2)),
    normal_equations = sqrt(sum((b_normal - beta)^2) / sum(beta^2)))
}))
knitr::kable(accuracy, digits = 12)
stopifnot(max(accuracy$reconstruction) < 1e-10)


## ----conditioning-plot, fig.cap=guide_caption("conditioning-plot", "Relative reconstruction and coefficient errors versus the prescribed condition number. Missing normal-equation values indicate that the reference solver rejected the system as numerically singular.")----
accuracy_long <- rbind(
  data.frame(condition = accuracy$condition, error = accuracy$reconstruction, metric = "QR reconstruction"),
  data.frame(condition = accuracy$condition, error = accuracy$QR_coefficients, metric = "QR coefficient solve"),
  data.frame(condition = accuracy$condition, error = accuracy$normal_equations, metric = "Normal-equation solve"))
ggplot2::ggplot(accuracy_long, ggplot2::aes(condition, pmax(error, .Machine$double.eps), colour = metric)) +
  ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(na.rm = TRUE) +
  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::labs(x = expression(kappa[2](X)), y = "Relative error", colour = NULL) +
  ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom")


## ----timing-----------------------------------------------------------------------------
set.seed(2033)
repetitions <- 5L
sizes <- c(500L, 1000L, 2000L)
p <- 40L
m <- 10L
elapsed_per_call <- function(fun, batch_size) {
  unname(system.time(for (i in seq_len(batch_size)) fun())["elapsed"]) / batch_size
}
timing <- do.call(rbind, lapply(sizes, function(n) {
  X <- matrix(rnorm(n * p), n, p)
  U <- matrix(rnorm(p * m), p, m)
  X_new <- rbind(X, t(U))
  R <- fastQR::qr(X, complete = FALSE)$R
  methods <- list(
    "R update" = function() fastQR::rupdate(X, R, U, type = "row", fast = FALSE),
    "Fresh thin QR" = function() fastQR::qr(X_new, complete = FALSE),
    "Fresh R only" = function() fastQR::rchol(X_new))
  updated <- methods[[1]]()
  gram_error <- norm(crossprod(updated) - crossprod(X_new), "F") / norm(crossprod(X_new), "F")
  stopifnot(gram_error < 1e-10)
  invisible(lapply(methods, function(fun) fun()))
  batches <- vapply(methods, function(fun) {
    b <- 1L
    repeat {
      elapsed <- unname(system.time(for (i in seq_len(b)) fun())["elapsed"])
      if (elapsed >= .05 || b >= 4096L) break
      b <- b * 2L
    }
    b
  }, integer(1))
  do.call(rbind, lapply(seq_len(repetitions), function(rep) {
    # Randomize execution order within each repetition.
    order <- sample(names(methods))
    data.frame(n = n, p = p, added_rows = m, repetition = rep, method = order,
      batch_calls = unname(batches[order]),
      seconds = vapply(order, function(name)
        elapsed_per_call(methods[[name]], batches[[name]]), numeric(1)))
  }))
}))
timing_summary <- aggregate(seconds ~ n + method, timing, median)
knitr::kable(timing_summary, digits = 6,
             col.names = c("Original rows", "Method", "Median seconds per call"))


## ----timing-plot, fig.cap=guide_caption("timing-plot", sprintf("Measured time to obtain factors after adding ten observations. Initial factors are assumed available for the update. Points and boxes summarize five repetitions; each repetition averages a calibrated batch of calls. Each design has %d columns and receives %d new rows.", p, m))----
ggplot2::ggplot(timing, ggplot2::aes(factor(n), seconds * 1000, colour = method)) +
  ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = .8), outlier.shape = NA) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = .8), size = 1.6) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = "Original number of observations", y = "Milliseconds per call (log scale)", colour = NULL) +
  ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom")


## ----benchmark-helpers------------------------------------------------------------------
full_benchmarks <- isTRUE(getOption("fastQR.guide.full_benchmarks", FALSE))
benchmark_reps <- if (full_benchmarks) 20L else 5L
benchmark_n <- if (full_benchmarks) 2000L else 500L
benchmark_p <- if (full_benchmarks) c(20L,40L,100L,200L,400L,800L,1200L,1500L,1800L) else c(20L,40L,80L)
measure_methods <- function(methods, repetitions) {
  invisible(lapply(methods, function(f) f()))
  batches <- vapply(methods, function(f) {
    b <- 1L
    repeat {
      seconds <- unname(system.time(for (i in seq_len(b)) f())["elapsed"])
      if (seconds >= .05 || b >= 4096L) break
      b <- b * 2L
    }
    b
  }, integer(1))
  do.call(rbind, lapply(seq_len(repetitions), function(rep) {
    order <- sample(names(methods))
    data.frame(repetition = rep, method = order, batch_calls = unname(batches[order]),
      seconds = vapply(order, function(k) elapsed_per_call(methods[[k]], batches[[k]]), numeric(1)))
  }))
}
set.seed(2041)
interface_n <- if (full_benchmarks) 1000L else 300L
interface_p <- if (full_benchmarks) 500L else 60L
interface_reps <- if (full_benchmarks) 10L else 5L
X_interface <- matrix(rnorm(interface_n * interface_p, 1), interface_n, interface_p)
interface_methods <- list(
  "base::qr (compact)" = function() base::qr(X_interface),
  "fastQR::qr (thin)" = function() fastQR::qr(X_interface, nb = 5, complete = FALSE),
  "qr_fast, no pivot" = function() fastQR::qr_fast(X_interface, pivot = FALSE),
  "qr_fast, pivot" = function() fastQR::qr_fast(X_interface, pivot = TRUE),
  "Matrix::qr (dense input)" = function() Matrix::qr(X_interface))
interface_timing <- measure_methods(interface_methods, interface_reps)
knitr::kable(aggregate(seconds ~ method, interface_timing, median), digits = 6,
  caption = sprintf("Median seconds per call: n = %d, p = %d, repetitions = %d.",
                    interface_n, interface_p, interface_reps))


## ----original-scaling-benchmarks--------------------------------------------------------
set.seed(2042)
benchmark_validation <- list()
scaling_timing <- do.call(rbind, lapply(benchmark_p, function(p) {
  n <- benchmark_n
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  G <- crossprod(X); h <- crossprod(X, y)
  block <- max(1L, p %/% 2L)
  b_reference <- base::qr.solve(G, h)
  errors <- numeric()
  for (nb in c(5L, block)) {
    f <- fastQR::qr(X, nb = nb, complete = FALSE)
    U <- fastQR::qrchol(X, nb = nb)
    b <- fastQR::qrsolve(A = G, b = h, nb = nb)
    errors <- c(errors, norm(X - f$Q %*% f$R, "F") / norm(X, "F"),
      norm(crossprod(U) - G, "F") / norm(G, "F"),
      sqrt(sum((b - b_reference)^2)) / max(1, sqrt(sum(b_reference^2))))
  }
  stopifnot(max(errors) < 1e-8)
  benchmark_validation[[as.character(p)]] <<- data.frame(p = p, max_relative_error = max(errors))
  workloads <- list(
    QR = list("fastQR, nb = 5" = function() fastQR::qr(X, nb = 5, complete = FALSE),
      "fastQR, nb = p/2" = function() fastQR::qr(X, nb = block, complete = FALSE),
      "base QR (compact)" = function() base::qr(X)),
    Cholesky = list("qrchol, nb = 5" = function() fastQR::qrchol(X, nb = 5),
      "qrchol, nb = p/2" = function() fastQR::qrchol(X, nb = block),
      "chol, Gram available" = function() base::chol(G),
      "crossprod + chol" = function() base::chol(crossprod(X))),
    Solve = list("qrsolve, nb = 5" = function() fastQR::qrsolve(G, h, nb = 5),
      "qrsolve, nb = p/2" = function() fastQR::qrsolve(G, h, nb = block),
      "base qr.solve" = function() base::qr.solve(G, h)))
  do.call(rbind, lapply(names(workloads), function(workload) {
    out <- measure_methods(workloads[[workload]], benchmark_reps)
    out$workload <- workload; out$n <- n; out$p <- p
    out
  }))
}))
benchmark_validation <- do.call(rbind, benchmark_validation)
knitr::kable(benchmark_validation, digits = 12,
             caption = "Largest relative validation error across the two block sizes and three workloads.")
knitr::kable(aggregate(seconds ~ workload + n + p + method, scaling_timing, median), digits = 6,
             caption = "Median seconds per call; dimensions and method names are reported explicitly.")
benchmark_plot <- function(workload_name) {
  dat <- subset(scaling_timing, workload == workload_name)
  ggplot2::ggplot(dat, ggplot2::aes(factor(p), seconds * 1000, colour = method)) +
    ggplot2::geom_boxplot(outlier.shape = NA, position = ggplot2::position_dodge(.8)) +
    ggplot2::geom_point(size = 1.3, position = ggplot2::position_dodge(.8)) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Number of columns", y = "Milliseconds per call (log scale)", colour = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}


## ----benchmark-reading, include=FALSE---------------------------------------------------
benchmark_reading <- function(workload_name) {
  z <- subset(scaling_timing, workload == workload_name)
  med <- aggregate(seconds ~ p + method, z, median)
  large <- subset(med, p == max(p))
  best <- large[which.min(large$seconds), ]
  sprintf("At the largest tested dimension (p = %d), the smallest median time is %.3g milliseconds per call for %s.",
          best$p, 1000 * best$seconds, best$method)
}


## ----restored-qr-plot, fig.cap=guide_caption("restored-qr-plot", sprintf("QR interface timings across column dimensions. Compact and explicit outputs are distinguished in the legend. Designs have %d rows; each setting uses %d timed repetitions.", benchmark_n, benchmark_reps))----
benchmark_plot("QR")


## ----restored-cholesky-plot, fig.cap=guide_caption("restored-cholesky-plot", sprintf("Triangular-factor timings, with separate baselines for an available Gram matrix and for its formation from X. Designs have %d rows; each setting uses %d timed repetitions.", benchmark_n, benchmark_reps))----
benchmark_plot("Cholesky")


## ----restored-solve-plot, fig.cap=guide_caption("restored-solve-plot", sprintf("Solution times for the same precomputed normal-equation system. Designs have %d rows; each setting uses %d timed repetitions.", benchmark_n, benchmark_reps))----
benchmark_plot("Solve")


## ----lapack-comparison-data, include=FALSE----------------------------------------------
lapack_root <- if (dir.exists("../results")) ".." else "."
initial_lapack <- read.csv(file.path(lapack_root, "results/bench_initial.csv"))
modifications_lapack <- read.csv(file.path(lapack_root, "results/bench_modifications.csv"))
lapack_cfg <- readRDS(file.path(lapack_root, "results/benchmark-settings.rds"))


## ----lapack-matched-tables--------------------------------------------------------------
initial_medians <- aggregate(seconds ~ n + p + workload + method, initial_lapack, median)
largest_thin <- subset(initial_medians, workload == "thin_unpivoted" & n == max(n) & p == max(p))
base_thin_seconds <- largest_thin$seconds[largest_thin$method == "LAPACK DGEQRF + DORGQR"]
largest_thin$milliseconds <- 1000 * largest_thin$seconds
largest_thin$LAPACK_over_method <- base_thin_seconds / largest_thin$seconds
knitr::kable(largest_thin[, c("n", "p", "method", "milliseconds", "LAPACK_over_method")],
             digits = 4, row.names = FALSE,
             caption = "Thin Q and R for the same largest design. A time ratio above one favors the named method over LAPACK.")

largest_R <- subset(initial_medians, workload == "R_unpivoted" & n == max(n) & p == max(p))
base_R_seconds <- largest_R$seconds[largest_R$method == "LAPACK DGEQRF + extract R"]
largest_R$milliseconds <- 1000 * largest_R$seconds
largest_R$LAPACK_over_method <- base_R_seconds / largest_R$seconds
knitr::kable(largest_R[, c("n", "p", "method", "milliseconds", "LAPACK_over_method")],
             digits = 4, row.names = FALSE,
             caption = "Triangular-factor-only requests on the same design, with no explicit Q returned.")


## ----lapack-thin-comparison, fig.cap=guide_caption("lapack-thin-comparison", sprintf("Thin-factor timings for a %d by %d design. Each point represents one calibrated repetition; elapsed time is shown on a logarithmic scale, with lower values indicating faster calls.", max(lapack_cfg$n), max(lapack_cfg$p)))----
plot_lapack <- subset(initial_lapack, workload == "thin_unpivoted" & n == max(n) & p == max(p))
ggplot2::ggplot(plot_lapack, ggplot2::aes(x = seconds * 1000, y = method)) +
  ggplot2::geom_boxplot(outlier.shape = NA, colour = "#526070") +
  ggplot2::geom_point(colour = "#b5121b", size = 1.5) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(x = "Milliseconds per call (log scale)", y = NULL) +
  ggplot2::theme_minimal(base_size = 11)


## ----lapack-update-tables---------------------------------------------------------------
modification_medians <- aggregate(seconds ~ operation + workload + method, modifications_lapack, median)
matched_updates <- do.call(rbind, lapply(split(modification_medians,
  interaction(modification_medians$operation, modification_medians$workload, drop = TRUE)), function(z) {
    lapack <- grepl("LAPACK", z$method)
    data.frame(operation = gsub("_", " ", z$operation[1]), workload = z$workload[1],
      fastQR_ms = 1000 * z$seconds[!lapack], LAPACK_ms = 1000 * z$seconds[lapack],
      LAPACK_over_fastQR = z$seconds[lapack] / z$seconds[!lapack])
}))
knitr::kable(matched_updates, digits = 4, row.names = FALSE,
             caption = sprintf("Single modifications of an initial %d by %d design, with batches of %d rows or columns. Ratios above one favor updating with fastQR.", max(lapack_cfg$n), max(lapack_cfg$p), lapack_cfg$batch))


## ----final-ridge-update-table-----------------------------------------------------------
final_ridge_results <- read.csv(file.path(lapack_root, "results/bench_ridge_summary.csv"))
knitr::kable(final_ridge_results, digits = 3, row.names = FALSE,
 caption = "Final ridge comparison for additions of ten observations or five predictors. Speedup is recomputation time divided by updating time; values above one favor updating.")


## ----final-ridge-speedup-summary, include=FALSE-----------------------------------------
final_ridge_largest <- subset(final_ridge_results, n == 300 & p == 450)
final_ridge_rows <- subset(final_ridge_largest, operation == "rows_fixed_alpha")
final_ridge_columns <- subset(final_ridge_largest, operation == "columns_fixed_alpha")


## ----final-ridge-speedup, fig.height=6.4, fig.cap=guide_caption("final-ridge-speedup", "Ridge updating speed relative to refactorization on the same target fit. Values to the right of one favor updating. Panels distinguish predictor additions at fixed total penalty, observation additions at fixed total penalty, and observation additions at fixed normalized penalty; all cases are retained.")----
final_ridge_results$design <- paste0("n = ", final_ridge_results$n, ", p = ", final_ridge_results$p)
final_ridge_results$operation_label <- factor(final_ridge_results$operation,
 levels = c("columns_fixed_alpha", "rows_fixed_alpha", "rows_fixed_lambda"),
 labels = c("Add predictors: fixed total penalty", "Add observations: fixed total penalty",
            "Add observations: fixed normalized penalty"))
ggplot2::ggplot(final_ridge_results, ggplot2::aes(speedup, design)) +
 ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
 ggplot2::geom_point(colour = "#b5121b", size = 2) +
 ggplot2::scale_x_log10() + ggplot2::facet_wrap(~operation_label, ncol = 1) +
 ggplot2::labs(x = "Recomputation time / updating time (log scale)", y = NULL) +
 ggplot2::theme_minimal(base_size = 11)


## ----environment------------------------------------------------------------------------
packageVersion("fastQR")
Sys.info()[c("sysname", "release", "machine")]
Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
sessionInfo()


## ----gram-schmidt-----------------------------------------------------------------------
mgs <- function(X) {
  n <- nrow(X); p <- ncol(X)
  Q <- matrix(0, n, p); R <- matrix(0, p, p)
  for (k in seq_len(p)) {
    v <- X[, k]
    if (k > 1L) for (j in seq_len(k - 1L)) {
      R[j, k] <- sum(Q[, j] * v)
      v <- v - R[j, k] * Q[, j]
    }
    R[k, k] <- sqrt(sum(v^2))
    if (R[k, k] <= .Machine$double.eps * sqrt(sum(X[, k]^2)))
      stop("Numerically dependent column: use a rank-revealing factorization.")
    Q[, k] <- v / R[k, k]
  }
  list(Q = Q, R = R)
}
set.seed(2040)
X_gs <- matrix(rnorm(80 * 8), 80, 8)
gs <- mgs(X_gs)
gs_checks <- c(reconstruction = norm(X_gs - gs$Q %*% gs$R, "F") / norm(X_gs, "F"),
               orthogonality = norm(crossprod(gs$Q) - diag(8), "F"))
gs_checks
stopifnot(max(gs_checks) < 1e-10)

