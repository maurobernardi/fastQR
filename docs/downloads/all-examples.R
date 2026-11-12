## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = TRUE, message = FALSE,
                     error = FALSE, fig.width = 8, fig.height = 4.8,
                     fig.align = "center", out.width = "100%", dpi = 160, cache = FALSE)
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
  Gram_identity = relative_error(crossprod(thin$R), crossprod(X)))
knitr::kable(data.frame(check = names(checks), relative_error = unname(checks)),
             digits = 16)
stopifnot(max(checks) < 1e-10)


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


## ----sliding-path, fig.cap="Coefficient estimates as the window advances. Solid lines use updated factors; dashed lines are fresh least-squares fits. The two computations agree to numerical precision."----
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
stopifnot(max(abs(ols$coeff - lm.fit(design, y)$coefficients)) < 1e-10)
knitr::kable(data.frame(
  term = c("Intercept", paste0("x", seq_len(p))),
  truth = c(1, beta_true), estimate = as.numeric(ols$coeff)), digits = 4)
ols_test_mse <- mean((y_test - drop(design_test %*% ols$coeff))^2)
ols_test_mse


## ----regression-predictions, fig.cap="Independent test responses versus least-squares predictions. The diagonal is a reference line, not a fitted trend."----
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
path_slopes <- sweep(path_standardized, 1, x_scale, "/")
path_intercept <- y_mean - drop(crossprod(x_mean, path_slopes))

# Independent augmented-design QR calculation at one grid value.
j <- 20L
reference <- qr.coef(base::qr(rbind(Z, sqrt(n * lambda[j]) * diag(p))),
                     c(y_centered, rep(0, p)))
stopifnot(max(abs(path_standardized[, j] - reference)) < 1e-10)


## ----ridge-coefficients, fig.cap="Ridge slopes on their original predictor scale, plotted against a logarithmic penalty grid. Coefficients shrink continuously; ridge does not generally create an exactly sparse fit."----
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


## ----ridge-cv-plot, fig.cap="Validation mean squared error along the penalty grid. The dashed line marks the minimum on this grid."----
ggplot2::ggplot(data.frame(lambda, MSE = cv_mse), ggplot2::aes(lambda, MSE)) +
  ggplot2::geom_line(colour = "#b5121b", linewidth = .9) +
  ggplot2::geom_vline(xintercept = lambda_selected, linetype = "dashed", colour = "#526070") +
  ggplot2::scale_x_log10() + ggplot2::labs(x = expression(lambda), y = "Cross-validated MSE") +
  ggplot2::theme_minimal(base_size = 12)


## ----ridge-df, fig.cap="Effective degrees of freedom for each fixed ridge penalty, including the unpenalized intercept."----
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
stopifnot(max(abs(multivariate$coeff - reference)) < 1e-10)
knitr::kable(as.matrix(multivariate$coeff), digits = 4,
             col.names = c("Response 1", "Response 2"))


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
  b_qr <- backsolve(f$R, crossprod(f$Q, y))
  b_normal <- tryCatch(solve(crossprod(X), crossprod(X, y)),
                        error = function(e) rep(NA_real_, p))
  data.frame(condition = kappa_target,
    reconstruction = norm(X - f$Q %*% f$R, "F") / norm(X, "F"),
    QR_coefficients = sqrt(sum((b_qr - beta)^2) / sum(beta^2)),
    normal_equations = sqrt(sum((b_normal - beta)^2) / sum(beta^2)))
}))
knitr::kable(accuracy, digits = 12)
stopifnot(max(accuracy$reconstruction) < 1e-10)


## ----conditioning-plot, fig.cap="Relative reconstruction and coefficient errors versus the prescribed condition number. Missing normal-equation values indicate that the reference solver rejected the system as numerically singular."----
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


## ----timing-plot, fig.cap="Measured time to obtain factors after adding ten observations. Initial factors are assumed available for the update. Points and boxes summarize five repetitions; each repetition averages a calibrated batch of calls."----
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
    ggplot2::labs(x = "Number of columns", y = "Milliseconds per call (log scale)", colour = NULL,
                  subtitle = paste("n =", benchmark_n, "; repetitions =", benchmark_reps)) +
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


## ----restored-qr-plot, fig.cap="QR interface timings across column dimensions. Compact and explicit outputs are distinguished in the legend."----
benchmark_plot("QR")


## ----restored-cholesky-plot, fig.cap="Triangular-factor timings, with separate baselines for an available Gram matrix and for its formation from X."----
benchmark_plot("Cholesky")


## ----restored-solve-plot, fig.cap="Solution times for the same precomputed normal-equation system."----
benchmark_plot("Solve")


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

