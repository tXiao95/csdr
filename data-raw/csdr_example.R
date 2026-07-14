set.seed(2026)

n <- 60L
C <- matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L)
colnames(C) <- paste0("C", seq_len(ncol(C)))

A <- cbind(
  A1 = 0.6 * C[, 1L] - 0.2 * C[, 2L] + stats::rnorm(n),
  A2 = 0.4 * C[, 2L] + 0.3 * C[, 3L] + stats::rnorm(n)
)
Y <- A[, 1L] + 0.5 * A[, 2L] + C[, 1L] + stats::rnorm(n)

observation_names <- paste0("obs", seq_len(n))
names(Y) <- observation_names
rownames(A) <- observation_names
rownames(C) <- observation_names

csdr_example <- list(Y = Y, A = A, C = C)
save(csdr_example, file = "data/csdr_example.rda", compress = "xz")
