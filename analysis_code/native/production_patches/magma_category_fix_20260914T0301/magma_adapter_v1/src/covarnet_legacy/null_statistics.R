# Explicit missing-statistic handling and exact fixed-margin conservation.
empirical_upper <- function(observed, null) {
  if (length(observed) != 1L || !is.finite(observed)) return(NA_real_)
  finite <- null[is.finite(null)]
  if (!length(finite)) return(NA_real_)
  (1 + sum(finite >= observed)) / (1 + length(finite))
}

validate_count_matrix <- function(x) {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0) ||
      any(x != floor(x)) || any(x > .Machine$integer.max)) {
    stop("Counts must be exact finite nonnegative integers before conversion")
  }
  storage.mode(x) <- "integer"
  x
}

checked_parent_null <- function(counts, parent) {
  stopifnot(length(parent) == nrow(counts), !anyNA(parent))
  result <- matrix(0L, nrow(counts), ncol(counts), dimnames = dimnames(counts))
  for (lineage in unique(parent)) {
    idx <- which(parent == lineage)
    block <- counts[idx, , drop = FALSE]
    if (length(idx) == 1L || sum(block) == 0L) {
      generated <- block
    } else {
      generated <- r2dtable(1L, rowSums(block), colSums(block))[[1L]]
    }
    stopifnot(all(rowSums(generated) == rowSums(block)),
              all(colSums(generated) == colSums(block)))
    result[idx, ] <- generated
  }
  result
}
