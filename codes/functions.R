make_quantile_basis <- function(data, ...) {
  sq_list <- list(...)
  
  cum_list <- list()
  int_list <- list()
  cuts_out <- list()
  
  for (sq in sq_list) {
    for (v in names(sq)) {
      x <- data[[v]]
      q <- sort(sq[[v]][, "quantile"])
      p <- as.numeric(rownames(sq[[v]]))[order(sq[[v]][, "quantile"])]
      K <- length(q)
      
      # cumulative: x1_0.25, x1_0.5, x1_0.75
      a <- vapply(q, function(qj) as.numeric(x <= qj), numeric(length(x)))
      colnames(a) <- paste0(v, "_", p)
      
      # interval: x1_q1, x1_q2, ...
      breaks <- c(-Inf, q, Inf)
      at <- matrix(0, nrow = length(x), ncol = K + 1)
      for (j in seq_len(K + 1)) {
        at[, j] <- as.numeric(x > breaks[j] & x <= breaks[j + 1])
      }
      colnames(at) <- paste0(v, "_q", seq_len(K + 1))
      
      cum_list[[v]] <- a
      int_list[[v]] <- at
      cuts_out[[v]] <- setNames(q, p)
    }
  }
  
  list(
    cumulative = do.call(cbind, cum_list),
    interval   = do.call(cbind, int_list),
    cutpoints  = cuts_out
  )
}


make_formulas <- function(basis, y_vars = NULL, x_vars = NULL, type = "cumulative") {
  # basis:  output of make_quantile_basis
  # y_vars: 
  #   - NULL:              returns ~ rhs (one-sided)
  #   - character vector:  returns paste(y_vars, "+") ~ rhs
  #   - named list of one-sided formulas: returns named list of formulas
  #     e.g. list(mod1 = ~ y11 + y12, mod2 = ~ y21 + y22)
  # x_vars: character vector of continuous x's to prepend, or NULL
  # type:   "cumulative" or "interval"
  
  mat <- basis[[type]]
  basis_names <- colnames(mat)
  
  if (!is.null(x_vars)) {
    rhs <- paste(c(x_vars, basis_names), collapse = " + ")
  } else {
    rhs <- paste(basis_names, collapse = " + ")
  }
  
  # no y: one-sided formula
  if (is.null(y_vars)) {
    return(as.formula(paste("~", rhs)))
  }
  
  # named list of formulas: return named list
  if (is.list(y_vars)) {
    out <- lapply(y_vars, function(f) {
      lhs <- deparse(f[[2]])
      as.formula(paste(lhs, "~", rhs))
    })
    return(out)
  }
  
  # character vector: single formula
  lhs <- paste(y_vars, collapse = " + ")
  as.formula(paste(lhs, "~", rhs))
}
