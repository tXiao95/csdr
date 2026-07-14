csdr_abort <- function(message, class, stage, variant = NULL, fold = NULL,
                       role = NULL, parent = NULL) {
  condition <- errorCondition(
    message = paste0("[csdr] ", message),
    class = c(class, "csdr_error"),
    stage = stage,
    variant = variant,
    fold = fold,
    role = role,
    parent = parent,
    call = NULL
  )
  stop(condition)
}

csdr_warn <- function(message, class = "csdr_warning", stage,
                      variant = NULL, fold = NULL, role = NULL) {
  condition <- warningCondition(
    message = paste0("[csdr] ", message),
    class = c(class, "csdr_warning"),
    stage = stage,
    variant = variant,
    fold = fold,
    role = role,
    call = NULL
  )
  warning(condition)
  invisible(condition)
}

csdr_rethrow <- function(parent, message, stage, variant = NULL, fold = NULL,
                         role = NULL) {
  if (inherits(parent, "csdr_error")) {
    stop(parent)
  }
  csdr_abort(
    message = paste0(message, " Original error: ", conditionMessage(parent)),
    class = "csdr_fit_error",
    stage = stage,
    variant = variant,
    fold = fold,
    role = role,
    parent = parent
  )
}
