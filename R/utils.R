`%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))
