# reads in spectral data (a vectorized function - multiple spectra at once)
# path: a character vector, length>=1. If type is set to "spectrum", then path to the folder containing "1r" and "1i" is expected
#       (for Bruker data this is usually something like X\Y\pdata\Z, where X,Y,Z correspond to different experiment, spectra, processed files etc.) )

#' @export 
read_spectrum <- function(path) {
  
    n_files <- length(path)
    data_out <- vector("list", n_files)
  
    params <- c("OFFSET", "SW_p", "SF", "SI", "BYTORDP", "NC_proc", "FTSIZE")
    info_out <- matrix(nrow = n_files, ncol = length(params))
    colnames(info_out) <- c(params)
    
    for (i in 1:n_files) {
      x1 <- readLines(paste0(path[i], "/procs"))
      x2 <- strsplit(x1, "=")
      x3 <- lapply(x2, function(x) {
        tmp <- x %>%
          stri_replace_all("", regex = "#*\\$*") %>%
          stri_trim_both()
      })
      
      params_all <- unlist(lapply(x3, function(x) x[1]))
      values_all <- unlist(lapply(x3, function(x) x[2]))
      params_ind <- match(params, params_all)
      info_out[i, ] <- as.numeric(values_all[params_ind])
      
      nspec <- info_out[i, "FTSIZE"]
      swp <- info_out[i, "SW_p"] / info_out[i, "SF"]
      dppm <- swp / (nspec - 1)
      ppm <- seq(info_out[i, "OFFSET"], (info_out[i, "OFFSET"] - swp), by = -dppm)
      
      spec_r <- readBin(paste0(path[i], "/1r"),
                        what = "int", n = nspec,
                        size = 4L, endian = "little"
      )
      spec_i <- readBin(paste0(path[i], "/1i"),
                        what = "int", n = nspec,
                        size = 4L, endian = "little"
      )
      spec <- complex(real = spec_r, imaginary = spec_i)
      names(spec) <- ppm
      data_out[[i]] <- spec
    }
    return(list(data = data_out, info = info_out))
}

# Same as read_spectrum but this version is only needed in case of client-server file interaction
#' @export 
read_spectrum2 <- function(filesDF) {
  
    n_files <- 1
    data_out <- vector("list", n_files)
    fileNames <- filesDF$name
    filePaths <- filesDF$datapath
  
    params <- c("OFFSET", "SW_p", "SF", "SI", "BYTORDP", "NC_proc", "FTSIZE")
    info_out <- matrix(nrow = n_files, ncol = length(params))
    colnames(info_out) <- c(params)
    
    for (i in 1:n_files) {
      x1 <- readLines(filePaths[grep("procs", fileNames)])
      x2 <- strsplit(x1, "=")
      x3 <- lapply(x2, function(x) {
        tmp <- x %>%
          stri_replace_all("", regex = "#*\\$*") %>%
          stri_trim_both()
      })
      
      params_all <- unlist(lapply(x3, function(x) x[1]))
      values_all <- unlist(lapply(x3, function(x) x[2]))
      params_ind <- match(params, params_all)
      info_out[i, ] <- as.numeric(values_all[params_ind])
      
      nspec <- info_out[i, "FTSIZE"]
      swp <- info_out[i, "SW_p"] / info_out[i, "SF"]
      dppm <- swp / (nspec - 1)
      ppm <- seq(info_out[i, "OFFSET"], (info_out[i, "OFFSET"] - swp), by = -dppm)
      
      spec_r <- readBin(filePaths[grep("1r", fileNames)],
                        what = "int", n = nspec,
                        size = 4L, endian = "little"
      )
      spec_i <- readBin(filePaths[grep("1i", fileNames)],
                        what = "int", n = nspec,
                        size = 4L, endian = "little"
      )
      spec <- complex(real = spec_r, imaginary = spec_i)
      names(spec) <- ppm
      data_out[[i]] <- spec
    }
    return(list(data = data_out, info = info_out))
}

#' @export 
zero_fill_apod <- function(x, size, LB, SW_h) {
  n <- length(x)
  x <- fft(x, inverse = TRUE) / n
  x <- c(x, rep(complex(real = 0, imaginary = 0), size - n))
  tt <- 0:(size - 1)
  x <- x * exp(-tt * LB * pi / (2 * SW_h))
  x <- fft(x)
  return(x)
}

#' @export 
ppm_zero_fill_apod <- function(ppm, info, spec_size){
  swp <- info[2] / info[3]
  dppm <- swp / (spec_size - 1)
  seq(info[1], (info[1] - swp), by = -dppm)
}

#' @export 
ph_corr <- function(x, ph) {
  x * exp(complex(real = 0, imaginary = ph))
}

#' @export 
get_ph_angle <- function(x, int, slope, pivot_point, n, ppm) {
  # the pivot_point value specified by the user may not be exactly present in the data
  # therefore take the closest point
  pivot_point_i <- which.min(abs(pivot_point-ppm))
  slope / n * x + int - slope * pivot_point_i / n
}

#' @export
shift_horizon <- function(x, y, delta) {
  approx(x = x, y = y, xout = x + delta)[[2]]
}

#' @export
norm_sum <- function(x) { 
  rex <- Re(x) 
  rex / sum(rex, na.rm = TRUE)
}

#' @export 
trim_values <- function(x, limits){
  if (x>=limits[2]) x <- limits[2] else if (x<=limits[1]) x <- limits[1]
  return(x)
}

#' @export 
obj_fun <- function(x, x_order, ppm_amo, ppm_cr, ppm_mix, amo, cr, mix, pivot_point, mode = "objective") {
  n_points <- length(mix)
  int_seq <- 1:n_points
  
  # create an index vector that correctly binds optimization parameter names and their positions
  ind <- purrr::map_dbl(c(prop_cr = "prop_cr", ph0_mix = "ph0_mix", ph1_mix = "ph1_mix", 
                          ppm_amo = "ppm_amo", ppm_mix = "ppm_mix"), 
                        ~ match(.x, x_order)) 
  x_amo <- Re(amo)
  if (!is.na(ind["ppm_amo"])) {
    x_amo <- shift_horizon(x = ppm_amo, y = x_amo, delta = x[ind["ppm_amo"]])
  }                  
  x_amo <- norm_sum(x_amo)
  
  x_cryst <- norm_sum(cr)
  
  y_mix <- mix
  if (all(!is.na(ind[c("ph0_mix", "ph1_mix")]))){
    lin_pred <- get_ph_angle(x = int_seq, int = x[ind["ph0_mix"]], slope = x[ind["ph1_mix"]], pivot_point = pivot_point, n = n_points, ppm=ppm_mix)
    y_mix <- ph_corr(y_mix, lin_pred)   
  }
  if (!is.na(ind["ppm_mix"])) {
    y_mix <- shift_horizon(x = ppm_mix, y = Re(y_mix), delta = x[ind["ppm_mix"]])
  }
  y_mix <- norm_sum(y_mix)
  y_mix <- y_mix - x_amo
  x_cryst <- x_cryst - x_amo
  fitted <- as.matrix(x_cryst) %*% x[ind["prop_cr"]]
  #fitted <- x[ind["prop_cr"]] * x_cryst
  residuals <- y_mix - fitted
  # if mode = "objective" then return RSS (used in nloptr); else return input data + fitted values + residuals
  if (mode!="objective"){
    return(data.frame(ppm_amo = ppm_amo, ppm_cr = ppm_cr, ppm_mix = ppm_mix, amo = x_amo, cr = x_cryst + x_amo, mix = y_mix + x_amo, fitted = fitted + x_amo, residuals = residuals))
  } else return(sum(residuals^2, na.rm = TRUE))
}

#' @export 
nloptr_wrapper <- function(data, x_order, obj_fun, param_start, param_constraints, optim_algorithm, pivot_point) {
  results <- nloptr::nloptr(
    x0 = param_start,
    eval_f = obj_fun,
    lb = param_constraints$lb,
    ub = param_constraints$ub,
    opts = nloptr_opts <- list(
      "algorithm" = optim_algorithm,
      "xtol_rel" = 1.0e-10,
      "xtol_abs" = 1.0e-10,
      "maxeval" = 1500,
      print_level = 0
    ),
    x_order = x_order,
    amo = data$amo,
    mix = data$mix,
    cr = data$cr,
    ppm_amo = data$ppm_amo, ppm_cr = data$ppm_cr, ppm_mix = data$ppm_mix,
    pivot_point = pivot_point,
    mode = "objective"
  )
  
  dat <- obj_fun(x = results$solution, x_order = x_order, ppm_amo = data$ppm_amo, ppm_cr = data$ppm_cr, ppm_mix = data$ppm_mix, amo = data$amo, cr = data$cr, mix = data$mix, pivot_point = pivot_point, mode="prediction")
  dat <- dat[complete.cases(dat), ]
  solution <- c(
    prop_cr = results$solution[1], rmse = sqrt(results$objective/length(data$mix)), ph0_mix = results$solution[2] * to_deg_const,
    ph1_mix = results$solution[3] * to_deg_const, ppm_amo = results$solution[4], ppm_mix = results$solution[5], pivot_point = pivot_point
  )
  return(list(dat = dat, solution = solution, start = rlang::set_names(results$x0, x_order)))
}

# ppm values of the amorphous spectrum are used for all traces to be consistent with my analysis notebooks (there, only ppm of ref spectrum used)
#' @export 
plot_model_fit <- function(model_fit) {
  p <- plot_ly(x = ~ model_fit$dat$ppm_amo, y = ~ (1 - model_fit$solution["prop_cr"]) * model_fit$dat$amo, type = "scatter", mode = "lines", name = "0% crystal reference", line = list(width = 2, color = "black"), source = "p1") %>%
    add_trace(x = ~ model_fit$dat$ppm_amo, y = ~ model_fit$solution["prop_cr"] * model_fit$dat$cr, name = "100% crystal reference", mode = "lines", line = list(color = "green")) %>%
    add_trace(x = ~ model_fit$dat$ppm_amo, y = ~ model_fit$dat$mix, name = "mixture spectrum", mode = "lines", line = list(color = "4682B4")) %>%
    add_trace(x = ~ model_fit$dat$ppm_amo, y = ~ model_fit$dat$residuals, name = "residuals", mode = "lines", line = list(color = "red")) %>%
    add_trace(x = ~ model_fit$dat$ppm_amo, y = ~ model_fit$dat$fitted, name = "fit", mode = "lines", line = list(color = "orange")) %>%
    layout(xaxis = list(zeroline = FALSE, title = "ppm"), yaxis = list(title = "intensity")) 
  return(p)
}