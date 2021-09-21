# functions for computations ----------------------------------------------

# reads in spectral data (a vectorized function - multiple spectra at once)
# Arguments:
# type: "spectrum", "fid", or "text".
# path: a character vector, length>=1. If type is set to "spectrum", then path to the folder contating "1r" and "1i" is expected
#       (for Bruker data this is usually something like X\Y\pdata\Z, where X,Y,Z correspond to different experiment, spectra, processed files etc.) )
#       If the type is set to "fid", the path should get you to the folder with the fid file of interest.
#       If the type is set to, the path should lead to the text file
# type

# Value: document the output later on
#' @export 
read_spectrum <- function(path, type) {

  type <- match.arg(type, c("spectrum", "fid", "text"))

  n_files <- length(path)
  data_out <- vector("list", n_files)

  if (type == "fid") {
    params <- c("TD", "BYTORDA", "DIGMOD", "DECIM", "DSPFVS", "SW_h", "SW", "O1", "GRPDLY")
    info_out <- matrix(nrow = n_files, ncol = length(params) + 1)
    colnames(info_out) <- c(params, "DW")

    for (i in 1:n_files) {
      x1 <- readLines(paste0(path[i], "/acqus"))
      x2 <- strsplit(x1, "=")
      x3 <- lapply(x2, function(x) {
        tmp <- x %>%
          stri_replace_all("", regex = "#*\\$*") %>%
          stri_trim_both()
      })

      params_all <- unlist(lapply(x3, function(x) x[1]))
      values_all <- unlist(lapply(x3, function(x) x[2]))
      params_ind <- match(params, params_all)
      values <- as.numeric(values_all[params_ind])
      values <- c(values, 1 / (2 * values[params == "SW_h"]))

      info_out[i, ] <- values
      endianness <- ifelse(info_out[, "BYTORDA"] != 0, "big", "little")

      fid <- readBin(paste0(path[i], "/fid"),
        what = "int", n = info_out[, "TD"],
        size = 4L, endian = endianness
      )

      fid <- complex(
        real = fid[seq(from = 1, to = info_out[i, "TD"], by = 2)],
        imaginary = fid[seq(from = 2, to = info_out[i, "TD"], by = 2)]
      )
      time_vec <- seq(0, (length(fid) - 1) * info_out[i, "DW"], by = info_out[i, "DW"])
      names(fid) <- time_vec
      data_out[[i]] <- fid
    }

    return(list(data = data_out, info = info_out))
  } else if (type == "spectrum") {
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
      spec <- complex(real = spec_r, imaginar = spec_i)
      names(spec) <- ppm
      data_out[[i]] <- spec
    }
    return(list(data = data_out, info = info_out))
  } else if (type == "text") {
    for (i in 1:n_files) {
      tmp_spec <- read.table(path[i], sep = " ", header = FALSE)
      data_out[[i]] <- tmp_spec[, 2]
      names(data_out[[i]]) <- tmp_spec[, 1]
    }
    return(list(data = data_out))
  }
}

# Same as read_spectrum but this version is only needed in case of client-server file interaction
#' @export 
read_spectrum2 <- function(filesDF, type) {
	
	type <- match.arg(type, c("spectrum", "fid", "text"))
	
	n_files <- 1
	data_out <- vector("list", n_files)
	fileNames <- filesDF$name
	filePaths <- filesDF$datapath
	
	if (type == "fid") {
		params <- c("TD", "BYTORDA", "DIGMOD", "DECIM", "DSPFVS", "SW_h", "SW", "O1", "GRPDLY")
		info_out <- matrix(nrow = n_files, ncol = length(params) + 1)
		colnames(info_out) <- c(params, "DW")
		
		for (i in 1:n_files) {
			x1 <- readLines(filePaths[grep("acqus", fileNames)])
			x2 <- strsplit(x1, "=")
			x3 <- lapply(x2, function(x) {
						tmp <- x %>%
								stri_replace_all("", regex = "#*\\$*") %>%
								stri_trim_both()
					})
			
			params_all <- unlist(lapply(x3, function(x) x[1]))
			values_all <- unlist(lapply(x3, function(x) x[2]))
			params_ind <- match(params, params_all)
			values <- as.numeric(values_all[params_ind])
			values <- c(values, 1 / (2 * values[params == "SW_h"]))
			
			info_out[i, ] <- values
			endianness <- ifelse(info_out[, "BYTORDA"] != 0, "big", "little")
			
			fid <- readBin(filePaths[grep("fid", fileNames)],
					what = "int", n = info_out[, "TD"],
					size = 4L, endian = endianness
			)
			
			fid <- complex(
					real = fid[seq(from = 1, to = info_out[i, "TD"], by = 2)],
					imaginary = fid[seq(from = 2, to = info_out[i, "TD"], by = 2)]
			)
			time_vec <- seq(0, (length(fid) - 1) * info_out[i, "DW"], by = info_out[i, "DW"])
			names(fid) <- time_vec
			data_out[[i]] <- fid
		}
		
		return(list(data = data_out, info = info_out))
	} else if (type == "spectrum") {
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
			spec <- complex(real = spec_r, imaginar = spec_i)
			names(spec) <- ppm
			data_out[[i]] <- spec
		}
		return(list(data = data_out, info = info_out))
	} else if (type == "text") {
		for (i in 1:n_files) {
			tmp_spec <- read.table(path[i], sep = " ", header = FALSE)
			data_out[[i]] <- tmp_spec[, 2]
			names(data_out[[i]]) <- tmp_spec[, 1]
		}
		return(list(data = data_out, info = NULL))
	}
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
  approx(x = x, y = y, xout = x + delta)$y
}

#' @export
norm_sum <- function(x) Re(x) / sum(Re(x), na.rm = TRUE)


#' @export 
trim_values <- function(x, limits){
  if (x>=limits[2]) x <- limits[2] else if (x<=limits[1]) x <- limits[1]
  return(x)
}

#' @export 
obj_fun <- function(x, x_order, ppm, amo, cr, mix, pivot_point, mode = "objective") {
  n_points <- length(mix)
  int_seq <- 1:n_points
  
  # create an index vector that correctly binds optimization parameter names and their positions
  ind <- purrr::map_dbl(c(prop_cr = "prop_cr", ph0_mix = "ph0_mix", ph1_mix = "ph1_mix", 
                        ppm_amo = "ppm_amo", ppm_mix = "ppm_mix"), 
                      ~ match(.x, x_order)) 
  x_amo <- Re(amo)
  if (!is.na(ind["ppm_amo"])) {
    x_amo <- shift_horizon(x = ppm, y = x_amo, delta = x[ind["ppm_amo"]])
  }                  
  x_amo <- norm_sum(x_amo)
  
  x_cryst <- norm_sum(cr)
  
  y_mix <- mix
  if (all(!is.na(ind[c("ph0_mix", "ph1_mix")]))){
    lin_pred <- get_ph_angle(x = int_seq, int = x[ind["ph0_mix"]], slope = x[ind["ph1_mix"]], pivot_point = pivot_point, n = n_points, ppm=ppm)
    y_mix <- ph_corr(y_mix, lin_pred)   
  }
  if (!is.na(ind["ppm_mix"])) {
    y_mix <- shift_horizon(x = ppm, y = Re(y_mix), delta = x[ind["ppm_mix"]])
  }
  y_mix <- norm_sum(y_mix)
  fitted <- x[ind["prop_cr"]] * x_cryst + (1 - x[ind["prop_cr"]]) * x_amo
  residuals <- y_mix - fitted
  # if mode = "objective" then return RSS (used in nloptr); else return input data + fitted values + residuals
  if (mode!="objective"){
    return(data.frame(ppm = ppm, amo = x_amo, cr = x_cryst, mix = y_mix, fitted = fitted, residuals = residuals))
  } else return(sum((y_mix - (x[ind["prop_cr"]] * x_cryst + (1 - x[ind["prop_cr"]]) * x_amo))^2, na.rm = TRUE))
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
      "maxeval" = 2000,
      print_level = 3
    ),
    x_order = x_order,
    amo = data$amo,
    mix = data$mix,
    cr = data$cr,
    ppm = data$ppm,
    pivot_point = pivot_point,
    mode = "objective"
  )

  dat <- obj_fun(x = results$solution, x_order = x_order, ppm = data$ppm, amo = data$amo, cr = data$cr, mix = data$mix, pivot_point = pivot_point, mode="prediction")
  dat <- dat[complete.cases(dat), ]
  solution <- c(
    prop_cr = results$solution[1], rmse = sqrt(results$objective/length(data$mix)), ph0_mix = results$solution[2] * to_deg_const,
    ph1_mix = results$solution[3] * to_deg_const, ppm_amo = results$solution[4], ppm_mix = results$solution[5], pivot_point = pivot_point
  )
  return(list(dat = dat, solution = solution, start = rlang::set_names(results$x0, x_order)))
}

# modular optim -----------------------------------------------------------

process_spectrum <- function(spec_cplx, proc_steps, proc_steps_order, x, acq_info){
  norm_ind <- FALSE
  n_points <- length(spec_cplx)
  if (any(str_detect(proc_steps, "norm"))) {
    norm_ind <- TRUE
    proc_steps <- proc_steps %>% str_remove("norm") %>% stri_remove_empty() %>% compact()
  }
  # apodization (exponential function)
  if ("apod" %in% proc_steps) {
    apod_match_ind <- match("apod", proc_steps)
    spec_cplx <- zero_fill_apod(x = spec_cplx, size = n_points, LB = x[proc_steps_order[apod_match_ind]] , SW_h = acq_info[2])
  }
  # phase correction
  if (any(proc_steps %in% c("ph0", "ph1"))) {
    ph0_match_ind <- match("ph0", proc_steps)
    ph1_match_ind <- match("ph1", proc_steps)
    ph0_val <- ifelse(!is.na(ph0_match_ind), x[proc_steps_order[ph0_match_ind]], 0)
    ph1_val <- ifelse(!is.na(ph1_match_ind), x[proc_steps_order[ph1_match_ind]], 0)
    lin_pred <- get_ph_angle(x = 1 : n_points, int = ph0_val, slope = ph1_val, pivot_point = ppm[which.max(Re(spec_cplx))[1]], n = n_points, ppm = ppm)
    spec_cplx <- ph_corr(spec_cplx, lin_pred) 
  }
  # shift
  if ("shift" %in% proc_steps) {
    shift_match_ind <- match("shift", proc_steps)
    spec_cplx <- shift_horizon(x = ppm, y = Re(spec_cplx), delta = x[proc_steps_order[shift_match_ind]])
  }
  # normalization
  if (norm_ind) spec_cplx <- norm_sum(spec_cplx)
  # multiplication factor
  if ("multiply" %in% proc_steps) {
    multiply_match_ind <- match("multiply", proc_steps) 
    spec_cplx <- spec_cplx * x[proc_steps_order[multiply_match_ind]]
  }
  return(spec_cplx)
}

get_param_index <- function(proc_steps, n_prop) {
  proc_steps <- proc_steps %>% map(~str_remove(.x, "norm") %>% stri_remove_empty())
  len <- proc_steps %>% map(length) %>% unlist()
  empty_ind <- len == 0
  proc_steps_flat <- proc_steps %>% unlist()
  x_seq <- (n_prop+1):(n_prop+1+length(proc_steps_flat))
  end_ind <- proc_steps %>% map_dbl(length) %>% cumsum()
  start_ind <- end_ind - (proc_steps %>% map_dbl(length) - 1)
  x_list <- proc_steps
  for (i in seq_along(x_list)) {
    if (!empty_ind[i]) {
      x_list[[i]] <- x_seq[start_ind[i] : end_ind[i]]
    } else {
      x_list[[i]] <- NA
    }
    }
  return(x_list)
}

obj_fun2 <- function(x, proc_steps, y, X, ppm, acq_info, ref_template_label, mode = "objective") {
  X_org <- X
  x_list <- get_param_index(proc_steps, length(X)-1)
  # run proc_spec for each spectrum in X and y
  y <- process_spectrum(spec_cplx = y, proc_steps = pluck(proc_steps, 1), proc_steps_order = pluck(x_list, 1), x = x, acq_info = acq_info)
  tmp_list <- list(X, proc_steps[-1], x_list[-1])
  X <- pmap_dfc(tmp_list, process_spectrum, x = x, acq_info = acq_info)
  X_ref_vector <- X[[ref_template_label]]
  # the linear model
  y <- y - X_ref_vector
  X_full <- X %>% mutate(across(everything(), .fns = ~(.x - X_ref_vector), .names = NULL)) 
  X <- X_full %>% select(-any_of(ref_template_label))
  fitted <- as.matrix(X) %*% x[1:length(X)] %>%  as.numeric()
  residuals <- y - fitted
  if (mode!="objective"){
    return(tibble(y = y + X_ref_vector, X_full + X_ref_vector , ppm = ppm, fitted = fitted + X_ref_vector, residuals = residuals) %>% drop_na())
  } else return(sum(residuals^2, na.rm = TRUE))
}

get_start_constraints <- function(proc_steps, n_prop){
  proc_steps <- proc_steps %>% map(~str_remove(.x, "norm") %>% stri_remove_empty()) %>% compact()
  if (length(proc_steps) >= 1){
    x_list <- get_param_index(proc_steps, n_prop)
  } else x_list <- list()
  param_details <- matrix(NA, nrow = max(c(unlist(x_list), 1)), ncol = 3)
  # prop
  param_details[1 : n_prop, ] <- matrix(rep(c(0, 0, 1), n_prop), byrow = TRUE, ncol = 3)
  if (length(proc_steps) >= 1){
    # other processing params
    default_details <- tribble(~name, ~start, ~lb, ~ub,
                               "apod", 0, 0, 1000,
                               "ph0", 0, 0, 2*pi,
                               "ph1", 0, 0, 2*pi,
                               "shift", 0, -0.2, 0.2,
                               "multiply", 1, 0, 1000)
    param_details[-(1 : n_prop), ] <- map2_dfr(unlist(proc_steps), unlist(x_list), ~{
      ind <- match(.x, default_details$name)
      default_details[ind, -1]
    }) %>% as.matrix()
  }
  colnames(param_details) <- c("start", "lb", "ub")
  return(as_tibble(param_details))
}

nloptr_wrapper2 <- function(dat, obj_fun, proc_steps, start_constraints, optim_algorithm) {
  results <- nloptr::nloptr(
    x0 = start_constraints$start,
    eval_f = obj_fun2,
    lb = start_constraints$lb,
    ub = start_constraints$ub,
    opts = list(
      "algorithm" = optim_algorithm,
      "xtol_rel" = 1.0e-10,
      "xtol_abs" = 1.0e-10,
      "maxeval" = 2000,
      print_level = 3
    ),
    proc_steps = proc_steps,
    y = dat$y,
    X = dat$X,
    ppm = dat$ppm,
    acq_info = dat$acq_info,
    ref_template_label = dat$ref_template_label,
    mode = "objective"
  )
  prop_names <- dat$X %>% select(-dat$ref_template_label) %>% names()
  dat_tmp <- obj_fun2(x = results$solution, proc_steps = proc_steps, y = dat$y, X = dat$X, ppm = dat$ppm, acq_info = dat$acq_info , ref_template_label = dat$ref_template_label, mode = "prediction") %>% drop_na()
  solution <- list(prop = set_names(results$solution[1 : (length(dat$X)-1)], prop_names), rmse = sqrt(results$objective/length(na.omit(dat_tmp$y))))
  return(list(dat = dat_tmp, solution = solution))
}

# functions for graphs ----------------------------------------------------

#' @export 
plot_model_fit <- function(model_fit) {
  p <- plot_ly(x = ~ model_fit$dat$ppm, y = ~ (1 - model_fit$solution["prop_cr"]) * model_fit$dat$amo, type = "scatter", mode = "lines", name = "0% crystal reference", line = list(width = 2, color = "black"), source = "p1") %>%
    add_trace(x = ~ model_fit$dat$ppm, y = ~ model_fit$solution["prop_cr"] * model_fit$dat$cr, name = "100% crystal reference", mode = "lines", line = list(color = "green")) %>%
    add_trace(x = ~ model_fit$dat$ppm, y = ~ model_fit$dat$mix, name = "mixture spectrum", mode = "lines", line = list(color = "4682B4")) %>%
    add_trace(x = ~ model_fit$dat$ppm, y = ~ model_fit$dat$residuals, name = "residuals", mode = "lines", line = list(color = "red")) %>%
    add_trace(x = ~ model_fit$dat$ppm, y = ~ model_fit$dat$fitted, name = "fit", mode = "lines", line = list(color = "orange")) %>%
    layout(xaxis = list(zeroline = FALSE, title = "ppm"), yaxis = list(title = "intensity")) 
  return(p)
}

add_deconvolution_palette <- function(gg_plot, return_palette = FALSE){
  my_pal <- c("#000000", brewer.pal(8, "Dark2"))
  if (return_palette) return(my_pal) else {
    gg_plot <- gg_plot + scale_colour_manual(values = my_pal)
    return(gg_plot)
  }
}

# x - numeric vector with xaxis values; y - numeric vector or matrix/data frame/tibble with intensities
# pass quoted commands via ... to modify ggplot object
# add feature to rescale template spectra according to the estimated proportion values
plot_spectrum <- function(x, y, ... , rev_xaxis = TRUE, interactive = FALSE) {
  gg_dat <- bind_cols(x = x,y) %>% pivot_longer(cols = -x, names_to = "spectrum", values_to = "intensity")
  p <- gg_dat %>% ggplot(aes(x = x, y = intensity, colour = spectrum)) +
    geom_line() +
    theme_bw() +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8))
  if (length(unique(gg_dat)) <=9) p <- add_deconvolution_palette(p)
  add_options <- list(...)  
  for (i in seq_along(add_options)) {
    p <- p + eval(add_options[[i]])
  }
  if (rev_xaxis) p <- p + scale_x_reverse(breaks = scales::pretty_breaks(n = 8)) else {
    p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 8))
  }
  if (!interactive) p else ggplotly(p)
}

# adjust ratio computation to handle >2 templates situations  
plot_calibration_curve <- function(x, y, scale_in, scale_out, ..., intercept = TRUE) {
  if (length(y) != length(x)) stop("x and y vectors are not equally sized")
  if (scale_in=='proportion' & scale_out=='ratio'){
    signal_est <- y/(1-y)
    signal_true <- x/(1-x)
  } else if (scale_in=='ratio' & scale_out=='proportion'){
    signal_est <- 1+1/y
    signal_true <- 1+1/x 
  } else {
    signal_est <- y
    signal_true <- x 
  }
  x_lab <- str_c("true", scale_out, sep = " ")
  y_lab <- str_c("estimated", scale_out, sep = " ")
  gg_dat <- bind_cols(signal_est, signal_true)
  p <- gg_dat %>% ggplot(aes(x = signal_true, y = signal_est)) +
    geom_point() +
    theme_bw() +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(x = x_lab, y = y_lab)
  if (intercept) mod <- lm(signal_est ~ signal_true) else mod <- lm(signal_est ~ -1+signal_true)
  cf <- coef(mod)
  summ <- summary(mod)
  legend_txt1 <- paste(paste0(quote('R^2='), round(summ$r.squared,4)), '\n', 
                    paste0(quote('sigma='), round(summ$sigma,5)), collapse=' ')
  cf_sign <- ifelse(sign(cf)<0, "-", "+")
  cf_abs <- abs(cf)
  legend_txt2 <- paste0('y = ', ifelse(intercept, 
                                     paste0(round(coef(mod)[2],3), 'x', " ", cf_sign[1], " ", round(cf_abs[1],3)),
                                     paste0(round(coef(mod)[1],3), 'x')))
  range_x <- range(signal_true)
  range_y <- range(signal_est)
  p <- p + geom_text(x = range_x[1] + diff(range_x) * 0.15, y = range_y[2] - diff(range_y) * 0.15 , label = legend_txt2) +
    geom_text(x = range_x[2] - diff(range_x) * 0.15, y = range_y[1] + diff(range_y) * 0.15, label = legend_txt1) + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE)
  
  add_options <- list(...)  
  for (i in seq_along(add_options)) {
    p <- p + eval(add_options[[i]])
  }
  p
}

# x - vector with ground truth
# y - vector or matrix with estimates (columns)
# add option to handle ratios instead of proportions
plot_estimate_errors <- function(x, y, error_type, ...) {
  y_org <- y
  y <- as_tibble(y)
  err_fun <- function(x,y, type){
    switch(type,
           "diff" = y-x,
           "relative diff" = ifelse(x != 0, (y-x)/x, NaN),
           "multiplicative err" = ifelse(x != 0, y/x, NaN))
  }
  error_match <- match.arg(error_type, c("diff", "relative diff", "multiplicative err"))
  line_y_intercept <- ifelse(error_match == "multiplicative err", 1, 0)
  gg_dat <- bind_cols(x = x, y) %>%  pivot_longer(cols = -x, names_to = "group", values_to = "y") %>% 
    mutate(err = err_fun(x, y, error_match))
  p <- gg_dat %>% ggplot(aes(x = x, y = err, colour = group)) +
    geom_line() + 
    geom_point(size = 2) +
    #geom_point(aes(fill = group), size=2, color='black', shape=21) +
    geom_hline(yintercept = line_y_intercept, linetype="dashed") +
    theme_bw() +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    labs(x = "true proportion", y = error_match)
  if (length(unique(gg_dat$group)) <=9) p <- add_deconvolution_palette(p)
  if (length(y) == 1) p <- p + theme(legend.position = "none")
  add_options <- list(...)  
  for (i in seq_along(add_options)) {
    p <- p + eval(add_options[[i]])
  }
  p
}

# functions for notebooks -------------------------------------------------

# input: 
  # xlsx with paths and spiked proportions (doe_file); "spectrum"/"fid"/"text" (data_type); 
  # column names to be excluded from labels (ref_template)
# output: 
  # list storing: named list with spectra, where vector name is ppm (data); acquistion parameters from TopSpin (acq_info)
  # infrom from the input xlsx plus labels (doe_info)
prep_data <- function(doe_file, data_type, ref_template) {
  doe <- read_excel(doe_file) 
  # create spectrum labels
  # the column given by ref_template will be excluded from label (assumption: props sum to 1)
  tmp <- doe %>% select(starts_with("p_") & !contains(identity(ref_template)))
  tmp2 <- matrix(names(tmp), nrow = 1) %>%  str_remove("_templ") %>% 
    t() %>% as_tibble() 
  tmp3 <- bind_cols(tmp, tmp2, id = str_c("id", 1 : nrow(tmp) ))
  templ_names <- names(tmp)
  prefix_names <- names(tmp3)
  prefix_names <- prefix_names[str_detect(prefix_names, pattern = "V")]
  glue_var_order <- matrix(c(prefix_names, templ_names), nrow = 2, byrow = TRUE) %>% as.character
  doe <- bind_cols(doe, unite(tmp3, col = "spec_label", all_of(glue_var_order), id))
  # create indicators of template spectra
  a <- doe %>% select(starts_with("p_"))
  b <- a; b[] <- FALSE
  # only the first spectrum with 100% concentration is marked as a template 
  a <- a %>% map(~(min(which(.x==1))))
  templ_ind = map2(a,b, ~{.y[.x] = TRUE; .y}) %>% as_tibble() %>% rename_with(~str_remove(string = .x, pattern = "p_"))
  doe <- bind_cols(doe, templ_ind)
  # get the reference template spectrum label, useful in modelling
  ref_template_label <- doe %>% filter(across(str_remove(ref_template, "p_"))) %>% pull(spec_label)
  # read in the data
  input_dat <- read_spectrum(path = doe$path, type = data_type) 
  input_dat$data <-  input_dat$data %>% purrr::set_names(doe$spec_label)
  return(list(data = input_dat$data, acq_info = input_dat$info, doe_info = doe, n_spec = nrow(doe), ref_template_label = ref_template_label))
}

# switch from list to tibble (more convenient), first column: ppm, other columns: labeled spectra
get_spectral_tibble <- function(dat) {
  #dat %>% map_dfc(~.x) %>% bind_cols(ppm = as.numeric(names(input_dat$data$p1_0_id1))) %>% select(ppm, everything())
  dat %>% map_dfc(~.x)
}

