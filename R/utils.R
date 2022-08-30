#' @export 
read_spectrum <- function(file_input) {
    
    int_ind <- stri_detect(file_input$name, fixed = "intensity")
    
    #intensity file
    path_intensity <- file_input[int_ind, "datapath"]
    intensity <- read.csv(path_intensity, header = TRUE)[,1]
    
    #param file
    path_param <- file_input[!int_ind, "datapath"]
    tmp <- read.csv(path_param, header = TRUE)
    param <- matrix(as.numeric(tmp), nrow = 1)
    colnames(param) <- names(tmp)
  
    info_out <- param
    
    # compute ppm values  
    nspec <- length(intensity)
    info_out <- cbind(info_out, FTSIZE = nspec)
    
    swp <- info_out[1, "SW_p"] / info_out[1, "SF"]
    dppm <- swp / nspec
    ppm <- seq(info_out[1, "OFFSET"], (info_out[1, "OFFSET"] - swp), by = -dppm)
    # the ppm of last point may not coincide with the sequence right limit
    ppm <- ppm[1 : nspec]
    
    names(intensity) <- ppm
    return(list(data = list(intensity), info = info_out))
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
zero_fill_apod2 <- function(x, size, LB, SW_h, GRPDLY = 0) {
  if (!(GRPDLY == 0)){
    n <- length(x)
    x2 <- ph_corr(x, -GRPDLY*2*pi*c((0:(n-1))/n))
    x3 <- zero_fill_apod(x2, size, LB, SW_h)
    x4 <- ph_corr(x3, GRPDLY*2*pi*c((0:(n-1))/n))
    return(x4)
  } else return(zero_fill_apod(x, size, LB, SW_h))
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
obj_fun <- function(x, x_order, ppm, form1, form2, mix, pivot_point, mode = "objective") {
  n_points <- length(mix)
  int_seq <- 1:n_points
  
  # create an index vector that correctly binds optimization parameter names and their positions
  ind <- purrr::map_dbl(c(prop_form2 = "prop_form2", ph0_mix = "ph0_mix", ph1_mix = "ph1_mix", 
                          ppm_form1 = "ppm_form1", ppm_mix = "ppm_mix"), 
                        ~ match(.x, x_order)) 
  x_form1 <- Re(form1)
  if (!is.na(ind["ppm_form1"])) {
    x_form1 <- shift_horizon(x = ppm, y = x_form1, delta = x[ind["ppm_form1"]])
  }                  
  x_form1 <- norm_sum(x_form1)
  
  x_form2 <- norm_sum(form2)
  
  y_mix <- mix
  if (all(!is.na(ind[c("ph0_mix", "ph1_mix")]))){
    lin_pred <- get_ph_angle(x = int_seq, int = x[ind["ph0_mix"]], slope = x[ind["ph1_mix"]], pivot_point = pivot_point, n = n_points, ppm = ppm)
    y_mix <- ph_corr(y_mix, lin_pred)   
  }
  if (!is.na(ind["ppm_mix"])) {
    y_mix <- shift_horizon(x = ppm, y = Re(y_mix), delta = x[ind["ppm_mix"]])
  }
  y_mix <- norm_sum(y_mix)
  y_mix <- y_mix - x_form1
  x_form2 <- x_form2 - x_form1
  fitted <- as.matrix(x_form2) %*% x[ind["prop_form2"]]
  #fitted <- x[ind["prop_form2"]] * x_form2
  residuals <- y_mix - fitted
  # if mode = "objective" then return RSS (used in nloptr); else return input data + fitted values + residuals
  if (mode!="objective"){
    return(data.frame(ppm = ppm, form1 = x_form1, form2 = x_form2 + x_form1, mix = y_mix + x_form1, fitted = fitted + x_form1, residuals = residuals))
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
      "xtol_rel" = 1.0e-6,
      "xtol_abs" = 1.0e-6,
      "maxeval" = 1500,
      print_level = 0
    ),
    x_order = x_order,
    form1 = data$form1,
    mix = data$mix,
    form2 = data$form2,
    ppm = data$ppm, 
    pivot_point = pivot_point,
    mode = "objective"
  )
  
  dat <- obj_fun(x = results$solution, x_order = x_order, ppm = data$ppm, form1 = data$form1, form2 = data$form2, mix = data$mix, pivot_point = pivot_point, mode="prediction")
  dat <- dat[complete.cases(dat), ]
  solution <- c(
    prop_form2 = results$solution[1], rmse = sqrt(results$objective/length(data$mix)), ph0_mix = results$solution[2] * to_deg_const,
    ph1_mix = results$solution[3] * to_deg_const, ppm_form1 = results$solution[4], ppm_mix = results$solution[5], pivot_point = pivot_point
  )
  return(list(dat = dat, solution = solution, start = rlang::set_names(results$x0, x_order)))
}

# separate ppm values for form1, form2, and mix (though residuals and the fit get mix's ppm values)
#' @export 
plot_model_fit <- function(model_fit) {
  p <- plot_ly(x = ~ model_fit$dat$ppm, y = ~ (1 - model_fit$solution["prop_form2"]) * model_fit$dat$form1, type = "scatter", mode = "lines", name = "form1 reference", line = list(width = 2, color = "black"), source = "p1") %>%
    add_trace(x = ~ model_fit$dat$ppm, y = ~ model_fit$solution["prop_form2"] * model_fit$dat$form2, name = "form2 reference", mode = "lines", line = list(color = "green")) %>%
    add_trace(x = ~ model_fit$dat$ppm, y = ~ model_fit$dat$mix, name = "mixture spectrum", mode = "lines", line = list(color = "4682B4")) %>%
    add_trace(x = ~ model_fit$dat$ppm, y = ~ model_fit$dat$residuals, name = "residuals", mode = "lines", line = list(color = "red")) %>%
    add_trace(x = ~ model_fit$dat$ppm, y = ~ model_fit$dat$fitted, name = "fit", mode = "lines", line = list(color = "orange")) %>%
    layout(xaxis = list(zeroline = FALSE, title = "ppm"), yaxis = list(title = "intensity")) 
  return(p)
}

#' @export 
getFilename <- function(path){
	x1 <- readLines(path, n=30)
	fileNameInd <- grep("procs", x1)
	fileName <- sub("\\$\\$ ", "", x1[fileNameInd])
	
	return(fileName)
}
