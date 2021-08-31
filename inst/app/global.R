library(shiny)
library(plotly)
library(tidyverse)
library(PepsNMR)
library(stringi)

to_deg_const <- 180/pi
to_rad_const <- pi/180

# this vector determines the order of parameters that will be estimated via nloptr (parameters, constrains, starting values etc.)
estim_order = c("prop_cr", "ph0_mix", "ph1_mix", "ppm_amo", "ppm_mix")

param_defaults <- list(
  path_amo = "",	
  path_cr = "",	
  path_mix = "",	
  optim_algorithm = "NLOPT_LN_SBPLX",	
  fit_type = "",	
  approve_fit = "",
  comment = "",	
  rmse = NA,
  prop_cr =	0,
  ph0_mix	= 0,
  ph1_mix	= 0,
  ppm_amo = 0,	
  ppm_mix	= 0,
  pivot_point = NA,	
  add_zeroes = 0,	
  lb_global = 0,	
  lb_cr	= 0,
  ppm_amo_lower	= -0.2,
  ppm_amo_upper	= 0.2,
  ppm_mix_lower	= -0.2,
  ppm_mix_upper	= 0.2,
  ph0_mix_lower	= 0,
  ph0_mix_upper	= 359.9999,
  ph1_mix_lower	= 0,
  ph1_mix_upper	= 359.9999,
  prop_cr_start = 0,
  ppm_amo_start	= 0,
  ppm_mix_start	= 0,
  ph0_mix_start	= 0,
  ph1_mix_start	= 0)

optim_algorithms_list <-
		c("NLOPT_GN_DIRECT"
				,"NLOPT_GN_DIRECT_L"
				,"NLOPT_GN_DIRECT_L_RAND"
				,"NLOPT_GN_DIRECT_NOSCAL"
				,"NLOPT_GN_DIRECT_L_NOSCAL"
				,"NLOPT_GN_DIRECT_L_RAND_NOSCAL"
				,"NLOPT_GN_ORIG_DIRECT"
				,"NLOPT_GN_ORIG_DIRECT_L"
				,"NLOPT_LN_PRAXIS"
				,"NLOPT_GN_CRS2_LM"
				,"NLOPT_LN_COBYLA"
				,"NLOPT_LN_NEWUOA"
				,"NLOPT_LN_NEWUOA_BOUND"
				,"NLOPT_LN_NELDERMEAD"
				,"NLOPT_LN_SBPLX"
				,"NLOPT_LN_BOBYQA"
				,"NLOPT_GN_ISRES")