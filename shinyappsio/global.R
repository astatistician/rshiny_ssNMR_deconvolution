library(shiny)
library(plotly)
library(tidyverse)
library(PepsNMR)
library(stringi)

load("sysdata.rda")
source("utils.R")

to_deg_const <- 180/pi
to_rad_const <- pi/180

# this vector determines the order of parameters that will be estimated via nloptr (parameters, constrains, starting values etc.)
#estim_order = c("prop_form2", "ph0_mix", "ph1_mix", "ppm_form1", "ppm_mix")
estim_order = c("prop_form2", "ph0_mix", "ph1_mix", "ppm_mix", "ppm_form1")

param_defaults <- list(
  description_form1 = "",	
  description_form2 = "",	
  description_mix = "",	
  optim_algorithm = "NLOPT_LN_SBPLX",	
  fit_type = "",	
  approve_fit = "",
  comment = "",	
  rmse = NA,
  prop_form2 =	0,
  ph0_mix	= 0,
  ph1_mix	= 0,
  ppm_form1 = 0,	
  ppm_mix	= 0,
  pivot_point = NA,	
  add_zeroes = 0,	
  lb_global = 0,	
  lb_form1	= 0,
  lb_form2	= 0,
  lb_mix = 0,
  ppm_form1_lower	= -0.2,
  ppm_form1_upper	= 0.2,
  ppm_mix_lower	= -0.2,
  ppm_mix_upper	= 0.2,
  ph0_mix_lower	= -180,
  ph0_mix_upper	= 180,
  ph1_mix_lower	= -10 * to_deg_const,
  ph1_mix_upper	= 10 * to_deg_const,
  ppm_range1 = NA,
  ppm_range2 =  NA,
  prop_form2_start = 0,
  ppm_form1_start	= 0,
  ppm_mix_start	= 0,
  ph0_mix_start	= 0,
  ph1_mix_start	= 0
  )

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

param_button_mapping <- tribble(~param, ~label,
                                "ppm_form1", "Chemical shift of form1",
                                "ppm_mix", "Chemical shift of mixture",
                                "ph0_mix", "PH0 of mixture",
                                "ph1_mix", "PH1 of mixture")