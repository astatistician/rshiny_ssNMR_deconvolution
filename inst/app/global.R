library(shiny)
library(plotly)
library(tidyverse)
library(PepsNMR)
library(stringi)

to_deg_const <- 180/pi
to_rad_const <- pi/180

# this vector determines the order of parameters that will be estimated via nloptr (parameters, constrains, starting values etc.)
estim_order = c("prop_cr", "ph0_mix", "ph1_mix", "ppm_amo", "ppm_mix")

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