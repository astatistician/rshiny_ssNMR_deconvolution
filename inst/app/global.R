library(shiny)
library(plotly)
library(tidyverse)
library(PepsNMR)
library(stringi)

to_deg_const <- 180/pi
to_rad_const <- pi/180

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