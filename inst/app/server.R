## Useful link for I/O server/client: https://ryouready.wordpress.com/2013/11/20/sending-data-from-client-to-server-and-back-using-shiny/

server <- function(input, output, session) {
	
	options(shiny.reactlog = TRUE)
	
	
##### Local variable & function definitions:
	legend_states <- list(TRUE, "legendonly")
	param_defaults <- list(
	path_amo = "",	
	path_cr = "",	
	path_mix = "",	
	optim_algorithm = "NLOPT_LN_SBPLX",	
	fit_type = "",	
	approve_fit = "not approved",
	comment = "NA",	
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
	prop_cr_start =	0,
	ppm_amo_lower	= -1,
	ppm_amo_upper	= 1,
	ppm_mix_lower	= -1,
	ppm_mix_upper	= 1,
	ph0_mix_lower	= 0,
	ph0_mix_upper	= 360,
	ph1_mix_lower	= 0,
	ph1_mix_upper= 360)
	
# create a data frame (empty for the time being) for storing the estimated parameters
	param_defaults_names <- names(param_defaults)
	preproc_param_names <- param_defaults_names[param_defaults_names %in% c("ppm_amo", "ppm_mix", "ph0_mix", "ph1_mix")]
	adv_param_names <- param_defaults_names[param_defaults_names %in% c("add_zeroes", "lb_global", "lb_cr", "prop_cr_start",  "optim_algorithm", "ppm_amo_lower", "ppm_amo_upper", "ppm_mix_lower" ,"ppm_mix_upper", "ph0_mix_lower", "ph0_mix_upper", "ph1_mix_lower", "ph1_mix_upper")] 
	tmp_names <- names(param_defaults)
	saved_params <- matrix(ncol = length(tmp_names) + 1, nrow = 0); colnames(saved_params) <- c("id", tmp_names)
	
	update_n_input_single <- function(id, value) updateNumericInput(session, inputId = id, value = value)
	
	update_n_input_multi <- function(params, dat){
		lapply(as.list(params), function(x) {
					match_ind <- match(x, names(dat))
					update_n_input_single(x, as.numeric(dat[match_ind]))
				})
	}
#####
	
	# load a xlsx file with paths to spectra; these paths will be accessible for choosing 
	# in the path_amo, path_cr, path_mix input fields
	paths <- reactive({
				req(input$file_paths_bn)
				inFile <- input$file_paths_bn
				filePaths_DF <- readxl::read_excel(inFile$datapath)
				colNames <- colnames(filePaths_DF)
				if(is.na(match("path", colNames))){
					showNotification("Provided .xlsx file does not have the right format", type = "error")
					filePaths_DF <- NULL
				}
				return(filePaths_DF)
			}, label = "load xlsx file with paths to spectra")
	
	# update the individual path fields from the xlsx above
	observe({
				updateSelectInput(session, "path_amo", choices = paths()[,1])
				updateSelectInput(session, "path_cr", choices = paths()[,1])
				updateSelectInput(session, "path_mix", choices = paths()[,1])
			}, label = "update path fields from the xlsx")
	
	# load parameter constraints for nloptr
	param_constraints <- reactive({
	      tmp_input <- reactiveValuesToList(input)
	      tmp_ind <- stri_detect(adv_param_names, fixed = "lower")
	      lowerVals <- unlist( tmp_input[adv_param_names[tmp_ind]] ) * ifelse(stri_detect(adv_param_names[tmp_ind], fixed = "ph"), to_rad_const, 1)
	      tmp_ind <- stri_detect(adv_param_names, fixed = "upper")
	      upperVals <- unlist( tmp_input[adv_param_names[tmp_ind]] ) * ifelse(stri_detect(adv_param_names[tmp_ind], fixed = "ph"), to_rad_const, 1)
				varNames <- stri_replace_all_fixed(names(lowerVals), "_lower", "")
				dat <- tibble::tibble(name = varNames, lb = lowerVals, ub = upperVals)
				dat <- rbind(c(NA, 0,1), dat); dat[1,1] <- "prop_cr"
				# make sure that constraints are in correct order
				dat <- dat[match(dat$name, estim_order), ]
				return(dat)      
			}, label = "load the user's parameter constraints for nloptr")
	
	# read in source data
	raw_data <- reactive({
				amo <- tryCatch(my_readin(input$path_amo, "spectrum"), error = function(err) return(err))
				if(class(amo)[1] == "simpleError"){
					showNotification(paste0("Error loading amorphous spectrum files: ", amo$message), type = "error")
					return(NULL)
				}
				info <- amo[[2]]
				spec_size <- info[length(info)] + input$add_zeroes
				amo <- amo$data[[1]]
				ppm <- as.numeric(names(amo))
				cr <- tryCatch(my_readin(input$path_cr, "spectrum"), error = function(err) return(err))
				if(class(cr)[1] == "simpleError"){
					showNotification(paste0("Error loading crystalline spectrum files: ", cr$message), type = "error")
					return(NULL)
				}
				cr <- cr$data[[1]]
				mix <- tryCatch(my_readin(input$path_mix, "spectrum"), error = function(err) return(err))
				if(class(mix)[1] == "simpleError"){
					showNotification(paste0("Error loading mixture spectrum files: ", mix$message), type = "error")
					return(NULL)
				}
				mix <- mix$data[[1]]
				return(list(data.frame(ppm = ppm, amo = amo, cr = cr, mix = mix), info, spec_size))
			}, label = "read in spectra")
	
	# update the pivot point field with the ppm value corresponding to 
	# the most abundant peak of the selected mixture spectrum
	observeEvent(input$path_mix, {
				tmp <- raw_data()
				req(tmp, cancelOutput = TRUE)
				updateNumericInput(session, inputId = "pivot_point", value = 
								tmp[[1]]$ppm[which.max(abs(Re(tmp[[1]]$mix)))])
			},ignoreInit = TRUE, label = "update pivot point to the largest peak of mix")
	
	# process all three spectra at once (i.e. zero filling and apodization)
	proc_data <- reactive({
				raw_data <- raw_data()
				
				req(raw_data, cancelOutput = TRUE)
				ppm <- raw_data[[1]]$ppm
				amo <- raw_data[[1]]$amo
				cr <- raw_data[[1]]$cr
				mix <- raw_data[[1]]$mix
				info <- raw_data[[2]]
				spec_size <- raw_data[[3]]
				
				if (input$add_zeroes > 0 || input$lb_global > 0) {
					swp <- info[2] / info[3]
					dppm <- swp / (spec_size - 1)
					ppm <- info[1]
					ppm <- seq(info[1], (info[1] - swp), by = -dppm)
					amo <- zero_fill_apod(amo, spec_size, input$lb_global, info[2])
					cr <- zero_fill_apod(cr, spec_size, input$lb_global, info[2])
					mix <- zero_fill_apod(mix, spec_size, input$lb_global, info[2])
				} 
				
				# additional line broadening of the crystalline template
				if (input$lb_cr > 0) {
					cr <- apod_only(cr, input$lb_cr, info[2])
				}
				
				return(list(data.frame(ppm = ppm, amo = amo, cr = cr, mix = mix), info, spec_size))
			}, label = "process three spectra (zero fill, apodization)")
	
	# create reactive values:
	# A and M, to control model_fit object - whether the manual or automated fitting procedure will be applied
	fit_bn_rv <- reactiveValues(A = 0, A.counter=0, M = 0, M.counter=0)
	observeEvent(input$manual_fit_bn, {
				fit_bn_rv$M <- 1
				fit_bn_rv$A <- 0
				fit_bn_rv$M.counter <- fit_bn_rv$M.counter+1
			}, label = "update the indicator and counter for manual fitting")
	observeEvent(input$automated_fit_bn, {
				fit_bn_rv$M <- 0
				fit_bn_rv$A <- 1
				fit_bn_rv$A.counter <- fit_bn_rv$A.counter+1
			}, label = "update the indicator and counter for automated fitting")
	
	# manual fitting: the crystalline proportion is estimated via NNLS based on the pre-processed data (according to the user input paramaters),
	model_fit_M <- eventReactive(input$manual_fit_bn, {
				validate(need(input$path_amo, "Amorphous spectrum path required"), 
						need(input$path_cr, "Crystalline spectrum path required"), 
						need(input$path_mix, "Mixture spectrum path required"))
				
				notify_id <- showNotification("Manual model fitting is running...", duration = NULL, closeButton = FALSE, type = "message")
				on.exit(removeNotification(notify_id), add = TRUE)
				
				proc_data <- proc_data()
				
				ppm <- proc_data[[1]]$ppm
				# process the amorphous spectrum
				amo <- proc_data[[1]]$amo
				amo <- ppm_shift(x = ppm, y = Re(amo), delta = input$ppm_amo)
				amo <- norm(amo)
				# process the crystalline spectrum
				cr <- proc_data[[1]]$cr
				cr <- norm(cr)
				# process mixture spectrum
				mix <- proc_data[[1]]$mix
				N <- length(ppm)
				lin_pred <- ph1_pred(1:N, int = input$ph0_mix * to_rad_const, slope = input$ph1_mix * to_rad_const, pivot_point = input$pivot_point, n = N, ppm)
				mix <- ph_corr(mix, lin_pred)
				mix <- ppm_shift(x = ppm, y = Re(mix), delta = input$ppm_mix)
				mix <- norm(mix)
				
				model_input <- data.frame(ppm = ppm, amo = amo, cr = cr, mix = mix)
				na_ind <- is.na(model_input$amo) | is.na(model_input$mix)
				model_input <- model_input[!na_ind, ]
				cr_mod <- model_input$cr - model_input$amo
				mix_mod <- model_input$mix - model_input$amo
				mat_cryst <- cbind(cr_mod)
				nnls_fit <- nnls::nnls(mat_cryst, mix_mod)
				dat <- cbind(model_input, fitted = nnls_fit$fitted + model_input$amo, residuals = nnls_fit$residuals)
				solution <- c(prop_cr = nnls_fit$x, rmse = sqrt(nnls_fit$deviance/nrow(nnls_fit$residuals)))
				return(list(dat = dat, solution = solution))
			}, label = "do manual fitting")
	
	# automated fitting: the crystalline proportion, phase correction, chemical shift are jointly optimized via NLOPTR
	model_fit_A <- eventReactive(input$automated_fit_bn, {
				validate(need(input$path_amo, "Amorphous spectrum path required"), 
						need(input$path_cr, "Crystalline spectrum path required"), 
						need(input$path_mix, "Mixture spectrum path required"))
				
				notify_id <- showNotification("Automated model fitting is running...", duration = NULL, closeButton = FALSE, type = "message")
				on.exit(removeNotification(notify_id), add = TRUE)
				
				proc_data <- proc_data()
				
				ppm <- proc_data[[1]]$ppm
				amo <- proc_data[[1]]$amo
				cr <- proc_data[[1]]$cr
				mix <- proc_data[[1]]$mix
				
				ph0_angle <- input$ph0_mix
				# if ph0_angle==0 then take pepsNMR PH0 as a starting value for PH0
				if (ph0_angle==0) {
					tmp <- matrix(mix, nrow = 1)
					colnames(tmp) <- ppm
					ph0_angle <- ZeroOrderPhaseCorrection(tmp, type.zopc = "max", returnAngle = TRUE)$Angle * to_deg_const
				}
				
				model_input <- data.frame(ppm = ppm, amo = amo, cr = cr, mix = mix)
				na_ind <- is.na(model_input$amo) | is.na(model_input$mix)
				model_input <- model_input[!na_ind, ]
				
				param_constraints_tmp <- param_constraints()
				
				# if user input values exceeds constraints specified for nloptr, then set starting values to the constraint values
				inputs_list <- reactiveValuesToList(input)
				param_start_tmp <- c(prop_cr = input$prop_cr_start, inputs_list[estim_order]) %>% 
				  keep(~!is.null(.x)) %>% 
				  data.frame() %>% 
				  pivot_longer(cols = everything(), names_to = "name", values_to = "start") %>% 
				  inner_join(param_constraints_tmp, by= "name") %>%
				  apply(1, function(x) {
				    x <- as.numeric(x[-1])
				    trim_values(x[1], x[2:3])
				    }
				  )
				
				mod <- nloptr_wrapper(
				      model_input
				    , x_order = estim_order
						, obj_fun = obj_fun
						, param_start = param_start_tmp
						, param_constraints = param_constraints_tmp
						, optim_algorithm = input$optim_algorithm
						, pivot_point = input$pivot_point)
				return(mod)
			}, label = "do automated fitting")
	
	# assign accordingly either the manual or automated fitting routine results
	model_fit <- eventReactive(c(input$manual_fit_bn, input$automated_fit_bn), {
				if (fit_bn_rv$M == 1 & fit_bn_rv$A == 0) {
					return(model_fit_M())
				} else
				if (fit_bn_rv$M == 0 & fit_bn_rv$A == 1) {
					return(model_fit_A())
				}
			}, label = "assign manual or automated fitting results to an object")
	
	## operations on the user inputs
	# after the automated fitting procedure, update the corresponding pre-processing inputs
	observeEvent(input$automated_fit_bn, {
				update_n_input_multi(preproc_param_names, model_fit()$solution)
			}, label = "update user inputs with nloptr estimates")
	
	# reset all input parameters to their default values
	observeEvent(input$reset_bn, {
				#lapply(as.list(proc_param_names[proc_param_names != "pivot_point"]), function(x) update_n_input_single(x, 0))
	       tmp_name <- c(preproc_param_names, adv_param_names[adv_param_names != "optim_algorithm"])
	       update_n_input_multi(tmp_name, unlist(param_defaults[tmp_name]))
	       updateTextInput(session, inputId = "optim_algorithm", value = param_defaults[["optim_algorithm"]])
			}, label = "reset all user inputs to the default values")
	
	# write the estimated parameters to a csv file
	track_inputs_rv <- reactiveValues(x = saved_params, count_rows=0)
	
	observeEvent(c(input$path_amo, input$path_cr, input$path_mix), {
				track_inputs_rv$count_rows <- 0
			}, label = "reset row counts to zero after changing a spectrum")
	
	# each time manual or automated fitting buton is pressed, update the underlying estimate tracking file
	observeEvent(c(input$manual_fit_bn, input$automated_fit_bn), {
				track_inputs_rv$count_rows <- track_inputs_rv$count_rows+1
				tmp_list <- reactiveValuesToList(input)
				next_row <- character(length(param_defaults) + 1); names(next_row) <- c("id", param_defaults_names)
				
				next_row[param_defaults_names] <- unlist(param_defaults)
				vals <- tmp_list[param_defaults_names] %>% purrr::keep(~ !is.null(.x)) %>%  unlist()
				next_row[names(vals)] <- vals
				ms <- model_fit()$solution
				next_row[c("id", "rmse", "prop_cr")] <- c(track_inputs_rv$count_rows, ms["rmse"], ms["prop_cr"])
				
				if (fit_bn_rv$M == 1 & fit_bn_rv$A == 0) {
				  next_row[c("fit_type", "optim_algorithm")] <- c("manual", "NA")
				} else if (fit_bn_rv$M == 0 & fit_bn_rv$A == 1) {
					ms <- model_fit()$solution
					next_row[c("fit_type", preproc_param_names)] <- c("auto", ms[preproc_param_names])
				}
				track_inputs_rv$x <- rbind(track_inputs_rv$x, next_row)
			}, ignoreInit = TRUE, label = "update the underlying estimate tracking file")
	
	# restore parameter values from the previous model fitting step
	observeEvent(input$undo_bn, {
				# this condition ensures that at least two records are in the results file
				if ((fit_bn_rv$M.counter+fit_bn_rv$A.counter)>=2) {
					update_n_input_multi(c(preproc_param_names, adv_param_names[adv_param_names != "optim_algorithm"]),
					                     track_inputs_rv$x[nrow(track_inputs_rv$x)-1,])
				} else{
					showNotification("No previous model fit available", type = "warning")
				}
			}, label = "update user inputs with the previous results")
	
	observeEvent(input$approve_bn, {
				track_inputs_rv$x[nrow(track_inputs_rv$x), "approve_fit"] <- "approved"
				showNotification("The fit has been approved", type = "message", duration = 2)
			})
	
	observeEvent(input$save_user_comments_bn, {
				track_inputs_rv$x[nrow(track_inputs_rv$x), "comment"] <- input$user_comments
				showNotification("The comment has been saved", type = "message", duration = 2)
			})
	##
	
	# event for retaining zoom between data recalulcations in plotly graph
	zoom <- reactive(event_data("plotly_relayout", source = "p1"))
	
	# event for retaining only selected lines between data recalulcations in plotly graph
	legend_click <- reactive(event_data("plotly_legendclick", source = "p1"))
	
	# register clicking event to update Chemical shift (ppm) input value
	plot_click <- reactive(event_data("plotly_click", source = "p1"))
	
	# showing/hiding lines by clicking on the legend 
	legend_items <- reactiveValues("0% crystal reference" = TRUE, "100% crystal reference" = TRUE, "mixture spectrum" = TRUE, "residuals" = TRUE, "fit" = TRUE)
	observe({
				if (!is.null(legend_click()$name)) {
					legend_items[[legend_click()$name]] <- ifelse(legend_click()$visible == TRUE, legend_states[[2]], legend_states[[1]])
				}
			})
	
	observeEvent(plot_click(), {
				updateNumericInput(session, inputId = "pivot_point", value = plot_click()$x)
			})
	
	
	output$fit_plot_out <- renderPlotly({
				req(model_fit())
				p <- plot_model_fit(model_fit())
				# the output of plotly_relayout changes depending on the action: if you zoom in
				# xaxis.range[0], xaxis.range[1], yaxis.range[0],	yaxis.range[1] values will be available
				# if you don't zoom in, the output will be NULL or other some other named vectors of length different than 4
				if (is.null(zoom()) | length(zoom()) != 4) {
					p <- layout(p, xaxis = list(zeroline = FALSE, autorange = "reversed", title = "ppm"), yaxis = list(title = "intensity"))
				} else {
					p <- layout(p, xaxis = list(zeroline = FALSE, range = c(zoom()$"xaxis.range[0]", zoom()$"xaxis.range[1]"), title = "ppm"), yaxis = list(title = "intensity", range = c(zoom()$"yaxis.range[0]", zoom()$"yaxis.range[1]")))
				}
				p <- plotly_build(p)
				for (i in seq_along(p$x$data)) {
					p$x$data[[i]]$visible <- legend_items[[p$x$data[[i]]$name]]
				}
				
				return(p)
			})
	
	
	output$fit_stats_out <- renderText({
				req(model_fit())
				paste0(
						"Estimated crystalline proportion (%): ", round(100 * model_fit()$solution["prop_cr"], 7),
						"\nrmse (x10^6):", round(10^6 * model_fit()$solution["rmse"], 7),
						"\nInitial spectrum size:", isolate(raw_data()[[2]][length(raw_data()[[2]])]),
						"\nFinal spectrum size:", isolate(raw_data()[[2]][length(raw_data()[[2]])] + input$add_zeroes)
				)	
			})
	
	# download a csv file with the estimated parameters
	output$download_params_bn <- downloadHandler(
			filename = function(){paste0("estimated_parameters_", strftime(Sys.time(), "%Y%m%d"), ".csv")},
			content = function(file){
				write.csv(x=track_inputs_rv$x, file=file, row.names=FALSE)
			})
	
	# load the spectra and estimated processing parameters from a file with previous results 
	# only the last record from this file will be read in
	observeEvent(input$file_load_bn, {
				req(input$file_load_bn)
				inFile <- input$file_load_bn
				tmp <- read.csv(inFile$datapath, stringsAsFactors = FALSE)
				
				# Since only the last row will be loaded in, the user has to manually remove other records.
				# The removal can be done in several ways: Delete, Backspace button on entire rows or selected cells.
				# As a consequence, sometimes the removed rows may still be loaded in, but they will contain NA values (of different types).
				# Heuristic: in such case remove rows with "empty" id, path_amo, path_cr, path_mix columns.
				var_check <- c("id", "path_amo", "path_cr", "path_mix")
				complete_cases_ind <- !apply(tmp, 1, function(x) all(is.na(x[var_check]) | x[var_check] == "" ))
				tmp <- tmp[complete_cases_ind,]
				dat <- tmp[nrow(tmp),]
				updateSelectInput(session, "path_amo", choices = dat[["path_amo"]])
				updateSelectInput(session, "path_cr", choices = dat[["path_cr"]])
				updateSelectInput(session, "path_mix", choices = dat[["path_mix"]])
				#update_n_input_multi(c("ppm_amo", "ppm_mix", "ph0_mix", "pivot_point", "ph1_mix", "add_zeroes", "lb_global", "lb_cr"), dat)
				param_selected <- c(preproc_param_names, adv_param_names[adv_param_names != "optim_algorithm"])
				update_n_input_multi(param_selected, dat[param_selected])
				updateTextInput(session, inputId = "optim_algorithm", value = dat["optim_algorithm"])
			}, label = "load in estimated results")
}



