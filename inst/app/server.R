## Useful link for I/O server/client: https://ryouready.wordpress.com/2013/11/20/sending-data-from-client-to-server-and-back-using-shiny/

# mapping: amo - form1, cr - form2

server <- function(input, output, session) {
	
	#options(shiny.reactlog = TRUE)
	
##### Local variable & function definitions:
	legend_states <- list(TRUE, "legendonly")
	
# create a data frame (empty for the time being) for storing the estimated parameters
	param_defaults_names <- names(param_defaults)
	preproc_param_names <- param_defaults_names[param_defaults_names %in% c("prop_form2", "ppm_form1", "ppm_mix", "ph0_mix", "ph1_mix")]
	adv_param_names <- param_defaults_names[param_defaults_names %in% 
	                                          c("add_zeroes", "lb_global", "lb_form2", "optim_algorithm", "ppm_form1_lower", "ppm_form1_upper", "ppm_mix_lower" ,"ppm_mix_upper", "ph0_mix_lower", "ph0_mix_upper", "ph1_mix_lower", "ph1_mix_upper")]
	start_val_names <- param_defaults_names[str_detect(param_defaults_names, "_start")]
	saved_params <- matrix(ncol = length(param_defaults_names) + 1, nrow = 0); colnames(saved_params) <- c("id", 	param_defaults_names)
	
	# List of reactive variables
	results <- reactiveValues(
			form1 = NULL,
			form2 = NULL,
			mix = NULL,
			resultsFileLoaded = FALSE
	)
	
	update_n_inputs <- function(params, dat){
		lapply(as.list(params), function(x) {
					match_ind <- match(x, names(dat))
					updateNumericInput(session, inputId = x, value = as.numeric(dat[match_ind]))
				})
	}
	
#####
	
	observeEvent(c(input$path_form1, input$files_form1), {
				
				if(localMode || input$inputSource == "Shared drive"){
					form1 <- tryCatch(read_spectrum(input$path_form1), error = function(err) return(err))
				} else {
					form1 <- tryCatch(read_spectrum2(input$files_form1), error = function(err) return(err))
				} 
				
				if(class(form1)[1] == "simpleError"){
					showNotification(paste0("Error loading form1 spectrum files: ", form1$message), type = "error")
				} else{ 
					results$form1 <- form1
					showNotification("form1 spectrum files successfully loaded", type = "message")
				}
			}, ignoreInit = TRUE)
	
	observeEvent(c(input$path_form2, input$files_form2), {
				
				if(localMode || input$inputSource == "Shared drive"){
					form2 <- tryCatch(read_spectrum(input$path_form2), error = function(err) return(err))
				} else {
					form2 <- tryCatch(read_spectrum2(input$files_form2), error = function(err) return(err))
				} 
				
				if(class(form2)[1] == "simpleError"){
					showNotification(paste0("Error loading form2 spectrum files: ", form2$message), type = "error")
				} else{ 
					results$form2 <- form2
					showNotification("form2 spectrum files successfully loaded", type = "message")
				}
			}, ignoreInit = TRUE)
	
	observeEvent(c(input$path_mix, input$files_mix), {
				
				if(localMode || input$inputSource == "Shared drive"){
					mix <- tryCatch(read_spectrum(input$path_mix), error = function(err) return(err))
				} else {
					mix <- tryCatch(read_spectrum2(input$files_mix), error = function(err) return(err))
				} 
				
				if(class(mix)[1] == "simpleError"){
					showNotification(paste0("Error loading mixture spectrum files: ", mix$message), type = "error")
				} else{ 
					results$mix <- mix
					showNotification("Mixture spectrum files successfully loaded", type = "message")
				}
			}, ignoreInit = TRUE)


	# load a xlsx file with paths to spectra; these paths will be accessible for choosing 
	# in the path_form1, path_form2, path_mix input fields
	paths <- reactive({
				req(input$file_paths_bn)
				results$resultsFileLoaded <- FALSE
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
				updateSelectInput(session, "path_form1", choices = paths()[,1])
				updateSelectInput(session, "path_form2", choices = paths()[,1])
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
				dat <- rbind(c(NA, 0,1), dat); dat[1,1] <- "prop_form2"
				# make sure that constraints are in correct order
				dat <- dat[match(estim_order, dat$name), ]
				return(dat)      
			}, label = "load the user's parameter constraints for nloptr")
	
	# read in source data
	raw_data <- reactive({
				validate(need(results$form1, message = FALSE),
						need(results$form2, message = FALSE),
						need(results$mix, message = FALSE))
				
				form1 <- results$form1
				form2 <- results$form2
				mix <- results$mix
				
				info <- form2[[2]]
				spec_size <- info[1, "FTSIZE"] + input$add_zeroes
				ppm <- as.numeric(names(form2$data[[1]]))
				initial_ppm_range <- range(ppm)
				
				form1 <- form1$data[[1]]
				form2 <- form2$data[[1]]
				mix <- mix$data[[1]]
				
				# linear interpolation of mix and form1 spectra onto forms2's ppm values
				ppm_tmp <- as.numeric(names(form1))
				form1 <- complex(real = approx(x = ppm_tmp, y = Re(form1), xout = ppm)$y, 
				                 imaginary = approx(x = ppm_tmp, y = Im(form1), xout = ppm)$y)
				ppm_tmp <- as.numeric(names(mix))
				mix <- complex(real = approx(x = ppm_tmp, y = Re(mix), xout = ppm)$y, 
				                 imaginary = approx(x = ppm_tmp, y = Im(mix), xout = ppm)$y)
				df <-  bind_cols(ppm = ppm, data.frame(form1 = unclass(form1), form2 = unclass(form2), mix = unclass(mix)))
				
				if (!any(is.na(c(input$ppm_range1, input$ppm_range2)))) {
					if (input$ppm_range1 >= input$ppm_range2) {
						showNotification("The lower ppm boundary must be smaller than the upper boundary. The model will be fitted on the entire initial ppm range.", type = "warning") 
					} else {
						df <- subset(df, ppm >= input$ppm_range1 & ppm <= input$ppm_range2)
					}
				}
				
				return(list(df, info, spec_size, initial_ppm_range))
			}, label = "read in spectra")
	
	# update the pivot point field with the ppm value corresponding to 
	# the most abundant peak of the selected mixture spectrum
	max_peak_of_mix <- reactive({
	  tmp <- raw_data()
	  tmp[[1]]$ppm[which.max(abs(Re(tmp[[1]]$mix)))]
	  })
	
	observeEvent(results$mix, {
				tmp <- raw_data()
				req(tmp, cancelOutput = TRUE)
				updateNumericInput(session, inputId = "pivot_point", value = max_peak_of_mix())
			}, label = "update pivot point to the largest peak of mix", ignoreInit = TRUE)
	
	# process all three spectra at once (i.e. zero filling and apodization)
	proc_data <- reactive({
				raw_data <- raw_data()
				req(raw_data, cancelOutput = TRUE)
				
				ppm <- raw_data[[1]]$ppm
				form1 <- raw_data[[1]]$form1
				form2 <- raw_data[[1]]$form2
				mix <- raw_data[[1]]$mix
				info <- raw_data[[2]]
				spec_size <- raw_data[[3]]
				
				if (input$add_zeroes > 0 || input$lb_global > 0) {
				  
				  ppm <- ppm_zero_fill_apod(ppm = ppm, info = info, spec_size)
					form1 <- zero_fill_apod(form1, spec_size, input$lb_global, info[2])
					form2 <- zero_fill_apod(form2, spec_size, input$lb_global, info[2])
					mix <- zero_fill_apod(mix, spec_size, input$lb_global, info[2])
				} 
				
				# additional line broadening of the form2 template
				if (input$lb_form2 > 0) {
					form2 <- zero_fill_apod(form2, spec_size, input$lb_form2, info[2])
				}
				
				return(list(data.frame(ppm = ppm, form1 = form1, form2 = form2, mix = mix), info, spec_size))
			}, label = "process three spectra (zero fill, apodization)")
	
	fit_bn_rv <- reactiveValues(counter=0)
	observeEvent(input$fit_bn, {
	  fit_bn_rv$counter <- fit_bn_rv$counter+1
	}, label = "update the counter for the fitting button")
	
	# model fitting with nloptr
	model_fit <- eventReactive(input$fit_bn, {
				validate(need(results$form1, "form1 spectrum path required"), 
						need(results$form2, "form2 spectrum path required"), 
						need(results$mix, "Mixture spectrum path required"))

				notify_id <- showNotification("Model fitting is running...", duration = NULL, closeButton = FALSE, type = "message")
				on.exit(removeNotification(notify_id), add = TRUE)
				
				proc_data <- proc_data()
				
				ppm <- proc_data[[1]]$ppm
				form1 <- proc_data[[1]]$form1
				form2 <- proc_data[[1]]$form2
				mix <- proc_data[[1]]$mix
				
				if (input$estim_mode %in% c("only_prop", "explicit")) {
				  form1 <- shift_horizon(x = ppm, y = Re(form1), delta = input$ppm_form1)
				  form1 <- norm_sum(form1)
				  form2 <- norm_sum(form2)
				  N <- length(ppm)
				  lin_pred <- get_ph_angle(1:N, int = input$ph0_mix * to_rad_const, slope = input$ph1_mix * to_rad_const, pivot_point = input$pivot_point, n = N, ppm)
				  mix <- ph_corr(mix, lin_pred)
				  mix <- shift_horizon(x = ppm, y = Re(mix), delta = input$ppm_mix)
				  mix <- norm_sum(mix)
				  
				  model_input <- data.frame(ppm = ppm, form1 = form1, form2 = form2, mix = mix)
				  na_ind <- is.na(model_input$form1) | is.na(model_input$mix)
				  model_input <- model_input[!na_ind, ]
				  
				  if (input$estim_mode == "explicit") {
				    fitted <- input$prop_form2 * model_input$form2 + (1-input$prop_form2) * model_input$form1
				    residuals <- model_input$mix - fitted
				    dat <- cbind(model_input, fitted = fitted, residuals = residuals)
				    solution <- c(prop_form2 = input$prop_form2, rmse = sqrt(mean(residuals^2)))
				    return(list(dat = dat, solution = solution, 
				                # add below to be consistent with the other two modes of estimation
				                start = rlang::set_names(rep(0, length(estim_order)), estim_order)))
				  } else {
				    param_constraints_tmp <- param_constraints() %>% filter(name == "prop_form2")
				    param_start_tmp <- trim_values(input$prop_form2, param_constraints_tmp[, 2:3] %>% as.numeric())
				    mod <- nloptr_wrapper(
				      model_input
				      , x_order = "prop_form2"
				      , obj_fun = obj_fun
				      , param_start = param_start_tmp
				      , param_constraints = as.data.frame(param_constraints_tmp)
				      , optim_algorithm = input$optim_algorithm
				      , pivot_point = input$pivot_point)
				    return(mod)
				  }
				# proportion AND pre-processing parameters estimation below
				} else { 
				  ph0_angle <- input$ph0_mix
				  # if ph0_angle==0 then take pepsNMR PH0 as a starting value for PH0
				  if (input$do_PH0_PepsNMR) {
				    tmp <- matrix(mix, nrow = 1)
				    colnames(tmp) <- ppm
				    ph0_angle <- ZeroOrderPhaseCorrection(tmp, type.zopc = "max", returnAngle = TRUE)$Angle * to_deg_const
				  }
				  
				  model_input <- data.frame(ppm = ppm, form1 = form1, form2 = form2, mix = mix)
				  na_ind <- is.na(model_input$form1) | is.na(model_input$mix)
				  model_input <- model_input[!na_ind, ]
				  
				  param_constraints_tmp <- param_constraints()
				  # if user input values exceeds constraints specified for nloptr, then set starting values to the constraint values
				  inputs_list <- reactiveValuesToList(input)
				  param_start_tmp <- inputs_list[estim_order] 
				  param_start_tmp["ph0_mix"] <- ph0_angle
				  param_start_tmp <- param_start_tmp %>%
				    keep(~!is.null(.x)) %>% 
				    data.frame() %>% 
				    pivot_longer(cols = everything(), names_to = "name", values_to = "start") %>% 
				    inner_join(as.data.frame(param_constraints_tmp), by= "name") %>%
				    mutate(multi_fac = ifelse(stri_detect(name, fixed = "ph"), to_rad_const, 1),
				           start = start * multi_fac) %>%
				    select(-c(name, multi_fac)) %>% 
				    as.matrix() %>%
				    apply(1, function(x) {
				      trim_values(x[1], x[-1])
				    })
				  
				  mod <- nloptr_wrapper(
				    model_input
				    , x_order = estim_order
				    , obj_fun = obj_fun
				    , param_start = param_start_tmp
				    , param_constraints = as.data.frame(param_constraints_tmp)
				    , optim_algorithm = input$optim_algorithm
				    , pivot_point = input$pivot_point)
				  return(mod)
				}
			
			}, label = "do the fitting")
	
	## operations on the user inputs
	
	# after fitting with nloptr, update the corresponding inputs
	observeEvent(input$fit_bn, {
	      if (input$estim_mode == "only_prop") updateNumericInput(session, inputId = "prop_form2", value = as.numeric(model_fit()$solution["prop_form2"]))
	      else if (input$estim_mode == "prop_preproc") update_n_inputs(preproc_param_names, model_fit()$solution)
			}, label = "update user inputs with nloptr estimates")
	
	# reset all input parameters to their default values
	observeEvent({
	  input$reset_bn 
	  req(results$form1, results$form2, results$mix)
	  }, {
		  if(!results$resultsFileLoaded){
			  tmp_name <- c(preproc_param_names, adv_param_names[adv_param_names != "optim_algorithm"])
			  update_n_inputs(tmp_name, param_defaults %>% discard(is.character) %>% unlist())
			  updateNumericInput(session, inputId = "pivot_point", value = max_peak_of_mix())
			  updateTextInput(session, inputId = "optim_algorithm", value = param_defaults[["optim_algorithm"]])
			  updateNumericInput(session, inputId = "ppm_range1", value = param_defaults$ppm_range1)
			  updateNumericInput(session, inputId = "ppm_range2", value = param_defaults$ppm_range2)
		  } else results$resultsFileLoaded <- FALSE
			}, label = "reset all user inputs to the default values")
	
	# write the estimated parameters to a csv file
	track_inputs_rv <- reactiveValues(x = saved_params, count_rows = 0)
	
	observeEvent(c(results$form1, results$form2, results$mix), {
				track_inputs_rv$count_rows <- 0
			}, label = "reset row counts to zero after changing a spectrum")
	
	# each time Calculate button is pressed, update the underlying estimate tracking file
	observeEvent(input$fit_bn, {
				track_inputs_rv$count_rows <- track_inputs_rv$count_rows+1
				curr_input_vals <- reactiveValuesToList(input)
				next_row <- character(ncol(saved_params)); names(next_row) <- colnames(saved_params)
				# first, assign default values
				
				next_row[param_defaults_names] <- unlist(param_defaults)
				vals <- curr_input_vals[na.omit(match(param_defaults_names, names(curr_input_vals)))] %>%
						compact() %>% unlist()
				# next, update them with input values
				next_row[names(vals)] <- vals
				# finally, assign model estimated values
				ms <- model_fit()$solution
				next_row[c("id", "rmse", "prop_form2")] <- c(track_inputs_rv$count_rows, ms["rmse"], ms["prop_form2"])
				next_row[paste0(names(model_fit()$start), "_start")] <- as.numeric(model_fit()$start) * ifelse(stri_detect(names(model_fit()$start), fixed = "ph"), to_deg_const, 1)
				
				if (input$estim_mode == "only_prop") {
				  next_row["fit_type"] <- "only proportion"
				} else if (input$estim_mode == "explicit") {
				  next_row[c("fit_type", "optim_algorithm")] <- c("none, explicit proportion", "")
				} else if (input$estim_mode == "prop_preproc") {
				  next_row[c("fit_type", estim_order)] <- c("proportion and pre-processing parameters", ms[estim_order])
				}
				track_inputs_rv$x <- rbind(track_inputs_rv$x, next_row)
			}, label = "update the underlying estimate tracking file", ignoreInit = TRUE)
	
	# restore parameter values from the previous model fitting step
	observeEvent(input$undo_bn, {
				# this condition ensures that at least two records are in the results file
				if (fit_bn_rv$counter>=2) {
				  
				  nam <- param_defaults %>% discard(~is.character(.x)) %>% enframe() %>% unnest(cols = c(value)) %>% 
				    filter(!(name %in% "rmse") & !str_detect(name, pattern = "start")) %>%  pull(name)
				  prev_dat <- track_inputs_rv$x[nrow(track_inputs_rv$x)-1,]
				  update_n_inputs(nam, prev_dat)
				  updateTextInput(session, inputId = "optim_algorithm", value = as.character(prev_dat["optim_algorithm"]))
				  
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
				updateTextAreaInput(session, inputId = "user_comments", value = "")
			})
	
	# Re-initialize datasets at NULL when switching input source
	observeEvent(input$inputSource, {
				results$form1 <- NULL
				results$form2 <- NULL
				results$mix <- NULL
			})
	##
	
	# event for retaining zoom between data recalculations in plotly graph
	zoom <- reactive({
	  req(model_fit())
	  event_data("plotly_relayout", source = "p1")
	  })
	
	# event for retaining only selected lines between data recalculcations in the plotly graph
	legend_click <- reactive({
	  # to avoid warnings (https://github.com/ropensci/plotly/issues/1528)
	  event_data("plotly_legendclick", source = "p1")
	  })
	
	# register clicking event to update pivot point value
	plot_click <- reactive({
	  req(model_fit())
	  event_data("plotly_click", source = "p1")
	  })
	
	# showing/hiding lines by clicking on the legend 
	legend_items <- reactiveValues("form1 reference" = TRUE, "form2 reference" = TRUE, "mixture spectrum" = TRUE, "residuals" = TRUE, "fit" = TRUE)
	observe({
				if (!is.null(legend_click()$name)) {
					legend_items[[legend_click()$name]] <- ifelse(legend_click()$visible == TRUE, legend_states[[2]], legend_states[[1]])
				}
			})
	
	observeEvent(plot_click(), {
				updateNumericInput(session, inputId = "pivot_point", value = plot_click()$x)
			})
	
	plot_object <- eventReactive(input$fit_bn, {
	  req(model_fit())
	  p <- plot_model_fit(model_fit()) 
	  # the output of plotly_relayout changes depending on the action: if you zoom in
	  # xaxis.range[0], xaxis.range[1], yaxis.range[0],	yaxis.range[1] values will be available
	  # if you don't zoom in, the output will be NULL or some other named vectors of length different than 4
	  if (is.null(zoom()) | length(zoom()) != 4) {
	    p <- layout(p, xaxis = list(zeroline = FALSE, autorange = "reversed", title = "ppm"), yaxis = list(title = "intensity"))
	  } else {
	    p <- layout(p, xaxis = list(zeroline = FALSE, range = c(zoom()$"xaxis.range[0]", zoom()$"xaxis.range[1]"), title = "ppm"), yaxis = list(title = "intensity", range = c(zoom()$"yaxis.range[0]", zoom()$"yaxis.range[1]")))
	  }
	  
	  p <- plotly_build(p) %>%
	    event_register("plotly_legendclick")
	  for (i in seq_along(p$x$data)) {
	    p$x$data[[i]]$visible <- legend_items[[p$x$data[[i]]$name]]
	  }
	  return(p)
	}, label = "call plot function and modify legend and axis range")
	
	output$fit_plot_out <- renderPlotly({
	  req(plot_object())
	  plot_object()
			})
	
	output$fit_stats_out <- renderText({
				req(model_fit())
				paste0(
						"Estimated form2 proportion [0-1 interval] ", round(model_fit()$solution["prop_form2"], 7),
						"\nrmse (x10^6):", round(10^6 * model_fit()$solution["rmse"], 7),
						"\nInitial spectrum size:", isolate(raw_data()[[2]][1, "FTSIZE"]),
						"\nFinal spectrum size:", isolate(raw_data()[[2]][1, "FTSIZE"] + input$add_zeroes)
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
				results$resultsFileLoaded <- TRUE
				inFile <- input$file_load_bn
				tmp <- read.csv(inFile$datapath, stringsAsFactors = FALSE)
				
				# Since only the last row will be loaded in, the user has to manually remove other records.
				# The removal can be done in several ways: Delete, Backspace button on entire rows or selected cells.
				# As a consequence, sometimes the removed rows may still be loaded in, but they will contain NA values (of different types).
				# Heuristic: in such case remove rows with "empty" id, path_form1, path_form2, path_mix columns.
				var_check <- c("id", "path_form1", "path_form2", "path_mix")
				complete_cases_ind <- !apply(tmp, 1, function(x) all(is.na(x[var_check]) | x[var_check] == "" ))
				tmp <- tmp[complete_cases_ind,]
				dat <- tmp[nrow(tmp),]
				if(localMode || input$inputSource == "Shared drive"){
					if(!is.na(dat[["path_form1"]])) updateSelectInput(session, "path_form1", choices = dat[["path_form1"]])
					if(!is.na(dat[["path_form2"]])) updateSelectInput(session, "path_form2", choices = dat[["path_form2"]])
					if(!is.na(dat[["path_mix"]])) updateSelectInput(session, "path_mix", choices = dat[["path_mix"]])
				} else {
					showNotification("Note that file paths in results file cannot be re-read from local computer.", type = "warning") 
				}
				param_selected <- c(preproc_param_names, adv_param_names[adv_param_names != "optim_algorithm"], "ppm_range1", "ppm_range2")
				update_n_inputs(param_selected, dat[param_selected])
				updateNumericInput(session, inputId = "pivot_point", value = dat[["pivot_point"]]) # Also take this value from result file (in case no file paths are available)
				updateTextInput(session, inputId = "optim_algorithm", value = as.character(dat["optim_algorithm"]))
			}, label = "load in estimated results")
	
	
	output$localMode <- reactive({
				return(localMode)
			})
	outputOptions(output, 'localMode', suspendWhenHidden=FALSE)
	
}



