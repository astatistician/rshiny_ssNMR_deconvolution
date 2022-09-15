ui <- fluidPage(
		headerPanel("Automated deconvolution of solid-state NMR mixture spectra"),
		fluidRow(
				column(
						width = 3,
						wellPanel(
										fileInput('files_form1', 'Choose form1 csv files (intensity and params)',
												multiple = TRUE),
										fileInput('files_form2', 'Choose form2 csv files (intensity and params)',
												multiple = TRUE),
										fileInput('files_mix', 'Choose mixture csv files (intensity and params)',
												multiple = TRUE),
										div(style="display:inline-block", actionButton("load_example_bn", "Load example data", width = 210)),
										div(style="display:inline-block", downloadButton("download_example_data_bn", "Download example data", width = 210))

						),
						wellPanel(
								div(style="display:inline-block", numericInput(inputId = "ppm_range1", label = "Select lower ppm boundary", value = param_defaults$ppm_range1, step = 1, width = 140) %>% 
						          bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
						                                           bsplus::bs_embed_tooltip(title = "By default, ppm boundaries are set to invisible NA values, meaning no zoom into particular ppm region.", placement ="left")
						          )),
				  				div(style="display:inline-block", h4(" ")),
								div(style="display:inline-block", h4(" ")),
								div(style="display:inline-block", h4(" ")),
				  				div(style="display:inline-block",numericInput(inputId = "ppm_range2", label = "Select upper ppm boundary", value = param_defaults$ppm_range2, step = 1, width = 140)),
								numericInput("ppm_form1", label = "Chemical shift (ppm) of form1 (>0 right, <0 left)", value = param_defaults$ppm_form1, min = -50, max = 50, step = 0.01),
								numericInput("ppm_mix", label = "Chemical shift (ppm) of mixture (>0 right, <0 left)", value = param_defaults$ppm_mix, min = -50, max = 50, step = 0.01),
								numericInput("ph0_mix", label = "PH0 (degrees) of mixture", value = param_defaults$ph0_mix, min = -180, max = 180, step = 0.01),
								numericInput("pivot_point", label = "Pivot point (ppm)", value = 0, min = -100000, max = 100000, step = 0.1),
								numericInput("ph1_mix", label = "PH1 (degrees) of mixture", value = param_defaults$ph1_mix, min = -180, max = 180, step = 0.01),
								numericInput("prop_form2", label = "form2 proportion value [0-1]", value = param_defaults$prop_form2, min = 0, max = 100, step = 0.001),
						),
						
						wellPanel(
								shinyBS::bsCollapse(id = "show_adv_opt",
										shinyBS::bsCollapsePanel("Advanced options",
												style = "primary",
												numericInput("add_zeroes", label = "Number of additional zeroes", value = param_defaults$add_zeroes, min = 0, max = +Inf, step = 100),
												numericInput("lb_global", label = "Line broadening for each spectrum (Hz)", value = param_defaults$lb_global, min = 0, max = +Inf, step = 0.1),
												numericInput("lb_form1", label = "Line broadening for form1 spectrum (Hz)", value = param_defaults$lb_form1, min = 0, max = +Inf, step = 0.1),
												numericInput("lb_form2", label = "Line broadening for form2 spectrum (Hz)", value = param_defaults$lb_form2, min = 0, max = +Inf, step = 0.1),
												numericInput("lb_mix", label = "Line broadening for mixture spectrum (Hz)", value = param_defaults$lb_mix, min = 0, max = +Inf, step = 0.1),
												selectInput("optim_algorithm", "Select an optimization algorithm", selected= param_defaults$optim_algorithm, choices=optim_algorithms_list,  selectize=FALSE),
												div(style="display:inline-block",numericInput(inputId="ph0_mix_lower", label="Lower PH0 (degrees)", value = param_defaults$ph0_mix_lower, min = 0, max = 360, step = 0.01, width = 145)),
												div(style="display:inline-block",numericInput(inputId="ph0_mix_upper", label="Upper PH0", value = param_defaults$ph0_mix_upper, min = 0, max = 360, step = 0.01, width = 145)),
												br(),
												div(style="display:inline-block",numericInput(inputId="ph1_mix_lower", label="Lower PH1 (degrees)", value = param_defaults$ph1_mix_lower, min = 2 * param_defaults$ph1_mix_lower, max = 0, step = 0.01, width = 145)),
												div(style="display:inline-block",numericInput(inputId="ph1_mix_upper", label="Upper PH1", value = param_defaults$ph1_mix_upper, min = 0, max = 2 * param_defaults$ph1_mix_upper, step = 0.01, width = 145)),
												br(),
												div(style="display:inline-block",numericInput(inputId="ppm_form1_lower", label="Lower ppm shift form1", value = param_defaults$ppm_form1_lower, min = -Inf, max = +Inf, step = 0.01, width = 145)),
												div(style="display:inline-block",numericInput(inputId="ppm_form1_upper", label="Upper ppm shift form1", value = param_defaults$ppm_form1_upper, min = -Inf, max = +Inf, step = 0.01, width = 145)),
												br(),
												div(style="display:inline-block",numericInput(inputId="ppm_mix_lower", label="Lower ppm shift mixture", value = param_defaults$ppm_mix_lower, min = -Inf, max = +Inf, step = 0.01, width = 145)),
												div(style="display:inline-block",numericInput(inputId="ppm_mix_upper", label="Upper ppm shift mixture", value = param_defaults$ppm_mix_upper, min = -Inf, max = +Inf, step = 0.01, width = 145))
										)
								)
						),
						br(),
						
						wellPanel(
						  div(style="display:inline-block", radioButtons("estim_mode", "Choose estimation mode", selected = "only_prop",
						                                                 c("none, explicit proportion" = "explicit",
						                                                   "only proportion" = "only_prop",
						                                                   "proportion and pre-processing parameters" = "prop_preproc")
						                                                 )
						      ),
						  conditionalPanel(
						    condition = "input.estim_mode=='prop_preproc'",
						    checkboxInput("do_PH0_PepsNMR", label = "Initial PH0 from PepsNMR (recommended for unprocessed mixture spec)", value = TRUE)
						  ),
						  actionButton("fit_bn", "Calculate", width = 210)
						  ),
						br(),br(),
						actionButton("reset_bn", "Reset parameters", width = 210),
						br(),br(),
						actionButton("undo_bn", "Restore the previous model fit", width = 210),
						br(),br(),
						downloadButton("download_params_bn", "Download results", width = 210),
						br(),br(),
						downloadButton("download_spectral_data_bn", "Download spectral data", width = 210),
						br(),br(),
						fileInput('file_load_bn', 'Load results from a file',accept = c(".csv"))
				),
				column(
						width = 9,
						plotlyOutput("fit_plot_out", height = 600),
						verbatimTextOutput("fit_stats_out", placeholder = FALSE),
						textAreaInput(inputId = "user_comments", label = "Enter comments about the fit", value = ""),
						actionButton(inputId = "save_user_comments_bn", label = "Save the comments"), 
						actionButton("approve_bn", "Approve the fit")
				)
		)
)