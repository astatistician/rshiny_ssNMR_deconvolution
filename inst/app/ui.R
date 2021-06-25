ui <- fluidPage(
		headerPanel("Semi-automated pre-processing and quantitative evaluation of ssNMR mixture spectra"),
		fluidRow(
				column(
						width = 3,
						wellPanel(
								fileInput('file_paths_bn', 'Choose a xlsx file with paths to spectra',
										accept = c(".xlsx")
								),
								selectInput("path_amo", "Select the amorphous spectrum path", choices="",  selectize=FALSE),
								selectInput("path_cr", "Select the crystalline spectrum path", choices="",  selectize=FALSE),
								selectInput("path_mix", "Select the mixture spectrum path", choices="",  selectize=FALSE)
						),
						wellPanel(
								numericInput("ppm_amo", label = "Chemical shift (ppm) of amorphous", value = param_defaults$ppm_amo, min = -50, max = 50, step = 0.01),
								numericInput("ppm_mix", label = "Chemical shift (ppm) of mixture", value = param_defaults$ppm_mix, min = -50, max = 50, step = 0.01),
								numericInput("ph0_mix", label = "PH0 (degrees) of mixture", value = param_defaults$ph0_mix, min = -180, max = 180, step = 0.01),
								numericInput("pivot_point", label = "Pivot point (ppm)", value = 0, min = -100000, max = 100000, step = 0.1),
								numericInput("ph1_mix", label = "PH1 (degrees) of mixture", value = param_defaults$ph1_mix, min = -180, max = 180, step = 0.01)
						),
						
						wellPanel(
								shinyBS::bsCollapse(id = "show_adv_opt",
										shinyBS::bsCollapsePanel("Advanced options",
												style = "primary",
												numericInput("add_zeroes", label = "Number of additional zeroes", value = param_defaults$add_zeroes, min = 0, max = +Inf, step = 100),
												numericInput("lb_global", label = "Line broadening for each spectrum (Hz)", value = param_defaults$lb_global, min = 0, max = +Inf, step = 0.1),
												numericInput("lb_cr", label = "Line broadening for crystalline spectrum (Hz)", value = param_defaults$lb_cr, min = 0, max = +Inf, step = 0.1),
												numericInput("prop_cr_start", label = "Starting value for crystalline proportion", value = param_defaults$prop_cr_start, min = 0, max = 1, step = 0.01),
												selectInput("optim_algorithm", "Select an optimization algorithm", selected= param_defaults$optim_algorithm, choices=optim_algorithms_list,  selectize=FALSE),
												div(style="display:inline-block",numericInput(inputId="ph0_mix_lower", label="Lower PH0 (degrees)", value = param_defaults$ph0_mix_lower, min = 0, max = 360, step = 0.01, width = 145)),
												div(style="display:inline-block",numericInput(inputId="ph0_mix_upper", label="Upper PH0", value = param_defaults$ph0_mix_upper, min = 0, max = 360, step = 0.01, width = 145)),
												br(),
												div(style="display:inline-block",numericInput(inputId="ph1_mix_lower", label="Lower PH1 (degrees)", value = param_defaults$ph1_mix_lower, min = 0, max = 360, step = 0.01, width = 145)),
												div(style="display:inline-block",numericInput(inputId="ph1_mix_upper", label="Upper PH1", value = param_defaults$ph1_mix_upper, min = 0, max = 360, step = 0.01, width = 145)),
												br(),
												div(style="display:inline-block",numericInput(inputId="ppm_amo_lower", label="Lower ppm shift amorphous", value = param_defaults$ppm_amo_lower, min = -50, max = 50, step = 0.01, width = 145)),
												div(style="display:inline-block",numericInput(inputId="ppm_amo_upper", label="Upper ppm shift amorphous", value = param_defaults$ppm_amo_upper, min = -50, max = 50, step = 0.01, width = 145)),
												br(),
												div(style="display:inline-block",numericInput(inputId="ppm_mix_lower", label="Lower ppm shift mixture", value = param_defaults$ppm_mix_lower, min = -50, max = 50, step = 0.01, width = 145)),
												div(style="display:inline-block",numericInput(inputId="ppm_mix_upper", label="Upper ppm shift mixture", value = param_defaults$ppm_mix_upper, min = -50, max = 50, step = 0.01, width = 145))
										)
								)
						),
						br(),
						actionButton("manual_fit_bn", "Manual fitting", width = 210),
						br(),br(),
						actionButton("automated_fit_bn", "Automated fitting", width = 210),
						br(),br(),
						actionButton("reset_bn", "Reset parameters", width = 210),
						br(),br(),
						actionButton("undo_bn", "Restore the previous model fit", width = 210),
						br(),br(),
						downloadButton("download_params_bn", "Download results", width = 210),
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



