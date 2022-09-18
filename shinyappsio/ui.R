ui <- fluidPage(
  headerPanel("Automated deconvolution of solid-state NMR mixture spectra"),
  fluidRow(
    column(
      width = 4,
      wellPanel(
        strong(HTML("<p style='color:#2596be; font-size:15px'> Test the app with example data </p>")),
        div(style = "display:inline-block", actionButton("load_example_bn", "Load example")),
        div(style = "display:inline-block", h4(" ")),
        div(style = "display:inline-block", downloadButton("download_example_data_bn", "Download input file template"))
      ),
      wellPanel(
        strong(HTML("<p style='color:#2596be; font-size:15px'> Load your own spectra </p>")),
        fileInput("files_form1", "Choose form1 csv files (intensity and params)",
          multiple = TRUE
        ),
        fileInput("files_form2", "Choose form2 csv files (intensity and params)",
          multiple = TRUE
        ),
        fileInput("files_mix", "Choose mixture csv files (intensity and params)",
          multiple = TRUE
        ) %>% 
          bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                           bsplus::bs_embed_tooltip(title = "NMR acquisition and processing of form1, form2 and mixture signals should be ideally conducted in the same way. The number of spectral data points must be identical.", placement = "left"))
      ),
      wellPanel(
        strong(HTML("<p style='color:#2596be; font-size:15px'> Spectral processing parameters</p>")),
        br(),
        HTML("<p style='font-size:15px'> Adjust the x-axis range:</p>"),
        div(style = "display:inline-block", numericInput(inputId = "ppm_range1", label = "Select lower ppm limit", value = param_defaults$ppm_range1, step = 1, width = 140)),
        div(style = "display:inline-block", h4(" ")),
        div(style = "display:inline-block", numericInput(inputId = "ppm_range2", label = "Select upper ppm limit", value = param_defaults$ppm_range2, step = 1, width = 140) %>% 
              bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                               bsplus::bs_embed_tooltip(title = "By default, ppm boundaries are set to invisible NA values, meaning no zoom into particular ppm region.", placement = "left"))),
        numericInput("ppm_form1", label = "Horizontal shift of form1 (ppm; >0 right, <0 left)", value = param_defaults$ppm_form1, min = -50, max = 50, step = 0.01),
        numericInput("ppm_mix", label = "Horizontal shift of mixture (ppm; >0 right, <0 left)", value = param_defaults$ppm_mix, min = -50, max = 50, step = 0.01),
        numericInput("ph0_mix", label = "PH0 of mixture (degrees)", value = param_defaults$ph0_mix, min = -180, max = 180, step = 0.01) %>% 
          bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                           bsplus::bs_embed_tooltip(title = "Zero-order phase correction", placement = "left")),
        numericInput("pivot_point", label = "Pivot point (ppm)", value = 0, min = -100000, max = 100000, step = 0.1) %>% 
          bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                           bsplus::bs_embed_tooltip(title = "Enter a value from keyboard or left-click on a certain point of a line on the graph.", placement = "left")),
        numericInput("ph1_mix", label = "PH1 of mixture (degrees)", value = param_defaults$ph1_mix, min = -180, max = 180, step = 0.01) %>% 
        bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                         bsplus::bs_embed_tooltip(title = "First-order phase correction", placement = "left")),
        numericInput("prop_form2", label = "form2 proportion value [0-1]", value = param_defaults$prop_form2, min = 0, max = 100, step = 0.001) %>% 
          bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                           bsplus::bs_embed_tooltip(title = "form1 proportion is then calculated as 1 - form2 proportion.", placement = "left")),
        
        shinyBS::bsCollapse(
          id = "show_adv_opt",
          shinyBS::bsCollapsePanel("Advanced options",
            style = "primary",
            numericInput("add_zeroes", label = "Number of additional zeroes", value = param_defaults$add_zeroes, min = 0, max = +Inf, step = 100) %>% 
              bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                               bsplus::bs_embed_tooltip(title = "Otherwise known as zero filling or zero padding.", placement = "left")),
            numericInput("lb_global", label = "Line broadening of all spectra (Hz)", value = param_defaults$lb_global, min = 0, max = +Inf, step = 0.1) %>% 
              bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                               bsplus::bs_embed_tooltip(title = "Exponential function multiplication.", placement = "left")),
            numericInput("lb_form1", label = "Line broadening of form1 spectrum (Hz)", value = param_defaults$lb_form1, min = 0, max = +Inf, step = 0.1),
            numericInput("lb_form2", label = "Line broadening of form2 spectrum (Hz)", value = param_defaults$lb_form2, min = 0, max = +Inf, step = 0.1),
            numericInput("lb_mix", label = "Line broadening of mixture spectrum (Hz)", value = param_defaults$lb_mix, min = 0, max = +Inf, step = 0.1),
            selectInput("optim_algorithm", "Select an optimization algorithm", selected = param_defaults$optim_algorithm, choices = optim_algorithms_list, selectize = FALSE),
            div(style = "display:inline-block", numericInput(inputId = "ph0_mix_lower", label = "PH0 lower", value = param_defaults$ph0_mix_lower, min = 0, max = 360, step = 0.01, width = 170)),
            div(style = "display:inline-block", numericInput(inputId = "ph0_mix_upper", label = "upper limit (degrees)", value = param_defaults$ph0_mix_upper, min = 0, max = 360, step = 0.01, width = 170) %>% 
              bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                               bsplus::bs_embed_tooltip(title = "Limits required by optimisation algorithm.", placement = "left"))),
            br(), br(),
            div(style = "display:inline-block", numericInput(inputId = "ph1_mix_lower", label = "PH1 lower", value = param_defaults$ph1_mix_lower, min = 2 * param_defaults$ph1_mix_lower, max = 0, step = 0.01, width = 170)),
            div(style = "display:inline-block", numericInput(inputId = "ph1_mix_upper", label = "upper limit (degrees)", value = param_defaults$ph1_mix_upper, min = 0, max = 2 * param_defaults$ph1_mix_upper, step = 0.01, width = 170) %>% 
              bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                               bsplus::bs_embed_tooltip(title = "Limits required by optimisation algorithm.", placement = "left"))),
            br(), br(),
            div(style = "display:inline-block", numericInput(inputId = "ppm_form1_lower", label = "Horizontal shift of form1 lower", value = param_defaults$ppm_form1_lower, min = -Inf, max = +Inf, step = 0.01, width = 170)),
            div(style = "display:inline-block", numericInput(inputId = "ppm_form1_upper", label = "upper limit (ppm)", value = param_defaults$ppm_form1_upper, min = -Inf, max = +Inf, step = 0.01, width = 170) %>% 
              bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                               bsplus::bs_embed_tooltip(title = "Limits required by optimisation algorithm.", placement = "left"))),
            br(), br(),
            div(style = "display:inline-block", numericInput(inputId = "ppm_mix_lower", label = "Horizontal shift of mixture lower", value = param_defaults$ppm_mix_lower, min = -Inf, max = +Inf, step = 0.01, width = 170)),
            div(style = "display:inline-block", numericInput(inputId = "ppm_mix_upper", label = "upper limit (ppm)", value = param_defaults$ppm_mix_upper, min = -Inf, max = +Inf, step = 0.01, width = 170) %>% 
                  bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                                   bsplus::bs_embed_tooltip(title = "Limits required by optimisation algorithm.", placement = "right")))
          )
        ),
        div(style = "display:inline-block", actionButton("reset_bn", "Reset parameters"))
      ),
      br(),
      wellPanel(
        div(style = "display:inline-block", radioButtons("estim_mode", "Choose the optimisation mode",
          selected = "only_prop",
          c(
            "no optimisation, apply the fixed processing values" = "explicit",
            "optimise only proportion" = "only_prop",
            "optimise proportion and other processing parameters" = "prop_preproc"
          )
        ) %>% 
          bsplus::shinyInput_label_embed(bsplus::shiny_iconlink() %>%
                                           bsplus::bs_embed_tooltip(title = 
                                          "1)'no optimisation, apply the fixed processing values': the input field values in the 'Spectral processing parameters' module are directly applied to process the spectra and calculate the fitted spectrum.
                                          \n2) 'optimise only proportion': the form2 proportion is estimated from data.
                                          \n3) 'optimise proportion and other processing parameters': form2 proportion, PH0, PH1 and horizontal shifts are estimated from data.
                                          \n* After using the 2) or 3) optimisation, the estimated values are immediately displayed in the corresponding input fields in 'Spectral processing parameters'.
                                          \n** Starting values for optimisation parameters in 2) and 3) are obtained from inputs in 'Spectral processing parameters'.", placement = "bottom"))),
        conditionalPanel(
          condition = "input.estim_mode=='prop_preproc'",
          checkboxInput("do_PH0_PepsNMR", label = "Initial PH0 from PepsNMR (recommended for mixture spectra with large phase errors)", value = TRUE)
        ),
        actionButton("fit_bn", "Fit the model!", width = 210)
      ),
      br(), br(),
      wellPanel(
        strong(HTML("<p style='color:#2596be;font-size:15px'> Misc </p>")),
        div(style = "display:inline-block", actionButton("undo_bn", "Restore the previous model fit")),
        br(), br(),
        div(style = "display:inline-block", downloadButton("download_params_bn", "Download estimated parameters")),
        div(style = "display:inline-block", h4(" ")),
        div(style = "display:inline-block", downloadButton("download_spectral_data_bn", "Download outcome spectra")),
        br(), br(),
        fileInput("file_load_bn", "Load estimated parameters from a file", accept = c(".csv"))
      ),
    ),
    column(
      width = 8,
      plotlyOutput("fit_plot_out", height = 600),
      verbatimTextOutput("fit_stats_out", placeholder = FALSE),
      textAreaInput(inputId = "user_comments", label = "Enter comments about the fit", value = ""),
      actionButton(inputId = "save_user_comments_bn", label = "Save the comments"),
      actionButton("approve_bn", "Approve the fit")
    )
  )
)
