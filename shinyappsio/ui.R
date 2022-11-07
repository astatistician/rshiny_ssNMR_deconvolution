ui <- fluidPage(
  headerPanel("Automated deconvolution of solid-state NMR mixture spectra"),
    shinyFeedback::useShinyFeedback(),
  fluidRow(
    column(
      width = 4,
      wellPanel(
        strong(HTML("<p style='color:#2596be; font-size:15px'> Test the app with example data </p>")) %>% 
          helper(type = "inline",
                 title = "",
                 content = c("Press 'Load example' to access three example 19F ssNMR spectra and test the app.", 
                             "Press 'Download input file template (CSV)' to see the required structure of the CSV input files (note that spectra can also be uploaded in a JCAMP-DX format)."),
                 size = "l"),
        div(style = "display:inline-block", actionButton("load_example_bn", "Load example")),
        div(style = "display:inline-block", h4(" ")),
        div(style = "display:inline-block", downloadButton("download_example_data_bn", "Download input file template (CSV)"))
        ),
      wellPanel(
        strong(HTML("<p style='color:#2596be; font-size:15px'> Load your own spectra </p>")),
        "Bruker data folder saved as JCAMP-DX (.dx) or manually prepared CSV files" %>% 
          helper(type = "inline",
                 title = "",
                 content = c("To export Bruker spectra from TopSpin, open a spectrum, press ctrl+S and select 'Save data set in a JCAMP-DX file' with default values.",
                             "The form1 and form2 spectra should correspond to the reference spectra of pure forms of the compound of interest.",
                             "NMR acquisition and processing of form1, form2 and mixture signals should be ideally conducted in the same way.",
                             "The number of spectral data points MUST be identical."),
                 size = "l"),
        br(), br(),
        fileInput("files_form1", "Choose form1 files",
            multiple = TRUE,
            accept = c(".dx", ".csv")
          ),
        fileInput("files_form2", "Choose form2 files",
            multiple = TRUE,
            accept = c(".dx", ".csv")
          ),
          fileInput("files_mix", "Choose mixture files",
            multiple = TRUE,
            accept = c(".dx", ".csv")
          ),
        div(style = "display:inline-block", textInput("description_form1", label = "Form1 label", value = "", width = 100)),
        div(style = "display:inline-block", textInput("description_form2", label = "Form2 label", value = "", width = 100)),
        div(style = "display:inline-block", textInput("description_mix", label = "Mix label", value = "", width = 100) %>% 
              helper(type = "inline", title = "", size = "l", 
                     content = c("Optional labels displayed in the upper-left corner of the interactive graph.",
                     "These can be a short description, name, or a file path on a local computer that will create a meaningful mapping between the generic terms like 'form1', 'form2', 'mixture' and data identifier for traceability.")
                     )
            )),
      wellPanel(
        strong(HTML("<p style='color:#2596be; font-size:15px'> Spectral processing parameters</p>")),
        br(),
        HTML("<p style='font-size:15px'> Adjust the x-axis range:</p>"),
        div(style = "display:inline-block", numericInput(inputId = "ppm_range1", label = "Select lower ppm limit", value = param_defaults$ppm_range1, step = 1, width = 140)),
        div(style = "display:inline-block", h4(" ")),
        div(style = "display:inline-block", numericInput(inputId = "ppm_range2", label = "Select upper ppm limit", value = param_defaults$ppm_range2, step = 1, width = 140) %>%
              helper(type = "inline", title = "", size = "m", 
                     content = c("By default, ppm boundaries are set to invisible NA values, meaning no zoom into particular ppm region."))
            ),
        numericInput("ppm_form1", label = "Horizontal shift of form1 (ppm; >0 right, <0 left)", value = param_defaults$ppm_form1, min = -50, max = 50, step = 0.01),
        numericInput("ppm_mix", label = "Horizontal shift of mixture (ppm; >0 right, <0 left)", value = param_defaults$ppm_mix, min = -50, max = 50, step = 0.01),
        numericInput("ph0_mix", label = "PH0 of mixture (degrees)", value = param_defaults$ph0_mix, min = -180, max = 180, step = 0.01) %>%
          helper(type = "inline", title = "", size = "s", 
                 content = c("Zero-order phase correction.")), 
        numericInput("pivot_point", label = "Pivot point (ppm)", value = 0, min = -100000, max = 100000, step = 0.1) %>%
          helper(type = "inline", title = "", size = "m", 
                 content = c("Enter a value from keyboard or left-click on a certain point of a line on the graph.")),
        numericInput("ph1_mix", label = "PH1 of mixture (degrees)", value = param_defaults$ph1_mix, min = -180, max = 180, step = 0.01) %>% 
          helper(type = "inline", title = "", size = "s", 
                 content = c("First-order phase correction.")),
        numericInput("prop_form2", label = "form2 proportion value [0-1]", value = param_defaults$prop_form2, min = 0, max = 100, step = 0.001) %>% 
          helper(type = "inline", title = "", size = "m", 
                 content = c("form1 proportion is calculated as 1 - form2 proportion.")),
        textOutput("prop_form2_check"),
        
        shinyBS::bsCollapse(
          id = "show_adv_opt",
          shinyBS::bsCollapsePanel("Advanced options",
            style = "primary",
            numericInput("add_zeroes", label = "Number of additional zeroes", value = param_defaults$add_zeroes, min = 0, max = +Inf, step = 100) %>% 
              helper(type = "inline", title = "", size = "s", 
                     content = c("Otherwise known as zero filling or zero padding.")), 
            numericInput("lb_global", label = "Line broadening of all spectra (Hz)", value = param_defaults$lb_global, min = 0, max = +Inf, step = 0.1) %>% 
              helper(type = "inline", title = "", size = "s", 
                     content = c("Exponential function multiplication.")), 
            numericInput("lb_form1", label = "Line broadening of form1 spectrum (Hz)", value = param_defaults$lb_form1, min = 0, max = +Inf, step = 0.1),
            numericInput("lb_form2", label = "Line broadening of form2 spectrum (Hz)", value = param_defaults$lb_form2, min = 0, max = +Inf, step = 0.1),
            numericInput("lb_mix", label = "Line broadening of mixture spectrum (Hz)", value = param_defaults$lb_mix, min = 0, max = +Inf, step = 0.1),
            selectInput("optim_algorithm", "Select an optimization algorithm", selected = param_defaults$optim_algorithm, choices = optim_algorithms_list, selectize = FALSE),
            div(style = "display:inline-block", numericInput(inputId = "ph0_mix_lower", label = "PH0 lower", value = param_defaults$ph0_mix_lower, min = 0, max = 360, step = 0.01, width = 170)),
            div(style = "display:inline-block", numericInput(inputId = "ph0_mix_upper", label = "upper limit (degrees)", value = param_defaults$ph0_mix_upper, min = 0, max = 360, step = 0.01, width = 170) %>% 
                  helper(type = "inline", title = "", size = "s", 
                         content = c("Limits required by optimisation algorithm."))), 
            br(), br(),
            div(style = "display:inline-block", numericInput(inputId = "ph1_mix_lower", label = "PH1 lower", value = param_defaults$ph1_mix_lower, min = 2 * param_defaults$ph1_mix_lower, max = 0, step = 0.01, width = 170)),
            div(style = "display:inline-block", numericInput(inputId = "ph1_mix_upper", label = "upper limit (degrees)", value = param_defaults$ph1_mix_upper, min = 0, max = 2 * param_defaults$ph1_mix_upper, step = 0.01, width = 170) %>% 
                  helper(type = "inline", title = "", size = "s", 
                         content = c("Limits required by optimisation algorithm."))), 
            br(), br(),
            div(style = "display:inline-block", numericInput(inputId = "ppm_form1_lower", label = "Horizontal shift of form1 lower", value = param_defaults$ppm_form1_lower, min = -Inf, max = +Inf, step = 0.01, width = 170)),
            div(style = "display:inline-block", numericInput(inputId = "ppm_form1_upper", label = "upper limit (ppm)", value = param_defaults$ppm_form1_upper, min = -Inf, max = +Inf, step = 0.01, width = 170) %>% 
                  helper(type = "inline", title = "", size = "s", 
                         content = c("Limits required by optimisation algorithm."))), 
            br(), br(),
            div(style = "display:inline-block", numericInput(inputId = "ppm_mix_lower", label = "Horizontal shift of mixture lower", value = param_defaults$ppm_mix_lower, min = -Inf, max = +Inf, step = 0.01, width = 170)),
            div(style = "display:inline-block", numericInput(inputId = "ppm_mix_upper", label = "upper limit (ppm)", value = param_defaults$ppm_mix_upper, min = -Inf, max = +Inf, step = 0.01, width = 170) %>% 
                  helper(type = "inline", title = "", size = "s", 
                         content = c("Limits required by optimisation algorithm."))) 
        )),
        div(style = "display:inline-block", actionButton("reset_bn", "Reset parameters"))
      ),
      br(),
      wellPanel(
        div(style = "display:inline-block", radioButtons("estim_mode", "Choose the optimisation mode",
          selected = "prop_preproc",
          c(
            "no optimisation, apply the fixed processing values" = "explicit",
            "optimise only proportion" = "only_prop",
            "optimise proportion and other processing parameters" = "prop_preproc"
          )
        ) %>% 
          helper(type = "inline", title = "", size = "l", 
                 content = c("1) 'no optimisation, apply the fixed processing values': the input field values in the 'Spectral processing parameters' module are directly applied to process the spectra and calculate the fitted spectrum.",
                                          "2) 'optimise only proportion': the form2 proportion is estimated from data; the other parameter values are utilised to process the spectra.",
                                          "3) 'optimise proportion and other processing parameters': form2 proportion, PH0, PH1 and horizontal shifts are estimated from data.",
                                          "*After using the 2) or 3) optimisation, the estimated values are immediately displayed in the corresponding input fields in 'Spectral processing parameters'.",
                                          "**Starting values for optimisation parameters in 2) and 3) are obtained from inputs in 'Spectral processing parameters'."))),
        conditionalPanel(
          condition = "input.estim_mode=='prop_preproc'",
          checkboxInput("do_PH0_PepsNMR", label = "Initial PH0 from PepsNMR (recommended for mixture spectra with large phase errors)", value = TRUE)
        ),
        actionButton("fit_bn", "Fit the model!", width = 210)
      ),
      br(), br(),
      wellPanel(
        strong(HTML("<p style='color:#2596be;font-size:15px'> Misc </p>")),
        div(style = "display:inline-block", actionButton("undo_bn", "Restore the previous model fit") %>% 
              helper(type = "inline", title = "", size = "l", 
                     content = c("This will update the respective input fields in 'Spectral processing parameters' with values obtained from the previous model fitting.",
                                 "Select 'no optimisation, apply the fixed processing values' and 'Fith the model! to recreate the graph from the previous model fitting."))),
        br(), br(),
        div(style = "display:inline-block", downloadButton("download_params_bn", "Download estimated parameters") %>% 
              helper(type = "inline", title = "", size = "l", 
                     content = c("Download the tracking CSV file containing inputs and outputs of each round of performed model fitting.",
                                 "The last row of this file can be loaded in via 'Load estimated parameters from a file'.", 
                                 "Then, these values are shown in the respective fields in 'Spectral processing parameters' and can be used to recreate the graph."))),
        div(style = "display:inline-block", h4(" ")),
        div(style = "display:inline-block", downloadButton("download_spectral_data_bn", "Download outcome spectra") %>% 
              helper(type = "inline", title = "", size = "m", 
                     content = c("Download the CSV file with ppm and intensity values of the form1, form2, mixture, fitted and residual spectra after the model fitting."))),
        br(), br(),
        fileInput("file_load_bn", "Load estimated parameters from a file", accept = c(".csv")) %>% 
          helper(type = "inline", title = "", size = "m", 
                 content = c("Load the earlier saved result file ('Download estimated parameters' button) to resume the previous analysis.",
                             "Note that you still have to manually supply the input spectra."))
      )
    ),
    column(
      width = 8,
      HTML("version 1.0.0 (2022-09-19) - Valkenborg Lab (contact: piotr.prostko@uhasselt.be)"),
      plotlyOutput("fit_plot_out", height = 600),
      verbatimTextOutput("fit_stats_out", placeholder = FALSE),
      # conditionalPanel(
      #   condition = "input.fit_bn_rv>0",
      #   br(),
      #   textAreaInput(inputId = "user_comments", label = "Enter comments about the fit", value = ""),
      #   br()
      # ),
      # conditionalPanel(
      #   condition = "input.fit_bn_rv>0",
      #   br(),
      #   actionButton(inputId = "save_user_comments_bn", label = "Save the comments"),
      #   br()
      # ),
      # conditionalPanel(
      #   condition = "input.fit_bn_rv>0",
      #   br(),
      #   actionButton("approve_bn", "Approve the fit"),
      #   br()
      # )
      uiOutput("user_comments"),
      uiOutput("save_user_comments_bn"),
      uiOutput("approve_bn")
    )
  )
)