ui <- fluidPage(
    titlePanel("SCCM: Single-Case Causal Mediation Analysis"),
    withMathJax(),
    sidebarLayout(
        sidebarPanel(
            fileInput("file_dat",
                      h3("choose a CSV, RDS, or XLSX file"),
                      multiple = FALSE,
                      accept = c(".csv", ".CSV", ".rds", ".RDS", ".xlsx", ".XLSX")),
            textOutput("out_warning"),
            tags$head(tags$style("#out_warning{color: red;}")),
            selectInput("select_x",
                      label = h3("X"),
                      choices = NULL),
            selectInput("select_m",
                        label = h3("M"),
                        choices = NULL),
            selectInput("select_y",
                        label = h3("Y"),
                        choices = NULL),
            numericInput("txt_lag",
                        label = h3("lag"),
                        value = 1,
                        min = 1,
                        max = 99,
                        step = 1),
            br(),
            h3("permutation test"),
            checkboxInput("check_perm",
                          label = "yes",
                          value = FALSE,
                          width = NULL),
            numericInput("txt_reps",
                         label = h3("permutation reps"),
                         value = 1000,
                         min = 100,
                         max = Inf,
                         step = 100),
            actionButton("btn_run", "run!")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel(
                    "results",
                    verbatimTextOutput("out_results")
                ),
                tabPanel(
                    "path diagram",
                    plotOutput("out_diagram")
                ),
                tabPanel(
                    "lavaan output",
                    verbatimTextOutput("out_lav_output")
                ),
                tabPanel(
                    "lavaan syntax",
                    verbatimTextOutput("out_lav_syntax")
                )
            )
        )
    )
)

#' @keywords internal
#' @importFrom tools file_ext
#' @importFrom openxlsx read.xlsx
server <- function(input, output, session) {

    dat <- NULL
    mod <- NULL
    warning <- NULL
    dats <- data_sets

    observe({
        in_file <- input$file_dat
        if (!is.null(in_file) && file.exists(file_path <- in_file$datapath)) {
            if (tolower(file_ext(file_path)) == "csv") {
                dat <<- read.csv(file_path)
                warning <<- NULL
            } else if (tolower(file_ext(file_path)) == "rds") {
                dat <- readRDS(file_path)
                warning <<- NULL
            } else if (tolower(file_ext(file_path)) == "xlsx") {
                dat <- read.xlsx(file_path)
                warning <<- NULL
            } else {
                dat <<- NULL
                warning <<- "Unknown file type."
            }
            updateSelectInput(session,
                              "select_x",
                              choices = names(dat),
                              selected = names(dat)[1])
            updateSelectInput(session,
                              "select_m",
                              choices = names(dat),
                              selected = names(dat)[2])
            updateSelectInput(session,
                              "select_y",
                              choices = names(dat),
                              selected = names(dat)[3])
        }

    })

    observe({
        input$btn_run

        if (!is.null(dat)) {
            mod <<- sccm(
                X = isolate(input$select_x),
                M = isolate(input$select_m),
                Y = isolate(input$select_y),
                dat = dat,
                lag = isolate(input$txt_lag),
                permutation = isolate(input$check_perm),
                perm_reps = isolate(input$txt_reps)
            )
        } else {
            mod <<- NULL
        }

    })

    output$out_results <- renderPrint({
        input$btn_run

        if (!is.null(mod)) {
            summary(mod)
        }
    })
    output$out_diagram <- renderPlot({
        input$btn_run

        if (!is.null(mod)) {
            get_plot(mod)
        }
    })
    output$out_lav_output <- renderPrint({
        input$btn_run

        if (!is.null(mod)) {
            print(mod)
        }
    })
    output$out_lav_syntax <- renderPrint({
        input$btn_run

        if (!is.null(mod)) {
            cat(mod$syntax)
        }
    })

    output$out_warning <- renderText({
        input$file_dat
        if (!is.null(warning)) {
            paste0("Warning: ", warning)
        } else {
            ""
        }
    })
}
