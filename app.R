# Load packages ----
#library(shiny)

# Source helpers ----
source('helper.R')

# Define UI for EigenGWAS Application
ui <- fluidPage(
  #Title
  titlePanel("EigenGWAS"),
  #Gobal parameters
  navbarMenu(
    title = "",
    fileInput('file_input', 'Genotype files', multiple = TRUE, accept = c("bed", "fam", "bim")),
    numericInput(
      'espace',
      'Eigen space',
      value = 5,
      min = 1,
      max = 10,
      step = 1
    ),
    radioButtons(
      'bred',
      'Population type',
      choices = list('Inbred' = 'inbred', 'Outbred' = 'outbred'),
      selected = 'inbred',
      inline = T
    ),
    numericInput(
      'threshhold',
      'p-value threshhold',
      value = 5.00,
      min = 1,
      max = 10,
      step = 0.1
    ),
    actionButton("run", "EigenGWAS, Go!")
  ),

  #TabPanels
  mainPanel(
    tabsetPanel(
      type = 'tabs',
      tabPanel('GRM',
               sidebarPanel(actionButton(
                 'GRM_save', 'Download'
               )),
               mainPanel(plotOutput('grm'))),
      tabPanel(
        'Eigenvalue',
        sidebarPanel(actionButton('Eigenvalue_save', 'Download')),
        mainPanel(plotOutput('eigenvalue'))
      ),
      tabPanel('Summary', 'SUMMARY'),
      tabPanel('Data table', 'DATATABLE'),
      tabPanel(
        'EigenGWAS',
        tabsetPanel(
          type = 'tabs',
          tabPanel(
            'miamiPlot',
            sidebarPanel(
              #fileUpload
              #fileInput('miamiPlot_integer_file',label = h4("Upload your file pls")),
              #args
              sliderInput(
                'miamiPlot_integer',
                "Integer:",
                min = 1,
                max = 5,
                value = 1,
                step = 1
              ),
              sliderInput(
                'miamiPlot_cex',
                'Cex:',
                min = 0.05,
                max = 1,
                value = 0.5,
                step = 0.05
              ),
              sliderInput(
                'miamiPlot_pch',
                'Pch:',
                min = 8,
                max = 24,
                value = 16,
                step = 2
              ),
              #image-save
              actionButton('miamiPlot_save', 'Download')
            ),
            mainPanel(plotOutput('miami'))
          ),
          tabPanel(
            'EigenGWASPlot',
            sidebarPanel(
              sliderInput(
                'EigenGWASPlot_integer',
                "Integer:",
                min = 1,
                max = 5,
                value = 1,
                step = 1
              ),
              actionButton('EigenGWASPlot_save', 'Download')
            ),
            mainPanel(plotOutput('eigengwas'))
          ),
          tabPanel('Summary',
                   mainPanel(plotOutput(
                     'eigengwas_summary'
                   )))
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #Plot on the web
  observeEvent(input$run, {
    cat(getwd())
    cat("\nShowing", input$file_input$name, " ", input$file_input$datapath)
    idx=0
    if (length(input$file_input$name)!=3) {
      return()
    } else {
      idx=grep(".bed$", input$file_input$datapath)
      if ( length(idx)==1 ) {
        rt=substr(input$file_input$datapath[idx], 1, nchar(input$file_input$datapath[idx])-4)
      }
      for (i in 1:3) {
        if (i != idx) {
          f1 = input$file_input$datapath[i]
          tl = substr(f1, nchar(f1)-2, nchar(f1))
          file.symlink(f1, paste0(rt, ".", tl))
        }
      }
    }
    froot = substr(input$file_input$datapath[idx], 1, nchar(input$file_input$datapath[idx])-4)
    RunEigenGWAS(froot, input$espace, inbred = ifelse(input$bred == 'inbred', T, F))
    cat("\ndone EigenGWAS")
    #GRM-plotOutput
    output$grm <- renderPlot({
      layout(matrix(1:2, 1, 2))
      grmStats(froot)
    })

    #Eigenvalue-plotOutput
    output$eigenvalue <- renderPlot({
      EigenValuePlot(froot, input$espace)
    })

    #miamiPlot-plotOutput
    output$miami <- renderPlot({
      if(input$miamiPlot_integer <= input$espace) {
        miamiPlot(
          froot,
          input$miamiPlot_integer,
          Log1 = TRUE,
          Log2 = F,
          cex = input$miamiPlot_cex,
          pch = input$miamiPlot_pch,
          bty = "l",
          genomewideline = input$threshhold
        )
      }
    })

    #EigenGWASPlot-plotOutput
    output$eigengwas <- renderPlot({
      if (input$EigenGWASPlot_integer <= input$espace) {
        EigenGWASPlot(froot, input$EigenGWASPlot_integer)
      }
    })
  })
  
  #popup
  
  #GRM-   ImageSave-information
  observeEvent(input$GRM_save, {
    showModal(
      modalDialog(
        title = 'Save',
        textInput('GRM_filenameInput', label = 'Input Your Filename please', value = ''),
        checkboxGroupInput(
          'GRM_savetype',
          label = 'Choose your filetype',
          choices = list(
            '.pdf' = 1,
            '.jpeg' = 2,
            '.bmp' = 3,
            '.png' = 4
          ),
          selected = 4
        ),
        actionButton('GRM_savebutton', 'Save'),
        easyClose = TRUE
      )
    )
  })
  #GRM-   ImageSave-savebutton
  observeEvent(input$GRM_savebutton, {
    #.pdf file
    if ("1" %in% input$GRM_savetype) {
      GRM_SaveAsPdf(input$GRM_filenameInput)
    }
    
    #.jpeg file
    if ("2" %in% input$GRM_savetype) {
      GRM_SaveAsJpeg(input$GRM_filenameInput)
    }
    
    #.bmp file
    if ("3" %in% input$GRM_savetype) {
      GRM_SaveAsBmp(input$GRM_filenameInput)
    }
    
    #.png file
    if ("4" %in% input$GRM_savetype) {
      GRM_SaveAsPng(input$GRM_filenameInput)
    }
  })
  
  #Eigenvalue-   ImageSave-information
  observeEvent(input$Eigenvalue_save, {
    showModal(
      modalDialog(
        title = 'Save',
        textInput(
          'Eigenvalue_filenameInput',
          label = 'Input Your Filename please',
          value = ''
        ),
        checkboxGroupInput(
          'Eigenvalue_savetype',
          label = 'Choose your filetype',
          choices = list(
            '.pdf' = 1,
            '.jpeg' = 2,
            '.bmp' = 3,
            '.png' = 4
          ),
          selected = 4
        ),
        actionButton('Eigenvalue_savebutton', 'Save'),
        easyClose = TRUE
      )
    )
  })
  #Eigenvalue-   ImageSave-savebutton
  observeEvent(input$Eigenvalue_savebutton, {
    #.pdf file
    if ("1" %in% input$Eigenvalue_savetype) {
      Eigenvalue_SaveAsPdf(input$Eigenvalue_filenameInput)
    }
    
    #.jpeg file
    if ("2" %in% input$Eigenvalue_savetype) {
      Eigenvalue_SaveAsJpeg(input$Eigenvalue_filenameInput)
    }
    
    #.bmp file
    if ("3" %in% input$Eigenvalue_savetype) {
      Eigenvalue_SaveAsBmp(input$Eigenvalue_filenameInput)
    }
    
    #.png file
    if ("4" %in% input$Eigenvalue_savetype) {
      Eigenvalue_SaveAsPng(input$Eigenvalue_filenameInput)
    }
  })
  
  #miamiPlot-   ImageSave-information
  observeEvent(input$miamiPlot_save, {
    showModal(
      modalDialog(
        title = 'Save',
        textInput(
          'miamiPlot_filenameInput',
          label = 'Input Your Filename please',
          value = ''
        ),
        checkboxGroupInput(
          'miamiPlot_savetype',
          label = 'Choose your filetype',
          choices = list(
            '.pdf' = 1,
            '.jpeg' = 2,
            '.bmp' = 3,
            '.png' = 4
          ),
          selected = 4
        ),
        actionButton('miamiPlot_savebutton', 'Save'),
        easyClose = TRUE
      )
    )
  })
  #miamiPlot-   ImageSave-savebutton
  observeEvent(input$miamiPlot_savebutton, {
    #.pdf file
    if ("1" %in% input$miamiPlot_savetype) {
      miamiPlot_SaveAsPdf(
        input$miamiPlot_filenameInput,
        input$miamiPlot_integer,
        input$miamiPlot_cex,
        input$miamiPlot_pch,
        input$threshhold
      )
    }
    
    #.jpeg file
    if ("2" %in% input$miamiPlot_savetype) {
      miamiPlot_SaveAsJpeg(
        input$miamiPlot_filenameInput,
        input$miamiPlot_integer,
        input$miamiPlot_cex,
        input$miamiPlot_pch,
        input$threshhold
      )
    }
    
    #.bmp file
    if ("3" %in% input$miamiPlot_savetype) {
      miamiPlot_SaveAsBmp(
        input$miamiPlot_filenameInput,
        input$miamiPlot_integer,
        input$miamiPlot_cex,
        input$miamiPlot_pch,
        input$threshhold
      )
    }
    
    #.png file
    if ("4" %in% input$miamiPlot_savetype) {
      miamiPlot_SaveAsPng(
        input$miamiPlot_filenameInput,
        input$miamiPlot_integer,
        input$miamiPlot_cex,
        input$miamiPlot_pch,
        input$threshhold
      )
    }
  })
  
  
  #EigenGWASPlot-   ImageSave-information
  observeEvent(input$EigenGWASPlot_save, {
    showModal(
      modalDialog(
        title = 'Save',
        textInput(
          'EigenGWASPlot_filenameInput',
          label = 'Input Your Filename please',
          value = ''
        ),
        checkboxGroupInput(
          'EigenGWASPlot_savetype',
          label = 'Choose your filetype',
          choices = list(
            '.pdf' = 1,
            '.jpeg' = 2,
            '.bmp' = 3,
            '.png' = 4
          ),
          selected = 4
        ),
        actionButton('EigenGWASPlot_savebutton', 'Save'),
        easyClose = TRUE
      )
    )
  })
  #EigenGWASPlot-   ImageSave-savebutton
  observeEvent(input$EigenGWASPlot_savebutton, {
    #.pdf file
    if ("1" %in% input$EigenGWASPlot_savetype) {
      EigenGWASPlot_SaveAsPdf(input$EigenGWASPlot_filenameInput,
                              input$EigenGWASPlot_integer)
    }
    
    #.jpeg file
    if ("2" %in% input$EigenGWASPlot_savetype) {
      EigenGWASPlot_SaveAsJpeg(input$EigenGWASPlot_filenameInput,
                               input$EigenGWASPlot_integer)
    }
    
    #.bmp file
    if ("3" %in% input$EigenGWASPlot_savetype) {
      EigenGWASPlot_SaveAsBmp(input$EigenGWASPlot_filenameInput,
                              input$EigenGWASPlot_integer)
    }
    
    #.png file
    if ("4" %in% input$EigenGWASPlot_savetype) {
      EigenGWASPlot_SaveAsPng(input$EigenGWASPlot_filenameInput,
                              input$EigenGWASPlot_integer)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
