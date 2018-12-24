# Load packages ----
library(shiny)
#library(shinyjs)
# Source helpers ----
source('helper.R')
options(shiny.maxRequestSize=200*1024^2, shiny.launch.browser=T)

plink2 = "./plink_mac"
if(length(grep("linux",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("linux")
  plink2 = "./plink_linux"
} else if(length(grep("apple",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("apple")
  plink2 = "./plink_mac"
}

# Define UI for EigenGWAS Application
ui <- fluidPage(
  #Title
  titlePanel("EigenGWAS"),
  hr(),
  #Gobal parameters
  fluidRow(
    column(8,
      fileInput('file_input', 
        '3 source files (*.bim, *.bed, *.fam) [< 200 MB]', 
        multiple = TRUE,
        accept = c("bed", "fam", "bim")
      )
    )
  ),
  fluidRow(
    column(2,
      radioButtons('bred',
        'Population type',
        choices = list('Outbred' = 'outbred', 'Inbred' = 'inbred'),
        selected = 'outbred',
        inline = T
      )
    ),
    column(2,
      numericInput('espace',
        'Eigen space',
        value = 2,
        min = 1,
        max = 10,
        step = 1
      )
    ),
    column(3,
      selectInput('threshold',
        'p-value cutoff',
        choices = c(0.1, 0.05, 0.01, 0.005, 0.001), 
        selected = 0.05
      )
    )
  ),

  fluidRow(
    column(4,
      actionButton('run', 
        'EigenGWAS, Go!'
      )
    )
  ),

  hr(),
  #TabPanels
  fluidRow(
    column(3,
           downloadButton('eReport', 
                          'Generate eReport'
           )
    )
  ),
  fluidRow(
    column(12, 
      mainPanel(
        tabsetPanel(type = 'tabs',

          tabPanel('Freq',
            plotOutput('freq')
          ),
          tabPanel('PCA',
            plotOutput('PCA')
          ),
          tabPanel('GRM',
            plotOutput('grm')
          ),

          tabPanel('Eigenvalue',
            plotOutput('eigenvalue')
          ),

          tabPanel('EigenGWAS',
            tabPanel('EigenGWAS visualization',
              sidebarPanel(
                sliderInput('EigenGWASPlot_espace',
                  "PC:",
                  min = 1,
                  max = 5,
                  value = 1,
                  step = 1
                )
              ),
              mainPanel(plotOutput('eigengwas'))
            )
          )
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #Plot on the web

  observeEvent(input$run, {

    FileLoad=0
    str=""
    if(length(which(grepl("*.bed", input$file_input$name)))  != 1) {
      str=paste(str, "No bed file found.") 
    } else {
      FileLoad=FileLoad+1
    }
    
    if(length(which(grepl("*.bim", input$file_input$name)))  != 1) {
      str=paste(str, "\nNo bim file found.")
    } else {
      FileLoad=FileLoad+1
    }

    if(length(which(grepl("*.fam", input$file_input$name)))  != 1) {
      str=paste(str, "\nNo fam file found.")
    } else {
      FileLoad=FileLoad+1
    }
    
    if (FileLoad < 3) {
      showNotification(str, duration = 5, type = "error")
      return()
    } else if (FileLoad > 3) {
      showNotification("More than 3 files selected", duration = 5, type="error")
    }

    idx=grep(".bed$", input$file_input$datapath)
    if (length(idx)==1) {
      rt=substr(input$file_input$datapath[idx], 1, nchar(input$file_input$datapath[idx])-4)
    }
    for (i in 1:3) {
      if (i != idx) {
        f1 = input$file_input$datapath[i]
        tl = substr(f1, nchar(f1)-2, nchar(f1))
        file.symlink(f1, paste0(rt, ".", tl))
      }
    }

    froot = substr(input$file_input$datapath[idx], 1, nchar(input$file_input$datapath[idx])-4)

    PC = input$espace
    ##plink
    withProgress(message="EigenGWAS:", value=0, {
      n = 5+PC
      incProgress(1/n, detail = paste0(" estimating freq ..."))
      frqCmd=paste(plink2, "--bfile ", froot, "--freq --out", froot)
      system(frqCmd)

      incProgress(3/n, detail = paste0(" making grm ..."))
      grmCmd=paste(plink2, "--bfile", froot, "--make-grm-gz --out", froot)
      system(grmCmd)

      gz=gzfile(paste0(froot, ".grm.gz"))
      grm=read.table(gz, as.is = T)

      incProgress(1/n, detail = paste0(" conducting PCA ..."))
      pcRun=input$space
      if(input$espace < 20) {
        pcRun = 20
      }
      pcaCmd=paste(plink2, "--bfile", froot, "--pca", pcRun, "--out", froot)
      system(pcaCmd)

    #EigenGWAS
      for(i in 1:PC) {
        incProgress(1/n, detail = paste0(" scanning eSpace ", i))
        outRoot=paste0(froot, ".", i)
        liCmd=paste0(plink2, " --linear --bfile ", froot, " --pheno ", froot, ".eigenvec --mpheno ", i," --out ", outRoot)
        system(liCmd)
      }
      incProgress(1/n, detail = paste0(" finishing EigenGWAS."))
    })

    sc=ifelse(input$bred == 'inbred', 2, 1)
    withProgress(message="Visualizing:", value=0, {
      n=2+2*PC

      incProgress(1/n, detail = paste0(" MAF plot ... "))
      output$freq <- renderPlot({
        fq=read.table(paste0(froot, ".frq"), as.is = T, header = T)
        hist(fq$MAF, main="Minor allele frequency", xlab="MAF", xlim=c(0, 0.5), breaks = 50)
      })

      incProgress(1/n, detail = paste0(" MAF plot ... "))
      output$PCA <- renderPlot({
        layout(matrix(1:2, 1, 2))
        evalF=read.table(paste0(froot, ".eigenval"), as.is = T)
        barplot(evalF[,1]/sc, border = F, main="Eigenvalue")
        abline(h=1, lty=2, col="black")
        pcF=read.table(paste0(froot, ".eigenvec"), as.is = T)
        plot(main="eSpace 1 vs 2", pcF[,3], pcF[,4], xlab="eSpace 1", ylab="eSpace 2", bty='n', pch=16, cex=0.5,
             col=ifelse(pcF[,3]<0, "red", "blue"))
      })

      incProgress(1/n, detail = paste0(" GRM plot ... "))
      output$grm <- renderPlot({
        layout(matrix(1:2, 1, 2))
        gz=gzfile(paste0(froot, ".grm.gz"))
        grm=read.table(gz, as.is = T)
        Ne=-1/mean(grm[grm[,1]!=grm[,2], 4]/sc)
        Me=1/var(grm[grm[,1]!=grm[,2], 4]/sc)
        hist(grm[grm[,1]!=grm[,2],4]/sc, main="Pairwise relatedness ", xlab="Relatedness score", breaks = 50)
        
        nn=nrow(read.table(paste0(froot, ".fam"), as.is = T))
        mm=nrow(read.table(paste0(froot, ".bim"), as.is = T))
        legend("topright", legend = c(paste0("ne=", format(Ne, digits=3, nsmall=2), ' [',nn, ']'), paste0("me=", format(Me, digits=3, nsmall=2), ' [',mm,']')), bty='n')

        hist(grm[grm[,1]==grm[,2],4]/sc, main="Diagonal relatedness", xlab="Relatedness score", breaks = 15)
      })

      incProgress(1/n, detail = paste0(" Eigenvalue plot ... "))
      output$eigenvalue <- renderPlot( {
        Evev=read.table(paste0(froot, ".eigenval"), as.is = T)
        GC=array(0, dim=PC)
        for(i in 1:PC) {
          eg = read.table(paste0(froot, ".", i, ".assoc.linear"), as.is = T, header = T)
          GC[i] = qchisq(median(eg$P), 1, lower.tail = F)/qchisq(0.5, 1)
        }

        egc=matrix(c(Evev[1:PC,1]/sc, GC), PC, 2, byrow = F)
        rownames(egc)=seq(1, PC)
        barplot(t(egc), beside = T, border = F, xlab="eSpace", ylim=c(0,max(egc)+2))
        abline(h=1, lty=2, lwd=2)
        legend("topright", legend = c("Eigenvalue", expression(paste(lambda[gc]))), pch=15, col=c("black", "grey"), bty='n')
      })

      incProgress(1/n, detail = paste0(" EigenGWAS visualization ... "))
      output$eigengwas <- renderPlot({
        #plot
        pcIdx=input$EigenGWASPlot_espace
        layout(matrix(1:2, 1, 2))
        EigenRes=read.table(paste0(froot, ".",pcIdx, ".assoc.linear"), as.is = T, header = T)
        EigenRes$Praw=EigenRes$P
        gc=qchisq(median(EigenRes$P), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
        print(paste("GC = ", format(gc, digits = 4)))
        EigenRes$P=pchisq(qchisq(EigenRes$Praw, 1, lower.tail = F)/gc, 1, lower.tail = F)
        manhattan(EigenRes, genomewideline = -log10(as.numeric(input$threshold)/nrow(EigenRes)), title=paste("eSpace ", pcIdx), pch=16, cex=0.3, bty='n')

        #QQplot
        chiseq=qchisq(seq(1/nrow(EigenRes), 1-1/nrow(EigenRes), length.out = nrow(EigenRes)), 1)
        qqplot(chiseq, qchisq(EigenRes$Praw, 1, lower.tail = F), xlab=expression(paste("Theoretical ", chi[1]^2)), ylab=expression(paste("Observed ", chi[1]^2)), bty="n", col="grey", pch=16, cex=0.5)
        points(sort(chiseq), sort(qchisq(EigenRes$P, 1, lower.tail = F)), col="black", pch=16, cex=0.5)
        legend("topleft", legend = c("Raw", "GC correction"), pch=16, cex=0.5, col=c("grey", "black"), bty='n')
        abline(a=0, b=1, col="red", lty=2)
      })
    })
    
    
  })

  observeEvent(input$EigenGWASPlot_espace, {
    output$eigengwas <- renderPlot({
      pcIdx=input$EigenGWASPlot_espace
      if(pcIdx > input$espace) {
        showNotification(paste0("eSpace ", pcIdx, " hasn't been scanned."), duration = 5, type = "message")
        return()
      }

      froot = substr(input$file_input$datapath[1], 1, nchar(input$file_input$datapath[1])-4)

      layout(matrix(1:2, 1, 2))
      EigenRes=read.table(paste0(froot, ".",pcIdx, ".assoc.linear"), as.is = T, header = T)
      EigenRes$Praw=EigenRes$P
      gc=qchisq(median(EigenRes$P), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
      print(paste("GC = ", format(gc, digits = 4)))
      EigenRes$P=pchisq(qchisq(EigenRes$Praw, 1, lower.tail = F)/gc, 1, lower.tail = F)
      manhattan(EigenRes, genomewideline = -log10(as.numeric(input$threshold)/nrow(EigenRes)), title=paste("eSpace ", pcIdx), pch=16, cex=0.3, bty='n')
    
    #QQplot
      chiseq=qchisq(seq(1/nrow(EigenRes), 1-1/nrow(EigenRes), length.out = nrow(EigenRes)), 1)
      qqplot(chiseq, qchisq(EigenRes$Praw, 1, lower.tail = F), xlab=expression(paste("Theoretical ", chi[1]^2)), ylab=expression(paste("Observed ", chi[1]^2)), bty="n", col="grey", pch=16, cex=0.5)
      points(sort(chiseq), sort(qchisq(EigenRes$P, 1, lower.tail = F)), col="black", pch=16, cex=0.5)
      legend("topleft", legend = c("Raw", "GC correction"), pch=16, cex=0.5, col=c("grey", "black"), bty='n')
      abline(a=0, b=1, col="red", lty=2)
    })
  })
  
  output$eReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"

      filename = "egReport.html",
      content = function(file) {
        
        isRunable = TRUE
        if ( !is.null(input$file_input$datapath)) {
          frt = substr(input$file_input$datapath[1], 1, nchar(input$file_input$datapath[1])-4)
          if ( !file.exists(paste0(frt, ".eigenval"))) {
            isRunable = FALSE
          } else if (!file.exists(paste0(frt, ".frq"))) {
            isRunable = FALSE
          } else if (!file.exists(paste0(frt, ".eigenvec"))) {
            isRunable = FALSE
          } else {
            for(i in 1:input$espace) {
              if (!file.exists(paste0(frt, ".", i, ".assoc.linear"))) {
                isRunable = FALSE
              }
            }
          }
        } else {
          isRunable = FALSE
        }
        if(!isRunable) {
          showNotification("No results generated yet", duration = 5, type = "message")
          return()
        }
        
        
        tempReport <- file.path(tempdir(), "egReport.Rmd")
        file.copy("egReport.Rmd", tempReport, overwrite = TRUE)


      # Set up parameters to pass to Rmd document
        params <- list(froot = substr(input$file_input$datapath[1], 1, nchar(input$file_input$datapath[1])-4),
                     espace = input$espace, 
                     sc = ifelse(input$bred == 'inbred', 2, 1),
                     pcut = as.numeric(input$threshold))

        rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
        )
      }
  )
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
                              input$EigenGWASPlot_espace)
    }
    
    #.jpeg file
    if ("2" %in% input$EigenGWASPlot_savetype) {
      EigenGWASPlot_SaveAsJpeg(input$EigenGWASPlot_filenameInput,
                               input$EigenGWASPlot_espace)
    }
    
    #.bmp file
    if ("3" %in% input$EigenGWASPlot_savetype) {
      EigenGWASPlot_SaveAsBmp(input$EigenGWASPlot_filenameInput,
                              input$EigenGWASPlot_espace)
    }
    
    #.png file
    if ("4" %in% input$EigenGWASPlot_savetype) {
      EigenGWASPlot_SaveAsPng(input$EigenGWASPlot_filenameInput,
                              input$EigenGWASPlot_espace)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
