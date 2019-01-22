# Load packages ----
conf=read.table("EigenGWAS.conf", as.is = T)
library(shiny)
#library(shinyjs)
# Source helpers ----
source('helper.R')

plink2 = "./plink_mac"
if(length(grep("linux",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("linux")
  plink2 = "./plink_linux"
} else if(length(grep("apple",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("apple")
  plink2 = "./plink_mac"

#  system("git rev-list head --max-count 1 > gitTag.txt")
} else {
  print("windows")
  plink2 = "./plink.exe"
}

options(shiny.maxRequestSize=conf[1,2]*1024^2, shiny.launch.browser=T)
gTag=read.table("gitTag.txt")

# Define UI for EigenGWAS Application
ui <- fluidPage(
  #Title
  titlePanel("EigenGWAS"),
  h6(paste("Git version:", gTag[1,1])),
  hr(),
  #Gobal parameters
  fluidRow(
    column(6,
      fileInput('file_input', 
        paste0('3 source files (.bim, .bed, .fam) [< ', conf[1,2],' MB]'), 
        multiple = TRUE,
        accept = c("bed", "fam", "bim")
      )
    ),
    column(6,
           numericInput('autosome', "Autosome number", value=22)
           )
    ),

  hr(),
  fluidRow(
    column(2,
      radioButtons('bred',
        'Population type',
        choices = list('Outbred' = 'outbred', 'Inbred' = 'inbred'),
          selected = 'outbred',
          inline = T
        )
    ),
    column(4,
      sliderInput('maf_cut',
        'MAF threshold',
        value = 0, min = 0, max = 0.4, step = 0.05
      )
    ),
    column(6,
      sliderInput('espace',
        'Eigen space',
        value = 2, min = 1, max = 5, step = 1
      )
    )
  ),

  fluidRow(
    column(4,
      actionButton('run', 
        'EigenGWAS, Go!')
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
                         
          tabPanel('MAF',
            plotOutput('freq')
          ),

          tabPanel("PCA",
            sidebarPanel(
              sliderInput("x", "xLab", 1, 8, 1, step = 1),
              sliderInput("y", "yLab", 1, 8, 2, step = 1)
            ),

            mainPanel(
              plotOutput("PCA")
              )
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
                    "eSpace", min = 1, max = 5, 
                    value = c(1,1), dragRange=TRUE, step = 1),
                  selectInput('threshold',
                    'p-value cutoff',
                    choices = c(0.1, 0.05, 0.01, 0.005, 0.001), 
                    selected = 0.05
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
  currentFile <- reactive( {
    print("Reading files...")
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

    if(as.numeric(input$maf_cut) == 0) {
      frootMAF=froot
    } else {
      frootMAF=paste0(froot, "_t") 
      frqMAF=paste(plink2, "--allow-no-sex --bfile ", froot, " --autosome-num ", input$autosome, "--maf ", input$maf_cut, " --make-bed --out", frootMAF)
      system(frqMAF)
    }

    return (frootMAF)
  })

  observeEvent(input$run, {
    ##plink
    withProgress(message="EigenGWAS:", value=0, {
      PC = input$espace
      froot = currentFile()

      n = 5+PC
      incProgress(1/n, detail = paste0(" estimating freq ..."))
      frqCmd=paste(plink2, "--bfile ", froot, " --autosome-num ", input$autosome, "--freq --out", froot)
      system(frqCmd)

      incProgress(3/n, detail = paste0(" making grm ..."))
      grmCmd=paste(plink2, "--allow-no-sex --bfile", froot, " --autosome-num ", input$autosome, "--make-grm-gz --out", froot)
      system(grmCmd)

      gz=gzfile(paste0(froot, ".grm.gz"))
      grm=read.table(gz, as.is = T)

      incProgress(1/n, detail = paste0(" conducting PCA ..."))
      pcRun=input$space
      if(input$espace < 20) {
        pcRun = 20
      }
      pcaCmd=paste(plink2, "--allow-no-sex --bfile", froot, " --autosome-num ", input$autosome, "--pca", pcRun, "--out", froot)
      system(pcaCmd)

      #EigenGWAS
      for(i in 1:PC) {
        incProgress(1/n, detail = paste0(" scanning eSpace ", i))
        outRoot=paste0(froot, ".", i)
        liCmd=paste0(plink2, " --allow-no-sex --linear --bfile ", froot, " --autosome-num ", input$autosome, " --pheno ", froot, ".eigenvec --mpheno ", i," --out ", outRoot)
        system(liCmd)
      }
      incProgress(1/n, detail = paste0(" finishing EigenGWAS."))
    })

    sc=ifelse(input$bred == 'inbred', 2, 1)
    withProgress(message="Visualizing:", value=0, {
      PC = input$espace
      n=2+2*PC
      
      incProgress(1/n, detail = paste0(" MAF plot ... "))
      output$freq <- renderPlot({
        froot = currentFile()
        fq=read.table(paste0(froot, ".frq"), as.is = T, header = T)
        hist(fq$MAF, main="Minor allele frequency", xlab="MAF", xlim=c(0, 0.5), breaks = 50)
      })
      
      incProgress(1/n, detail = paste0(" MAF plot ... "))
      output$PCA <- renderPlot({
        froot = currentFile()
        layout(matrix(1:2, 1, 2))
        evalF=read.table(paste0(froot, ".eigenval"), as.is = T)
        barplot(evalF[,1]/sc, border = F, main="Eigenvalue")
        abline(h=1, lty=2, col="black")
        pcF=read.table(paste0(froot, ".eigenvec"), as.is = T)
        
        plot(main=paste0("eSpace ", input$x, " vs ", input$y), pcF[,input$x+2], pcF[,input$y+2], xlab=paste0("eSpace ", input$x), ylab=paste0("eSpace ", input$y), bty='n', pch=16, cex=0.5,
             col=ifelse(pcF[,3]<0, "red", "blue"))
        abline(v=0, h=0, col="grey", lty=2)
      })

      incProgress(1/n, detail = paste0(" GRM plot ... "))
      output$grm <- renderPlot({
        froot = currentFile()

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
        froot = currentFile()
        n=2+2*PC
        Evev=read.table(paste0(froot, ".eigenval"), as.is = T)
        GC=array(0, dim=PC)
        for(i in 1:PC) {
          eg = read.table(paste0(froot, ".", i, ".assoc.linear"), as.is = T, header = T)
          GC[i] = qchisq(median(eg$P, na.rm = T), 1, lower.tail = F)/qchisq(0.5, 1)
        }

        egc=matrix(c(Evev[1:PC,1]/sc, GC), PC, 2, byrow = F)
        rownames(egc)=seq(1, PC)
        barplot(t(egc), beside = T, border = F, xlab="eSpace", ylim=c(0,max(egc)+2))
        abline(h=1, lty=2, lwd=2)
        legend("topright", legend = c("Eigenvalue", expression(paste(lambda[gc]))), pch=15, col=c("black", "grey"), bty='n')
      })
      
      incProgress(1/n, detail = paste0(" EigenGWAS visualization ... "))
      output$eigengwas <- renderPlot({
        froot = currentFile()
        pcIdx=input$EigenGWASPlot_espace[1]
        print(paste0("ePC ", pcIdx))
        if (pcIdx > input$espace) {
          showNotification(paste0("eSpace ", pcIdx, " hasn't been scanned yet."), type = "warning")
          return()
        }

        egResF=paste0(froot, ".",pcIdx, ".assoc.linear")
        layout(matrix(1:2, 1, 2))
        EigenRes=read.table(egResF, as.is = T, header = T)
        EigenRes=EigenRes[which(!is.na(EigenRes$P)),]
        EigenRes$Praw=EigenRes$P
        gc=qchisq(median(EigenRes$P, na.rm = T), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
        print(paste0("GC = ", format(gc, digits = 4)))
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

  output$eReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    
    filename = "egReport.html",
    content = function(file) {
      
      isRunable = TRUE
      if ( !is.null(input$file_input$datapath)) {
        frt=currentFile()
#        frt = substr(input$file_input$datapath[1], 1, nchar(input$file_input$datapath[1])-4)
        print(frt)
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
      params <- list(froot = frt,
                     espace = input$espace, 
                     sc = ifelse(input$bred == 'inbred', 2, 1),
                     pcut = as.numeric(input$threshold))
      
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
