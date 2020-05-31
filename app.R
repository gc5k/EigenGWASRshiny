# Load packages ----
conf=read.table("EigenGWAS.conf", as.is = T)

library(shiny)
#library(shinyjs)
# Source helpers ----
source('helper.R')

unzip = "unzip -o -q plink.zip"
plink2 = "./plink_mac"
if(length(grep("linux",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("linux")
  system(paste0(unzip," plink_linux"))
  plink2 = "./plink_linux"
} else if(length(grep("apple",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("apple")
  system(paste0(unzip," plink_mac"))
  plink2 = "./plink_mac"
  
  #  system("git rev-list head --max-count 1 > gitTag.txt")
} else {
  print("windows")
  system("expand plink.cab plink_win.exe")
  plink2 = "plink_win.exe"
}

options(shiny.maxRequestSize=conf[1,2]*1024^2, shiny.launch.browser=T)
gTag=read.table("gitTag.txt")
multipleSpace = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
# Define UI for EigenGWAS Application
ui <- fluidPage(
  theme = "style.css",
  div(style = "padding: 1px 0px; width: '100%'",
      titlePanel(
        title = "",
        windowTitle = "EigenGWAS"
      )
  ),
  navbarPage(
    title = div(
      span(
        #img(src = "logo.png"),
        strong("EigenGWAS"),
        style = "position: relative; top: 50%; transform: translateY(-50%);"
      )
    ),
    id = "inNavbar",
    tabPanel(
      title = "Data Input",
      value = "datainput",
      fluidRow(
        column(
          6,
          fileInput(
            "file_input",
            paste0('3 source files (.bim, .bed, .fam) [< ', conf[1,2],' MB]'),
            multiple = TRUE,
            accept = c("bed", "fam", "bim")
          )
        ),
        column(
          6,
          numericInput('autosome', "Autosome number", value=22)
        )
      ),
      hr(),
      fluidRow(
        column(
          2,
          radioButtons(
            'bred',
            'Population type',
            choices = list('Outbred' = 'outbred', 'Inbred' = 'inbred'),
            selected = 'outbred',
            inline = T
          )
        ),
        column(
          4,
          sliderInput(
            'maf_cut',
            'MAF threshold',
            value = 0.01, min = 0, max = 0.4, step = 0.05
          )
        ),
        column(
          4,
          sliderInput(
            'espace',
            'Eigen space',
            value = 2, min = 1, max = 5, step = 1
          )
        )
      ),
      fluidRow(
        column(
          4,
          actionButton(
            'run',
            'EigenGWAS, Go!'
          )
        )
      )
    ),
    tabPanel(
      title = "Visualization",
      value = "visualization",
      fluidRow(
        column(12, 
               mainPanel(
                 tabsetPanel(
                   type = 'tabs',
                   tabPanel(
                     'MAF',
                     plotOutput('freq')
                   ),
                   tabPanel(
                     'GRM',
                     plotOutput('grm')
                   ),
                   tabPanel(
                     "PCA",
                     sidebarPanel(
                       sliderInput("x", "xLab", 1, 8, 1, step = 1),
                       sliderInput("y", "yLab", 1, 8, 2, step = 1)
                     ),
                     mainPanel(
                       plotOutput("PCA")
                     )
                   ),          					
                   tabPanel(
                     'Eigenvalue',
                     plotOutput('eigenvalue')
                   ),
                   tabPanel(
                     'EigenGWAS',
                     tabPanel(
                       'EigenGWAS visualization',
                       sidebarPanel(
                         sliderInput(
                           'EigenGWASPlot_espace',
                           "eSpace", 
                           min = 1, max = 5, value = c(1,1), dragRange=TRUE, step = 1),
                         selectInput(
                           'threshold',
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
      ),
      hr(),
      fluidRow(
        column(
          12,
          mainPanel(
            column(4,
                   downloadButton(
                     'eReport', 
                     'Generate eReport'
                   )
            ),
            column(
              6,
              DT::dataTableOutput('topsnp')
            )
          )
        )
      )
    ),
    tabPanel(
      title = "Power",
      value = "power",
      fluidPage(
        # Application title
        titlePanel("EigenGWAS Power Calculator"),
        hr(),
        # Sidebar with a slider input for number of bins
        sidebarLayout(
          sidebarPanel(
            #        numericInput("n", "Sample size",
            #                     value = 500, min=100),
            numericInput("m", "Marker number",
                         value = 10000, min=1),
            sliderInput("w1",
                        "Subpop 1 proportion",
                        min = 0.05, max = 0.95, value = 0.5),
            sliderInput("p1", "Frequency at population 1",
                        min = 0.01, max=0.99, value=0.3),
            sliderInput("p2", "Frequency at population 2",
                        min = 0.01, max=0.99, value=0.5),
            selectInput("alpha", "Alpha",
                        choices = c("0.001", "0.005", "0.01", "0.05", "0.1"), 
                        selected="0.05"),
            selectInput("beta", "Beta",
                        choices = c("0.5", "0.6", "0.7", "0.8", "0.9"),
                        selected="0.8"),
            actionButton(
              'powerupdate',
              'Update power',
              icon = icon("refresh")
            )
            #submitButton(text="Update power", icon("refresh"))
          ),
          
          # Show a plot of the generated distribution
          mainPanel(
            plotOutput("distPlot")
          )
        )
      )
    ),
    tabPanel(
      title = "About",
      value = "about",
      tags$br(),
      tags$h3("Citation"),
      tags$p(HTML("<a href=\"https://www.nature.com/articles/hdy201625\">Chen, G.B. et al, EigenGWAS: finding loci under selection through genome-wide association studies of eigenvectors in structured populations, Heredity, 2016, 117:51-61.</a>")),
      tags$br(),
      tags$p(HTML(paste("Git version:", gTag[1,1])))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
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
    updateTabsetPanel(session, "inNavbar", selected = "visualization")
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
      # 20200520 the next two lines of code for data reading is duplicated
      #gz=gzfile(paste0(froot, ".grm.gz"))
      #grm=read.table(gz, as.is = T)
      
      incProgress(1/n, detail = paste0(" conducting PCA ..."))
      #edit 20200520
      #original code, but I cannot find 'space' in input 
      #pcRun=input$space 
      #espace is not likely reach 20, so the following if() is quite unreasonable 
      #if(input$espace < 20) {
      #pcRun = 20
      #}
      pcRun=10 # 20200520: directly set pcRun to 10
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
      
      incProgress(1/n, detail = paste0(" PCA plot ... "))
      output$PCA <- renderPlot({
        froot = currentFile()
        layout(matrix(1:2, 1, 2))
        evalF=read.table(paste0(froot, ".eigenval"), as.is = T)
        pcF=read.table(paste0(froot, ".eigenvec"), as.is = T)
        
        barplot(evalF[,1]/sc, border = F, main="Eigenvalue")
        abline(h=1, lty=2, col="black")
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

        off_diagnol = grm[grm[,1]!=grm[,2], 4]
        Ne=-1/mean(off_diagnol/sc)
        Me=1/var(off_diagnol/sc)

        hist(off_diagnol/sc, main="Pairwise relatedness", xlab="Relatedness score", breaks = 50)
        
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
          GC[i] = qchisq(median(eg$P, na.rm = T), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
        }
        
        egc=matrix(c(Evev[1:PC,1]/sc, GC), PC, 2, byrow = F)
        #egc=matrix(c(Evev[1:PC,1], GC), PC, 2, byrow = F)
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
        EigenRes$P=pchisq((EigenRes$STAT)^2/gc, 1, lower.tail = F)
        manhattan(EigenRes, genomewideline = -log10(as.numeric(input$threshold)/nrow(EigenRes)), title=paste("eSpace ", pcIdx), annotatePvalue = TRUE, pch=16, cex=0.3, bty='n')
        
        #QQplot
        chiseq=qchisq(seq(1/nrow(EigenRes), 1-1/nrow(EigenRes), length.out = nrow(EigenRes)), 1)
        qqplot(chiseq, qchisq(EigenRes$Praw, 1, lower.tail = F), xlab=expression(paste("Theoretical ", chi[1]^2)), ylab=expression(paste("Observed ", chi[1]^2)), bty="n", col="grey", pch=16, cex=0.5)
        points(sort(chiseq), sort(qchisq(EigenRes$P, 1, lower.tail = F)), col="black", pch=16, cex=0.5)
        legend("topleft", legend = c("Raw", "GC correction"), pch=16, cex=0.5, col=c("grey", "black"), bty='n')
        abline(a=0, b=1, col="red", lty=2)
        
        output$topsnp <- DT::renderDataTable(
          DT::datatable(
            EigenRes[order(EigenRes$P),][c(1:5),c(1:3,9,10)],
            caption = htmltools::tags$caption(
              style = 'caption-side:bottom;text-align:center;','Table: ',htmltools::em(paste0('Top hits in EigenGWAS, eSpace',pcIdx))),
            options = list(pageLength = 5, 
                           searching = FALSE, 
                           dom = ''
            ),
            rownames = FALSE
          )
        )
      })
    })
  })
  
  observeEvent(input$powerupdate, {
    output$distPlot <- renderPlot({
      m=input$m
      alpha=as.numeric(input$alpha)
      pcut=alpha/m
      chiT=qchisq(pcut, 1, lower.tail = F)
      
      n=c(100,  200,  500,   1000,  1500,  2000,
          5000, 7500, 10000, 15000, 20000, 50000)
      PW=matrix(0, 2, length(n))
      
      w1=input$w1
      w2=1-w1
      
      p1=input$p1
      h1=2*p1*(1-p1)
      p2=input$p2
      h2=2*p2*(1-p2)
      
      p=w1*p1+w2*p2
      H=w1*h1+w2*h2
      
      for(i in 1:length(n)) {
        ncpA=4*n[i]*w1*w2*(p1-p2)^2/(2*p*(1-p))
        ncpD=n[i]*w1*w2*(h1-h2)^2/H
        
        PW[1,i]=pchisq(chiT, 1, ncp=ncpA, lower.tail = F)
        PW[2,i]=pchisq(chiT, 1, ncp=ncpD, lower.tail = F)
      }
      colnames(PW)=n
      par(las=2)
      barplot(PW, beside = T, border = F, ylab="Statistical power", xlab="Sample size")
      abline(h=as.numeric(input$beta), lty=2, col="grey")
      legend("topleft", legend=c("Add", "Dom"), pch=15, col=c("black", "grey"), bty='n')
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
  




