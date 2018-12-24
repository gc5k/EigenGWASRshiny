# Source helpers ----
source('EigenGWAS_Friends.R')
FN="arab";
PC=5;
inbred=F;
#Function below controls dataLoading
#RunEigenGWAS(FN, PC, inbred, "gear.jar");

RunPlink <- function(dat, PC, inbred=T, plink2) {
  layout(matrix(1:6, 2, 3))
  #make-grm
  grmCmd=paste(plink2, "--bfile ", dat, "--make-grm-gz --out ", dat)
  system(grmCmd)
  
  gz=gzfile(paste0(dat, ".grm.gz"))
  grm=read.table(gz, as.is = T)
  Ne=-1/mean(grm[grm[,1]!=grm[,2], 4]/2)
  Me=1/var(grm[grm[,1]!=grm[,2], 4]/2)
  print(paste("Ne=", format(Ne, digits = 2), "Me=", format(Me, digits = 2)))
  ## [1] "Ne= 293 Me= 395"
  hist(grm[grm[,1]!=grm[,2],4]/2, main="Pairwise relatedness", xlab="Relatedness score", breaks = 50)
  
  #pca
  pcaCmd=paste(plink2, "--bfile ", dat, "--pca 10 --out ", dat)
  system(pcaCmd)
  barplot(main="Top 10 eigenvalue", read.table(paste0(dat, ".eigenval"), as.is = T)[,1]/2, border = F)
  abline(h=1, col="red", lty=2, lwd=3)
  
  pc=read.table(paste0(dat, ".eigenvec"), as.is = T)
  plot(pc[,3], pc[,4], xlab="Eigenvector 1", ylab="Eigenvector 2", bty="n", main="Eigenspace", bty="n", col=ifelse(pc[,3]>0, "red", "blue"), pch=16, cex=0.5)
  
  #EigenGWAS
  for(i in 1:PC) {
    outRoot=paste0(dat, ".", i)
    liCmd=paste0(plink2, " --linear --bfile ", dat, " --pheno ", dat, ".eigenvec --out ", outRoot)
    system(liCmd)
  }

  # #plot
  # EigenRes=read.table(paste0(dat, ".assoc.linear"), as.is = T, header = T)
  # EigenRes$Praw=EigenRes$P
  # gc=qchisq(median(EigenRes$P), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
  # print(paste("GC = ", format(gc, digits = 4)))
  # ## [1] "GC =  9.047"
  # EigenRes$P=pchisq(qchisq(EigenRes$Praw, 1, lower.tail = F)/gc, 1, lower.tail = F)
  # manhattan(EigenRes, title="EigenGWAS 1", pch=16, cex=0.3, bty='n')
  # 
  # #QQplot
  # chiseq=rchisq(nrow(EigenRes), 1)
  # qqplot(chiseq, qchisq(EigenRes$Praw, 1, lower.tail = F), xlab=expression(paste("Theoretical ", chi[1]^2)), ylab=expression(paste("Observed ", chi[1]^2)), bty="n", col="grey", pch=16, cex=0.5)
  # points(sort(chiseq), sort(qchisq(EigenRes$P, 1, lower.tail = F)), col="black", pch=16, cex=0.5)
  # legend("topleft", legend = c("Raw", "GC correction"), pch=16, cex=0.5, col=c("grey", "black"), bty='n')
  # abline(a=0, b=1, col="red", lty=2)
  
  
}
#FileSave functions

manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=NULL, suggestiveline=-log10(1e-5), genomewideline=NULL, title="", annotate=NULL, ...) {

  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  if (!is.null(limitchromosomes)) {
    d=d[d$CHR %in% limitchromosomes, ]
  }
  
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0 #  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  colors <- rep(colors,max(d$CHR))[1:length(unique(d$CHR))]
  
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    Uchr=unique(d$CHR)
    for (i in 1:length(Uchr)) {
      if (i==1) {
        d[d$CHR==Uchr[i], ]$pos=d[d$CHR==Uchr[i], ]$BP
      } else {
        lastbase=lastbase+tail(subset(d, CHR==Uchr[i-1])$BP, 1)
        d[d$CHR==Uchr[i], ]$pos=d[d$CHR==Uchr[i], ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==Uchr[i], ]$pos[floor(length(d[d$CHR==Uchr[i], ]$pos)/2)+1])
    }
  }
  if (numchroms==1) {
    with(d, plot(main=title, pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  } else {
    with(d, plot(main=title, pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    Uchr=unique(d$CHR)
    for (i in 1:length(Uchr)) {
      with(d[d$CHR==Uchr[i], ], points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="green3", ...))
  }
  #  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (!is.null(genomewideline)) {
    abline(h=genomewideline, col="gray")
  } else {
    abline(h=-log10(0.05/nrow(d)), col="gray")    
  }
}

#GRM
GRM_SaveAsPdf <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.pdf')
  pdf(imagefile)
  layout(matrix(1:2, 1, 2))
  grmStats(FN)
  dev.off()
  list(src=imagefile)
}

GRM_SaveAsJpeg <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.jpeg')
  jpeg(imagefile)
  layout(matrix(1:2, 1, 2))
  grmStats(FN)
  dev.off()
  list(src=imagefile)
}

GRM_SaveAsBmp <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.bmp')
  bmp(imagefile)
  layout(matrix(1:2, 1, 2))
  grmStats(FN)
  dev.off()
  list(src=imagefile)
}

GRM_SaveAsPng <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.png')
  png(imagefile)
  layout(matrix(1:2, 1, 2))
  grmStats(FN)
  dev.off()
  list(src=imagefile)
}

#Eigenvalue
Eigenvalue_SaveAsPdf <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.pdf')
  pdf(imagefile)
  EigenValuePlot(FN, PC)
  dev.off()
  list(src=imagefile)
}

Eigenvalue_SaveAsJpeg <- function(path){
  imagefile <- paste('/Users/Shared/',path,'.jpeg')
  jpeg(imagefile)
  EigenValuePlot(FN, PC)
  dev.off()
  list(src=imagefile)
}

Eigenvalue_SaveAsBmp <- function(path,arg1,arg2,arg3){
  imagefile <- paste('/Users/Shared/',path,'.bmp')
  bmp(imagefile)
  EigenValuePlot(FN, PC)
  dev.off()
  list(src=imagefile)
}

Eigenvalue_SaveAsPng <- function(path,arg1,arg2,arg3){
  imagefile <- paste('/Users/Shared/',path,'.png')
  png(imagefile)
  EigenValuePlot(FN, PC)
  dev.off()
  list(src=imagefile)
}

#miamiPlot
miamiPlot_SaveAsPdf <- function(path,arg1,arg2,arg3,arg4){
  imagefile <- paste('/Users/Shared/',path,'.pdf')
  pdf(imagefile)
  miamiPlot(FN, arg1 , Log1 = TRUE, Log2 = F, cex=arg2, pch=arg3, bty="l",genomewideline = arg4)
  dev.off()
  list(src=imagefile)
}

miamiPlot_SaveAsJpeg <- function(path,arg1,arg2,arg3,arg4){
  imagefile <- paste('/Users/Shared/',path,'.jpeg')
  jpeg(imagefile)
  miamiPlot(FN, arg1 , Log1 = TRUE, Log2 = F, cex=arg2, pch=arg3, bty="l",genomewideline = arg4)
  dev.off()
  list(src=imagefile)
}

miamiPlot_SaveAsBmp <- function(path,arg1,arg2,arg3,arg4){
  imagefile <- paste('/Users/Shared/',path,'.bmp')
  bmp(imagefile)
  miamiPlot(FN, arg1 , Log1 = TRUE, Log2 = F, cex=arg2, pch=arg3, bty="l",genomewideline = arg4)
  dev.off()
  list(src=imagefile)
}

miamiPlot_SaveAsPng <- function(path,arg1,arg2,arg3,arg4){
  imagefile <- paste('/Users/Shared/',path,'.png')
  png(imagefile)
  miamiPlot(FN, arg1 , Log1 = TRUE, Log2 = F, cex=arg2, pch=arg3, bty="l",genomewideline = arg4)
  dev.off()
  list(src=imagefile)
}

#EigenGWASPlot
EigenGWASPlot_SaveAsPdf <- function(path,arg){
  imagefile <- paste('/Users/Shared/',path,'.pdf')
  pdf(imagefile)
  EigenGWASPlot(FN, arg)
  dev.off()
  list(src=imagefile)
}

EigenGWASPlot_SaveAsJpeg <- function(path,arg){
  imagefile <- paste('/Users/Shared/',path,'.jpeg')
  jpeg(imagefile)
  EigenGWASPlot(FN, arg)
  dev.off()
  list(src=imagefile)
}

EigenGWASPlot_SaveAsBmp <- function(path,arg){
  imagefile <- paste('/Users/Shared/',path,'.bmp')
  bmp(imagefile)
  EigenGWASPlot(FN, arg)
  dev.off()
  list(src=imagefile)
}

EigenGWASPlot_SaveAsPng <- function(path,arg){
  imagefile <- paste('/Users/Shared/',path,'.png')
  png(imagefile)
  EigenGWASPlot(FN, arg)
  dev.off()
  list(src=imagefile)
}