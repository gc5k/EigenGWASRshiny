# Source helpers ----
source('EigenGWAS_Friends.R')
FN="arab";
PC=5;
inbred=F;
#Function below controls dataLoading
#RunEigenGWAS(FN, PC, inbred, "gear.jar");


#FileSave functions

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