# R script for plotting SSE number histogram
# Alex Stivala
# August 2008
# $Id: plotssehistogram.r 1946 2008-10-04 05:14:43Z astivala $

# EPS suitable for inserting into LaTeX
postscript('ssehistogram.eps',onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)
ssenums <- read.table('ssenums.list')
hist(ssenums$V1,freq=TRUE,density=20,main='Histogram of tableau sizes',xlab='Tableau size (number of SSEs in domain)',breaks=100)
dev.off()
