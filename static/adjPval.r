##adjust p values for multiple comparisons
mainpath=getwd()
filepath=paste(getwd(),'/static/OUT',sep='')
setwd(filepath)
mainpath=filepath
csvfiles <- list.files(path=mainpath, pattern="output.csv")
for (singlecsv in csvfiles) {
	fullfilename <-paste(mainpath,"/",singlecsv,sep="")
	generesult <- read.csv(file=fullfilename, header=FALSE, sep=",")
	pvalues <- generesult$V4
	fdrs <- p.adjust(pvalues,method="BH",n=length(pvalues))
	generesult$V7 <- formatC(fdrs,digits=5)
	oldfilename <- gsub("-output.csv","",fullfilename)
	newfilename <- paste(oldfilename,"-adj.csv",sep="")
	write.table(generesult,"adj.csv",append=T, col.names = T, row.names = F, sep=",")
}
