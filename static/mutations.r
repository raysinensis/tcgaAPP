##Sig mutations files
library(FirebrowseR)
cohorts = Metadata.Cohorts(format = "csv")
myArgs2 <- commandArgs(trailingOnly = TRUE)
mainpath<-getwd()
outputfile<-paste(mainpath,"/static/OUT/final.csv",sep='')
#file.remove(outputfile)
genequery=myArgs2
mutations=Analyses.Mutation.SMG(format="tsv",gene=genequery)
sigmuts=mutations$cohort[mutations$p<=0.05]
if (length(sigmuts)==0){sigmuts=list("holder")}

##copy number changes
deletions=tryCatch(Analyses.CopyNumber.Genes.Deleted(format="tsv",gene=genequery),error=function(e) {errorin=list()})
if (length(deletions)==0){sigdels=list("holder")
} else {sigdels=deletions$cohort}

insertions=tryCatch(Analyses.CopyNumber.Genes.Amplified(format="tsv",gene=genequery),error=function(e) {errorin=list()})
if (length(insertions)==0){sigins=list("holder")
}else{sigins=insertions$cohort}

filepath=paste(getwd(),'/static/OUT',sep='')
csvfiles <-paste(filepath,'/adj.csv',sep='')
finalpath=paste(filepath,'/final.csv',sep='')

##methylation files
methpath=paste(mainpath,'/static/methylation.csv',sep='')
generM <- read.csv(file=methpath, header=TRUE, sep=",")
for (singlecsv in csvfiles) {
	generesult <- read.csv(file=singlecsv, header=FALSE, sep=",")
	generesult=as.data.frame(generesult)
	generesult[,8:10]<-NA
	for (i in 1:(length(sigmuts))){
		generesult[generesult[2]==sigmuts[i],8]<-"MUT"}
	for (i in 1:(length(sigins))){
		generesult[generesult[2]==sigins[i],9]<-"INS"}
	for (i in 1:(length(sigdels))){
		generesult[generesult[2]==sigdels[i],10]<-"DEL"}
	#for (i in 1:(nrow(generesult))){
		#cancer=as.character(generesult[i,2])
		#generesult[i,11]<-as.numeric(as.character(generM[[cancer]][generM$gene==genequery]))}
	generesult[is.na(generesult)]<-"-"
	write.table(generesult,finalpath,append=F, col.names = T, row.names = F, sep=",")
}
