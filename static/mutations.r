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
setwd(filepath)
mainpath=filepath
csvfiles <- list.files(path=mainpath, pattern="adj.csv")

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
	write.table(generesult,"final.csv",append=F, col.names = T, row.names = F, sep=",")
}
