##adjust p values for multiple comparisons
mainpath=getwd()
filepath=paste(getwd(),'/static/OUT',sep='')
setwd(filepath)
mainpath=filepath
csvfiles <- list.files(path=mainpath, pattern="output.csv")
cancerlist=c()
for (singlecsv in csvfiles) {
	fullfilename <-paste(mainpath,"/",singlecsv,sep="")
	generesult <- read.csv(file=fullfilename, header=FALSE, sep=",")
	pvalues <- generesult$V4
	fdrs <- p.adjust(pvalues,method="BH",n=length(pvalues))
	generesult$V7 <- formatC(fdrs,digits=5)
	oldfilename <- gsub("-output.csv","",fullfilename)
	newfilename <- paste(oldfilename,"-adj.csv",sep="")
	write.table(generesult,"adj.csv",append=F, col.names = F, row.names = F, sep=",")
	cancerlist<-(generesult$V2)
}

diff.Exp.Genes = as.character(generesult$V1[1])
df=data.frame(matrix(NA, nrow = length(cancerlist), ncol = 4))
i=1

for (cancer in cancerlist){
cancer.Type=as.character(cancer)
library(FirebrowseR)
cohorts = Metadata.Cohorts(format = "csv")

##looking for cancer type
cancertype.Pats = Samples.Clinical(cohort = cancer.Type, format="tsv")


##dim(cancertype.Pats)

##pulling all patients
all.Received = F
page.Counter = 1
page.size = 150
cancertype.Pats = list()
while(all.Received == F){
  cancertype.Pats[[page.Counter]] = Samples.Clinical(format = "csv",
                                               cohort = cancer.Type,
                                               page_size = page.size,
                                               page = page.Counter)
  if(page.Counter > 1)
    colnames(cancertype.Pats[[page.Counter]]) = colnames(cancertype.Pats[[page.Counter-1]])

  if(nrow(cancertype.Pats[[page.Counter]]) < page.size){
    all.Received = T
  } else{
    page.Counter = page.Counter + 1
  }
}

cancertype.Pats = do.call(rbind, cancertype.Pats)
##dim(cancertype.Pats)

##pulling gene expression info

all.Found = F
page.Counter = 1
mRNA.Exp = list()
page.Size = 2000 # using a bigger page size is faster
while(all.Found == F){
  mRNA.Exp[[page.Counter]] = Samples.mRNASeq(format = "tsv",
                                             gene = diff.Exp.Genes,
                                             cohort = cancer.Type,
                                             page_size = page.Size,
                                             page = page.Counter)
  if(nrow(mRNA.Exp[[page.Counter]]) < page.Size){
    all.Found = T
  } else {
    page.Counter = page.Counter + 1
	}
}
mRNA.Exp = do.call(rbind, mRNA.Exp)
dim(mRNA.Exp)

##normal.Tissue.Pats = which(mRNA.Exp$sample_type[1] == "N")
##cancer.Tissue.Pats = which(mRNA.Exp$sample_type[1] == "T")
##patient.Barcodes = mRNA.Exp$tcga_participant_barcode[normal.Tissue.Pats]
##mRNA.Exp = mRNA.Exp[which(mRNA.Exp$tcga_participant_barcode %in% patient.Barcodes &
##                            mRNA.Exp$sample_type %in% c("NT", "TP")), ]

normal_med=median(mRNA.Exp$expression_log2[substr(mRNA.Exp$sample_type,start=1,stop=1)=="N"])
cancer_med=median(mRNA.Exp$expression_log2[substr(mRNA.Exp$sample_type,start=1,stop=1)=="T"])
if (is.na(normal_med)){
df$"X1"[i]=0}
else{
df$"X1"[i]=normal_med}
df$"X2"[i]=cancer_med
df$"X3"[i]=generesult$"V5"[i]
df$"X4"[i]=cancer.Type
i=i+1
}

widthcal=300*length(cancerlist)
colnames(df)=c("normal_median","tumor_median","surv_cutoff","tumor_type")
write.table(df,"level.csv",append=F, col.names = T, row.names = F, sep=",")
library(ggplot2)
library(reshape)
df2 <- melt(df, id.vars = "tumor_type")
colnames(df2)=c("tumor_type","condition","RSEM_log2")
p<-ggplot(data=df2, aes(x=tumor_type, y=RSEM_log2, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=sprintf("%0.2f",RSEM_log2)), vjust=1.5, color="black",
            position = position_dodge(1), size=2.5)+
  scale_fill_brewer(palette="Paired")+
  theme_classic(base_size = 10, base_family = "Helvetica")
png("RSEM.png",width = widthcal, height = 600, res = 200)
p
dev.off()
