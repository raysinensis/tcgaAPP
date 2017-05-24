##download data w/ cmd: firehose_get -tasks Clinical_vs_Methylation analyses latest
import os
import tarfile
import csv
import pandas
path=os.getcwd()
cancerlist=os.listdir(".")

#get genes list
folderpath2=path+"/"+"LAML"+"/"
csvpath=folderpath2+"gdac.broadinstitute.org_"+"LAML"+"-TB.Correlate_Clinical_vs_Methylation.Level_4.2016012800.0.0/supp.table1.txt"
csv1 = pandas.read_csv(csvpath,sep='\t')
df=pandas.DataFrame()
df["gene"]=csv1.iloc[:,0]
#pull 3rd column from each
for cancer in cancerlist:
	folderpath=path+"/"+cancer+"/20160128/"
	folderpath2=path+"/"+cancer+"/"
	csvpath=folderpath2+"gdac.broadinstitute.org_"+cancer+"-TP.Correlate_Clinical_vs_Methylation.Level_4.2016012800.0.0/supp.table1.txt"
	try:
		csv1 = pandas.read_csv(csvpath,sep='\t')
	except IOError,e:
		csvpath=folderpath2+"gdac.broadinstitute.org_"+cancer+"-TB.Correlate_Clinical_vs_Methylation.Level_4.2016012800.0.0/supp.table1.txt"
		try:
			csv1 = pandas.read_csv(csvpath,sep='\t')
		except IOError,e2:
			csvpath=folderpath2+"gdac.broadinstitute.org_"+cancer+"-TM.Correlate_Clinical_vs_Methylation.Level_4.2016012800.0.0/supp.table1.txt"
			csv1 = pandas.read_csv(csvpath,sep='\t')
	
	df[cancer]=csv1.iloc[:,1]
df.to_csv("methylation.csv")
	
