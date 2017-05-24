with open('/home/rf/Desktop/TCGA_alt_spl.csv','r') as f:
	readd=f.read()
lines=readd.split("\n")
genel=[]
for line in lines:
	allgene=line.split(",")
	for onegene in allgene:
		if onegene!='':
			genel.append(onegene)
genes=set(genel)
genes=list(genes)
with open('/home/rf/Desktop/TCGA_alt_spl_g.txt','w') as f:
	for item in genes:
		f.write("%s\n" %item)
