import datetime
import requests
import json
import pandas
import csv
import subprocess
import PythonMagick
import os
import glob
import shutil
import zipfile
import numpy
from flask import Flask, render_template, request, redirect, current_app
from bokeh.charts import Line, show, output_file, save, Dot, ColumnDataSource
from bokeh.models import Label,HoverTool
from bokeh.plotting import figure, show, output_file

from bokeh.embed import components
from collections import OrderedDict


app = Flask(__name__)
app.debug = True

app.vars=''

@app.route('/')
def main():
  return redirect('/index')

@app.route('/index',methods=['GET','POST'])
def index():
	if request.method == 'GET':
  		return render_template('entry.html')
	else:
		app.vars=request.form['name']
		f=open('history.txt','a')
		f.write(app.vars)
		f.write("\n")
		f.close()
		return redirect('/trying')

@app.route('/trying',methods=['GET','POST'])
def trying():
	##files = glob.glob('.static/OUT/*')
	##print files
	##for f in files:
    	##	os.remove(f)
	gname=app.vars
	
	if os.path.exists('./static/zips/'+app.vars+'.zip'):
		with zipfile.ZipFile("./static/zips/"+app.vars+".zip","r") as zip_ref:
			zip_ref.extractall("./static/OUT/")
		filelink="/static/zips/"+app.vars+".zip"
		pnglist=glob.glob('./static/OUT/*km.png')
		object_list = get_csv()
		script,div=plot_levels()
		return render_template('trying.html', splicing=get_splice(gname), object_list=object_list, filelink=filelink, pnglist=pnglist,gname=gname,script=script, div=div)

	command = 'Rscript'
	path2script = './static/TCGAkm.r'
	args = [app.vars]
	cmd = [command, path2script] + args
	x = subprocess.check_output(cmd, universal_newlines=True)

	img = PythonMagick.Image()
	img.density("300")
	pdflist=glob.glob('./static/OUT/*.pdf')
	pnglist=[]
	for item in pdflist:
		img.read(item)
		newf=item[:-3]+'png'
		img.write(newf)
		pnglist.append(newf)

	path2script = './static/expr_med.r'
	cmd2 = [command, path2script]
	y=subprocess.check_output(cmd2, universal_newlines=True)

	path2script = './static/mutations.r'
	args = [app.vars]
	cmd3 = [command, path2script] + args
	z=subprocess.check_output(cmd3, universal_newlines=True)

	shutil.make_archive("./static/zips/"+app.vars, 'zip', "./static/OUT/")
	filelink="/static/zips/"+app.vars+".zip"
		
	object_list = get_csv()
	script,div=plot_levels()

	return render_template('trying.html', splicing=get_splice(gname), object_list=object_list, filelink=filelink, pnglist=pnglist, gname=gname,script=script, div=div)

def get_csv():
	p = './static/OUT/final.csv'
	f = open(p, 'r')
	return list(csv.DictReader(f))

def plot_levels():
	#read csv
	levelp='./static/OUT/level.csv'
	df=pandas.read_csv(levelp,header='infer')
	df2=df.set_index(df["tumor_type"])
	df2=df2.iloc[:,:3]
	df3=df2.stack()
	df3.index.names=['tumor_type','condition']
	df4=pandas.DataFrame(df3)
	df4.columns=["RSEM"]
	df5=df4.reset_index()
	df5=df5.replace(0,numpy.nan)


	#plot w/ bokeh
	fold=[]
	for n in range(len(df5)/3):
		for i in range(3):
			fold.append(str(df5.iloc[n*3+i,2]/df5.iloc[n*3,2]))
	df5["folds"]=fold
	source=ColumnDataSource(df5)
	hover = HoverTool(tooltips=[
		("RSEM", "@RSEM"),
		("fold/normal", "@folds"),
    		("condition", "@condition"),
    		("type:", "@tumor_type"),
		])
	p=figure(tools=[hover])

	colormap = {'surv_cutoff': 'red', 'normal_median': "green", 'tumor_median': 'blue'}
	colors = [colormap[x] for x in df5["condition"]]
	p.circle("folds","RSEM",color=colors, fill_alpha=0.2, size=15, source=source)
	#li=Dot(df5, values="RSEM", label='tumor_type',
        #   group='condition',legend='bottom_right', ylabel="RSEM(log2)",fill_alpha=1,tools=[hover])

	script, div = components(p)
	return script,div

def get_splice(gname):
	with open('/home/rf/TCGAapp/static/TCGA_alt_spl_g.txt','r') as slist:
		rlist=slist.read().split('\n')
		if app.vars in rlist:
			splice="Significant Alternative Splicing"
		else:
			splice="No Significant Alternative Splicing"
		return splice

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=33507)
