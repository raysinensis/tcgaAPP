import datetime
import requests
import json
import pandas
import csv
import subprocess
#import PythonMagick
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


app = Flask(__name__,static_url_path='/static')
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

	#img = PythonMagick.Image()
	#img.density("300")
	#pdflist=glob.glob('./static/OUT/*.pdf')
	#pnglist=[]
	#for item in pdflist:
	#	img.read(item)
	#	newf=item[:-3]+'png'
	#	img.write(newf)
	#	pnglist.append(newf)

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
	levelp='./static/OUT/levels.csv'
	df=pandas.read_csv(levelp,header='infer')
	df["cancers"] = df["tumor_type"] + '_' + df["condition"]
	cancers=list(df["cancers"])
	upperscore=list(df.iloc[:,4])
	lowerscore=list(df.iloc[:,0])
	q1score=list(df.iloc[:,1])
	q2score=list(df.iloc[:,2])
	q3score=list(df.iloc[:,3])
	cutscore=list(df.iloc[:,5])
	
	#zip?

	p = figure(tools="save", background_fill_color="#EFE8E2", title="", x_range=cancers)
	p.segment(cancers, upperscore, cancers, q3score, line_color="black")
	p.segment(cancers, lowerscore, cancers, q1score, line_color="black")
	p.vbar(cancers, 0.7, q2score, q3score, fill_color="#E08E79", line_color="black")
	p.vbar(cancers, 0.7, q1score, q2score, fill_color="#3B8686", line_color="black")
	p.rect(cancers, cutscore, 0.7, 0.01, line_color="red")
	p.rect(cancers, lowerscore, 0.2, 0.01, line_color="black")
	p.rect(cancers, upperscore, 0.2, 0.01, line_color="black")

	p.xgrid.grid_line_color = None
	p.ygrid.grid_line_color = "white"
	p.grid.grid_line_width = 2
	p.xaxis.major_label_text_font_size="12pt"

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
