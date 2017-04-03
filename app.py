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
from flask import Flask, render_template, request, redirect, current_app
from bokeh.charts import Line, show, output_file, save
from bokeh.models import Label
from bokeh.embed import components

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
	if os.path.exists('./static/zips/'+app.vars+'.zip'):
		with zipfile.ZipFile("./static/zips/"+app.vars+".zip","r") as zip_ref:
			zip_ref.extractall("./static/OUT/")
		filelink="/static/zips/"+app.vars+".zip"
		object_list = get_csv()
		return render_template('trying.html', object_list=object_list, filelink=filelink)

	command = 'Rscript'
	path2script = './static/TCGAkm.r'
	args = [app.vars]
	cmd = [command, path2script] + args
	x = subprocess.check_output(cmd, universal_newlines=True)

	img = PythonMagick.Image()
	img.density("300")
	img.read('./static/OUT/km.pdf')
	img.write('./static/OUT/km.png')

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
	return render_template('trying.html', object_list=object_list, filelink=filelink)

def get_csv():
	p = './static/OUT/final.csv'
	f = open(p, 'r')
	return list(csv.DictReader(f))

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=33507)
