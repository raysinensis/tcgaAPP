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

from flask import Flask, render_template, request, redirect, current_app, url_for
from werkzeug.utils import secure_filename

from bokeh.charts import Line, show, output_file, save, Dot, ColumnDataSource
from bokeh.models import Label,HoverTool,Range1d,TapTool,OpenURL
from bokeh.plotting import figure, show, output_file
from bokeh.layouts import widgetbox,column,row
from bokeh.models.widgets import CheckboxGroup,Slider
from bokeh.models.callbacks import CustomJS
from bokeh.embed import components

from collections import OrderedDict


app = Flask(__name__,static_url_path='/static')
app.debug = True

UPLOAD_FOLDER = './static/upload'
ALLOWED_EXTENSIONS = set(['txt', 'csv'])

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

app.vars=''

@app.route('/')
def main():
  return redirect('/index')

@app.route('/<string:page_name>/')
def static_page(page_name):
    return render_template('%s.html' % page_name)

@app.route('/index',methods=['GET','POST'])
def index():
	if request.method == 'GET':
  		return render_template('entry.html')
	else:
		app.vars=request.form['name']
		if len(app.vars)>2:
			f=open('history.txt','a')
			f.write(app.vars)
			f.write("\n")
			f.close()
			return redirect('/trying')
		else:
			return render_template('entry.html')

@app.after_request
def add_header(r):
    r.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    r.headers["Pragma"] = "no-cache"
    r.headers["Expires"] = "0"
    r.headers['Cache-Control'] = 'public, max-age=0'
    return r

@app.route('/listquery',methods=['GET','POST'])
def listquery():
	if request.method == 'POST':
        # check if the post request has the file part
		if 'file' not in request.files:
		    result=[]
        file = request.files['file']
        # if user does not select file, browser also
        # submit a empty part without filename
        if file.filename == '':
            result=[]
        elif file and allowed_file(file.filename):
        	filename = secure_filename(file.filename)
        	file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
		readpath=os.path.join(app.config['UPLOAD_FOLDER'], filename)
		df=pandas.read_csv(readpath,header=None)
		result=list(df[0])
        
	makejson(result)
	return render_template("listquery.html")
		
@app.route('/genemap',methods=['GET','POST'])
def genemap():
	return render_template('genemap.html')

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
	pdflist=glob.glob('./static/OUT/*.pdf')
	pnglist=[]
	for item in pdflist:
		#img.read(item)
		#newf=item[:-3]+'png'
		#img.write(newf)
		#pnglist.append(newf)
		os.rename(item,item[:-7]+'_tumor-km.pdf')

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

def allowed_file(filename):
    """Does filename have the right extension?"""
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def get_csv():
	p = './static/OUT/final.csv'
	f = open(p, 'r')
	return list(csv.DictReader(f))

def plot_levels():
	#read csv
	levelp='./static/OUT/levels.csv'
	df=pandas.read_csv(levelp,header='infer')
	df=df.replace(numpy.nan,0)
	df["cancers"] = df["tumor_type"] + '_' + df["condition"]
	cancers=list(df["cancers"])
	maxy=max(df["max"])+0.5
	
	#upperscore=list(df.iloc[:,4])
	#lowerscore=list(df.iloc[:,0])
	#q1score=list(df.iloc[:,1])
	#q2score=list(df.iloc[:,2])
	#q3score=list(df.iloc[:,3])
	#cutscore=list(df.iloc[:,5])
	#folds=list(df.iloc[:,5]/df.iloc[:,2])
	
	#zip?
	q2ratio=[]
	for i in range(len(df)):
		q2ratio.append(2**(df["median"][i]-df["median"][i//2*2]))
	minratio=[]
	for i in range(len(df)):
		minratio.append(2**(df["min"][i]-df["min"][i//2*2]))
	maxratio=[]
	for i in range(len(df)):
		maxratio.append(2**(df["max"][i]-df["max"][i//2*2]))
	cutratio=[]
	for i in range(len(df)):
		cutratio.append(2**(df["surv_cutoff"][i]-df["median"][i//2*2]))
	source=ColumnDataSource(
		data=dict(
		cancers=list(df["cancers"]),
		upperscore=list(df.iloc[:,4]),
		lowerscore=list(df.iloc[:,0]),
		q1score=list(df.iloc[:,1]),
		q2score=list(df.iloc[:,2]),
		q3score=list(df.iloc[:,3]),
		cutscore=list(df.iloc[:,5]),
		ratio=q2ratio
		)
	)
	source2=ColumnDataSource(
		data=dict(
		cancers=list(df["cancers"]),
		upperscore=list(df.iloc[:,4]),
		lowerscore=list(df.iloc[:,0]),
		q1score=list(df.iloc[:,1]),
		q2score=list(df.iloc[:,2]),
		q3score=list(df.iloc[:,3]),
		cutscore=list(df.iloc[:,5]),
		ratio=cutratio,
		)
	)
	source3=ColumnDataSource(
		data=dict(
		cancers=list(df["cancers"]),
		upperscore=list(df.iloc[:,4]),
		lowerscore=list(df.iloc[:,0]),
		q1score=list(df.iloc[:,1]),
		q2score=list(df.iloc[:,2]),
		q3score=list(df.iloc[:,3]),
		cutscore=list(df.iloc[:,5]),
		ratio=minratio,
		)
	)
	source4=ColumnDataSource(
		data=dict(
		cancers=list(df["cancers"]),
		upperscore=list(df.iloc[:,4]),
		lowerscore=list(df.iloc[:,0]),
		q1score=list(df.iloc[:,1]),
		q2score=list(df.iloc[:,2]),
		q3score=list(df.iloc[:,3]),
		cutscore=list(df.iloc[:,5]),
		ratio=maxratio,
		)
	)
	
	p = figure(tools=["save","hover",'tap'], background_fill_color="white", title="", x_range=cancers,plot_width=880)
	hover = p.select(dict(type=HoverTool))
	hover.tooltips = [('condition','@cancers'),('ratio/normal','@ratio{1.11}')]
	hover.point_policy='snap_to_data'
	p.segment('cancers', 'upperscore', 'cancers', 'q3score', line_color="black",source=source,legend="max")
	p.segment('cancers', 'lowerscore', 'cancers', 'q1score', line_color="black",source=source,legend="min")
	p.vbar('cancers', 0.35, 'q1score', 'q3score', fill_color=["#9ACD32","#FF4500"]*(len(df)/2), line_color="black",source=source,legend="quartiles")
	p.rect('cancers', 'q2score', 0.35, 0.045, line_color="black",fill_color="black",source=source,legend="median")
	p.rect('cancers', 'cutscore', 0.35, 0.045, line_color=["#FFFFFF","blue"]*(len(df)/2),fill_color="blue",source=source2,legend="cutpoint--link to KM plot",name='cut')
	p.rect('cancers', 'lowerscore', 0.1, 0.023, line_color="black",fill_color="black",source=source3,legend="min")
	p.rect('cancers', 'upperscore', 0.1, 0.023, line_color="black",fill_color="black",source=source4,legend="max")

	p.xgrid.grid_line_color = None
	p.ygrid.grid_line_color = "white"
	p.grid.grid_line_width = 2
	p.xaxis.major_label_text_font_size="12pt"
	p.yaxis.major_label_text_font_size="12pt"
	p.yaxis.axis_label_text_font_size="12pt"
	p.y_range=Range1d(0,maxy)
	p.yaxis.axis_label="RSEM(log2)"
	p.yaxis.axis_label_text_font_style = "normal"
	p.xaxis.major_label_orientation=numpy.pi/4

	p.legend.location = "bottom_left"
	p.legend.click_policy="hide"

	url = "/static/OUT/@cancers"
	url2= url+'-km.pdf'
	taptool = p.select(type=TapTool)
	taptool.names=['cut']
	taptool.callback = OpenURL(url=url2)

	script, div = components(p)
	return script,div

def get_splice(gname):
	with open('./static/TCGA_alt_spl_g.txt','r') as slist:
		rlist=slist.read().split('\n')
		if app.vars in rlist:
			splice="Significant Alternative Splicing"
		else:
			splice="No Significant Alternative Splicing"
		return splice

def makejson(listof):
	listcosmic=pandas.read_csv('./static/COSMICcsv.csv',header=None)
	listcosmic=list(listcosmic[0])
	listpred=pandas.read_csv('./static/pred_pos.txt',header=None)
	listpred=list(listpred[0])
	if listof!=[]:
		list1=list(set(listof).intersection(listcosmic))
		list2=list(set(listof).intersection(listpred))
		list3=list(set(listof)-set(list1)-set(list2))
	else:
		list1=listcosmic
		list2=listpred
		list3=[]
	d1=[]
	for gene in list1:
		d1.append(("name",gene))
	d2=[]
	for gene in list2:
		d2.append(("name",gene))
	d3=[]
	for gene in list3:
		d3.append(("name",gene))
	count1=len(d1)
	count2=len(d2)
	count3=len(d3)
	
	dict1='%s' % ',\n'.join(['{{"{}": {}}}'.format(action, json.dumps(dictionary)) for action, dictionary in d1])
	dict2='%s' % ',\n'.join(['{{"{}": {}}}'.format(action, json.dumps(dictionary)) for action, dictionary in d2])
	dict3='%s' % ',\n'.join(['{{"{}": {}}}'.format(action, json.dumps(dictionary)) for action, dictionary in d3])
	c1='''{
 "name": "From List",
 "children": [
  {
   "name": "Curated Genes",
   "children": ['''.replace("Genes","Genes -- "+str(count1))
	c2=''']
  },
  {
   "name": "Predicted Genes",
   "children": ['''.replace("Genes","Genes -- "+str(count2))
	c3='''   ]
  },
  {
   "name": "Other Genes",
   "children": ['''.replace("Genes","Genes -- "+str(count3))
	c4='''   ]
  }
 ]
}'''
	final=c1+dict1+c2+dict2+c3+dict3+c4
	with open('./static/flare.json','w') as f:
		f.write(final)
	with open('./static/listresult.txt','w') as f:
		f.write('\n'.join(["curated:",'\n'.join(list1),"predicted:",'\n'.join(list2),"other:",'\n'.join(list3)]))
		

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=33507)
