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
#from flask_weasyprint import HTML, render_pdf
from werkzeug.utils import secure_filename
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import text
from config import ADMINS, MAIL_SERVER, MAIL_PORT, MAIL_USERNAME, MAIL_PASSWORD

from bokeh.charts import Line, show, output_file, save, Dot, ColumnDataSource
from bokeh.models import Label,HoverTool,Range1d,TapTool,OpenURL
from bokeh.plotting import figure, show, output_file
from bokeh.layouts import widgetbox,column,row
from bokeh.models.widgets import CheckboxGroup,Slider
from bokeh.models.callbacks import CustomJS
from bokeh.embed import components

from collections import OrderedDict
import networkx as nx
from networkx.readwrite import json_graph

app = Flask(__name__,static_url_path='/static')
#app.debug = True

UPLOAD_FOLDER = './static/upload'
ALLOWED_EXTENSIONS = set(['txt', 'csv'])

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///static/database/tcga.db'
db = SQLAlchemy(app)

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
		app.vars=app.vars.upper()
		if len(app.vars)>2:
			f=open('history.txt','a')
			f.write(app.vars)
			f.write("\n")
			f.close()
			return redirect(url_for('trying',name=app.vars))
		else:
			return render_template('entry.html')

def create_pdf(pdf_data):
    pdf = StringIO()
    pisa.CreatePDF(StringIO(pdf_data), pdf)
    return pdf

def validgene(gname):
	sqlstr = 'select * from cox where gene=\"'+gname+'\"'
	sqlcmd = text(sqlstr)
	result = db.engine.execute(sqlcmd).fetchall()
	if len(result)!=0:
		return True
	else:
		return False

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
	app.vars = request.args.get('name').upper()
	gname=app.vars
	if not validgene(app.vars):
		return render_template('entry-err.html')
	if os.path.exists('./static/zips/'+app.vars+'.zip'):
		with zipfile.ZipFile("./static/zips/"+app.vars+".zip","r") as zip_ref:
			zip_ref.extractall("./static/OUT/")
		filelink="/static/zips/"+app.vars+".zip"
		pnglist=glob.glob('./static/OUT/*km.png')
		object_list = get_csv(gname)
		script,div=plot_levels()
		getsubgraph(app.vars)
		rend=render_template('trying.html', splicing=get_splice(gname), object_list=object_list, filelink=filelink, pnglist=pnglist,gname=gname,script=script, div=div)
		#HTML(string=rend).write_pdf('./static/OUT/gene.pdf')
		return rend

	#files = glob.glob('.static/OUT/*')
	##print files
	#for f in files:
    	#	os.remove(f)
	
	command = 'Rscript'
	path2script = './static/TCGAkm.r'
	args = [app.vars]
	cmd = [command, path2script] + args
	x = subprocess.check_output(cmd, universal_newlines=True)

	#img = PythonMagick.Image()
	#img.density("300")
	pdflist=glob.glob('./static/OUT/*km.pdf')
	pnglist=[]
	for item in pdflist:
		#img.read(item)
		#newf=item[:-3]+'png'
		#img.write(newf)
		#pnglist.append(newf)
		os.rename(item,item[:-7]+'_tumor.pdf')

	path2script = './static/expr_med.r'
	cmd2 = [command, path2script]
	y=subprocess.check_output(cmd2, universal_newlines=True)

	path2script = './static/mutations.r'
	args = [app.vars]
	cmd3 = [command, path2script] + args
	z=subprocess.check_output(cmd3, universal_newlines=True)

	shutil.make_archive("./static/zips/"+app.vars, 'zip', "./static/OUT/")
	filelink="/static/zips/"+app.vars+".zip"
		
	object_list = get_csv(gname)
	script,div=plot_levels()
	getsubgraph(app.vars)

	rend=render_template('trying.html', splicing=get_splice(gname), object_list=object_list, filelink=filelink, pnglist=pnglist,gname=gname,script=script, div=div)
	#HTML(string=rend).write_pdf('./static/OUT/gene.pdf')
	return rend

def allowed_file(filename):
    """Does filename have the right extension?"""
    return '.' in filename and filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def get_csv(gname):
	gene = gname
	qcol = ['BLCA','BRCA','CESC','COAD','GBM','HNSC','KIRC','KIRP','LIHC','LUAD','LUSC','SKCM','STAD']
	qcolstr = ','.join(qcol)
	##methylation
	sqlstr = 'select '+qcolstr+ ' from methyl where gene=\"'+gene+'\"'
	sqlcmd = text(sqlstr)
	result = db.engine.execute(sqlcmd).fetchall()
	if len(result)==0:
		result= [[1]*len(qcol)]
	##cox coeff
	sqlstr = 'select '+qcolstr+ ' from cox where gene=\"'+gene+'\"'
	sqlcmd = text(sqlstr)
	result2 = db.engine.execute(sqlcmd).fetchall()
	if len(result2[0])==0:
		result2 = [['NaN']*len(qcol)]
	p = './static/OUT/final.csv'
	with open(p, 'r') as f:
		count=0
		lines=list(csv.DictReader(f))
		newlines=[]
		for line in lines:
			if result[0][count]<=0.1:
				metp='MET'
			else:
				metp='-'
			line.update({'V11':metp})
			line.update({'V5':result2[0][count]})
			newlines.append(line)
			count+=1
		return newlines

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
	
	p = figure(tools=["save","hover",'tap'], background_fill_color="white", title="", x_range=cancers,plot_width=920)
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
	url2= url+'.pdf'
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

def getsubgraph(gene):
	listcosmic=pandas.read_csv('./static/COSMICcsv.csv',header=None)
	listcosmic=list(listcosmic[0])
	listpred=pandas.read_csv('./static/pred_pos.txt',header=None)
	listpred=list(listpred[0])

	G=nx.read_edgelist("./static/edgelist",delimiter='\t')
	try:
		listof=list(G[gene])
		listof.append(gene)
		H=G.subgraph(listof)
	except:
		listof=[gene]
		H=nx.Graph()
		H.add_node(gene)
	for n in H:
		H.node[n]['name'] = n
	for n in H:
		if H.node[n]['name'] in listcosmic:
			H.node[n]['group'] =2
		elif H.node[n]['name'] in listpred:
			H.node[n]['group'] =1
		else:
			H.node[n]['group'] =3
	d = json_graph.node_link_data(H)
	json.dump(d, open('./static/onegene.json','w'))


import smtplib
import logging
from logging.handlers import SMTPHandler
class TlsSMTPHandler(SMTPHandler):
    def emit(self, record):
        try:
            import string # for tls add this line
            try:
                from email.utils import formatdate
            except ImportError:
                formatdate = self.date_time
            port = self.mailport
            if not port:
                port = smtplib.SMTP_PORT
            smtp = smtplib.SMTP(self.mailhost, 587)
            msg = self.format(record)
            msg = "From: %s\r\nTo: %s\r\nSubject: %s\r\nDate: %s\r\n\r\n%s" % (
                            self.fromaddr,
                            string.join(self.toaddrs, ","),
                            self.getSubject(record),
                            formatdate(), msg)
            if self.username:
                smtp.ehlo() # for tls add this line
                smtp.starttls() # for tls add this line
                #smtp.ehlo() # for tls add this line
                smtp.login(self.username, self.password)
            smtp.sendmail(self.fromaddr, self.toaddrs, msg)
            smtp.quit()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)
if not app.debug:
#email alert
    credentials = None
    if MAIL_USERNAME or MAIL_PASSWORD:
        credentials = (MAIL_USERNAME, MAIL_PASSWORD)
    mail_handler = TlsSMTPHandler(mailhost=(MAIL_SERVER, MAIL_PORT), fromaddr='no-reply@raysinensis.com', toaddrs=ADMINS, subject='TCGAapp failure', credentials=credentials)
    mail_handler.setLevel(logging.ERROR)
    app.logger.addHandler(mail_handler)
#logging to file
    from logging.handlers import RotatingFileHandler
    file_handler = RotatingFileHandler('static/logs/errors.log', 'a', 1 * 1024 * 1024, 10)
    file_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'))
    app.logger.setLevel(logging.INFO)
    file_handler.setLevel(logging.INFO)
    app.logger.addHandler(file_handler)
    app.logger.info('TCGAapp startup')

@app.route('/static/about.html',methods=['POST'])
def feedback():
	textfb = '\n'+request.form['text']
	server = smtplib.SMTP(MAIL_SERVER, MAIL_PORT)
        server.ehlo()
        server.starttls()
        server.login(MAIL_USERNAME, MAIL_PASSWORD)
        server.sendmail(MAIL_USERNAME, ADMINS, textfb)
        server.close()
	return ('', 204)

@app.errorhandler(500)
def page_not_found(e):
	return redirect('/sorry')
@app.route('/sorry',methods=['GET', 'POST'])
def sorry():
	if request.method == 'GET':
		return render_template('sorry.html'), 500
	else:
		textfb = '\n'+request.form['text']
		server = smtplib.SMTP(MAIL_SERVER, MAIL_PORT)
		server.ehlo()
		server.starttls()
		server.login(MAIL_USERNAME, MAIL_PASSWORD)
		server.sendmail(MAIL_USERNAME, ADMINS, textfb)
		server.close()
		return redirect('/index')

if __name__ == '__main__':
  app.run(host='0.0.0.0', port=33507)
