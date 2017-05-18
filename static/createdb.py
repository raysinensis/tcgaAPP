from sqlalchemy import create_engine,MetaData
from sqlalchemy.ext.automap import automap_base
import pandas as pd
##csv to sql
engine = create_engine('sqlite:///static/database/methyl.db')
df = pd.read_csv('./static/methylation db.csv')
df.to_sql(name='methyl', con=engine)
##query from db
metadata = MetaData(engine)
Base = automap_base()
Base.prepare(engine, reflect=True)
#engine.table_names() ##checking table names
#methyldb = Table('methyl',metadata, autoload=True)
#print methyldb.c ##to see column names

Session = sessionmaker(bind=engine)
methyldb = Session()

from sqlalchemy import text
from sqlalchemy.orm import sessionmaker
engine = create_engine('sqlite:///static/database/methyl.db')
Session = sessionmaker(bind=engine)
methyldb = Session()
gene = 'ZFP36L1'
qcol = ['BRCA','COAD','GBM','KICH','LUAD','PAAD','SARC','STAD']
qcolstr = ','.join(qcol)
sqlstr = 'select '+qcolstr+ ' from methyl where gene=\"'+gene+'\"'
sqlcmd = text(sqlstr)
result = methyldb.execute(sqlcmd).fetchall()

##cox coeff from xlsx
trying=pd.read_excel('/home/rf/Downloads/peerj-03-1499-s001.xlsx',sheetname=None )
tempgenes=[]
for cancer in trying.keys():
	tempgenes.extend((trying[cancer])['Gene Name'].tolist())
dfcox=pd.DataFrame()
dfcox['Gene Name']=list(set(tempgenes))
for cancer in trying.keys():
	tempdf=(trying[cancer])[['Gene Name','Raw Cox Coefficient']].drop_duplicates('Gene Name')
	dfcox=dfcox.merge(tempdf,on='Gene Name',how='left')
	dfcox.rename(columns={'Raw Cox Coefficient':cancer},inplace=True)
dfcox.rename(columns={'Gene Name':'Gene'},inplace=True)
engine = create_engine('sqlite:///static/database/methyl.db')
dfcox.to_sql(name='cox', con=engine)
