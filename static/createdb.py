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


