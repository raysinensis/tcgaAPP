from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
#from sklearn.metrics import precision_recall_fscore_support
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
from puAdapter import PUAdapter
#from sklearn_pandas import DataFrameMapper


df=pd.read_csv('/home/rf/Desktop/tcga9-final_Clabels.csv')
#fill missing data
df.iloc[:,1:5]=df.iloc[:,1:5].fillna(0)
df.iloc[:,17:21]=df.iloc[:,17:21].fillna(1)
##convert to 0 and 1
df['methyl.BRCA'][df['methyl.BRCA']<=0.1]=-1
df['methyl.SKCM'][df['methyl.SKCM']<=0.1]=-1
df['methyl.COAD'][df['methyl.COAD']<=0.1]=-1
df['methyl.LUAD'][df['methyl.LUAD']<=0.1]=-1

>>> df['methyl.BRCA'][df['methyl.BRCA']!=-1]=0
>>> df['methyl.SKCM'][df['methyl.SKCM']!=-1]=0
>>> df['methyl.COAD'][df['methyl.COAD']!=-1]=0
>>> df['methyl.LUAD'][df['methyl.LUAD']!=-1]=0

>>> df['methyl.LUAD'][df['methyl.LUAD']==-1]=1
>>> df['methyl.COAD'][df['methyl.COAD']==-1]=1
>>> df['methyl.SKCM'][df['methyl.SKCM']==-1]=1
>>> df['methyl.BRCA'][df['methyl.BRCA']==-1]=1

#getting y and X
y=df['CANCER']
X=df.drop(['CANCER',"Gene.Name"],axis=1)
y=y.tolist()
y=[-1 if i==0 else i for i in y]
y=np.asarray(y)
#y=pd.DataFrame(y)
#X_pos=df[df['CANCER']==1].drop(['CANCER','Gene.Name'],axis=1)#[:300]
#y_pos=y[y==1]#[:300]
#shuffle and sample
#shuffled = np.random.permutation(len(df))
#y = y[shuffled]
#X = X.iloc[shuffled]
#split = int(0.75*len(y))
#X_train = X.iloc[:split]
#y_train = y[:split]
#X_test = X.iloc[split:]
#y_test = y[split:]
#X_1500=X.iloc[:1500]
#y_1500=y[:1500]
#X_train=X_1500.append(X_pos,ignore_index=True)
#y_train=y_3000.append(y_pos,ignore_index=True)
X_t=pd.DataFrame.as_matrix(X)
#y_t=pd.DataFrame.as_matrix(y_train)
#y_t=y_t.tolist()
#y_t=[-1 if y==0 else y for y in y_t]
#y_t=np.asarray(y_t)

#y_pos=y_pos.tolist()
#y_1500=y_1500.tolist()
#y_train=y_1500+y_pos
y_t=np.array(y)
X_s=pd.DataFrame.as_matrix(X_test)
y_s=np.array(y_test)


estimator = RandomForestClassifier(n_estimators=55,
                                           criterion='gini', 
                                           bootstrap=False,
                                           n_jobs=2)
pu_estimator = PUAdapter(estimator, hold_out_ratio=0.2)
pu_estimator.fit(X_t,y_t)
y_pred = pu_estimator.predict(X_s)
diff=y_s[y_pred!=y_s]
print n,list(diff).count(-1),list(diff).count(1)
falseneg=list(diff).count(1)
posfind=list(diff).count(-1)

#searching for parameters
#136/4365 in testgroup,33 false negative,2764 new positive
#569/17458 in total, 69 false negative, 1343 new positive
for n in np.arange(1,20,2):
	estimator = RandomForestClassifier(n_estimators=n,#9?
                                           criterion='entropy', 
                                           bootstrap=True,
                                           n_jobs=6,
					   #min_samples_leaf=n,
					   oob_score=True,
					   #max_features='log2',
					   class_weight='balanced',)
	pu_estimator = PUAdapter(estimator, hold_out_ratio=0.2)
	pu_estimator.fit(X_t,y_t)
	y_pred = pu_estimator.predict(X_t)
	diff=y_t[y_pred!=y_t]
	print n,list(diff).count(-1),list(diff).count(1),estimator.oob_score_
estimator.feature_importances_

#without puadapter, 98,0	
for n in np.arange(1,20,2):
	estimator = RandomForestClassifier(n_estimators=n,#9?
                                           criterion='entropy', 
                                           bootstrap=True,
                                           n_jobs=6,
					   #min_samples_leaf=n,
					   oob_score=True,
					   #max_features='log2',
					   class_weight='balanced',)
	estimator.fit(X_t,y_t)
	y_pred = estimator.predict(X_t)
	diff=y_t[y_pred!=y_t]
	print n,list(diff).count(-1),list(diff).count(1),estimator.oob_score_

#actual run
estimator = RandomForestClassifier(n_estimators=15,
                                           criterion='entropy', 
                                           bootstrap=True,
                                           n_jobs=6,
					   #min_samples_leaf=n,
					   oob_score=True,
					   #max_features='log2',
					   class_weight='balanced',)
pu_estimator = PUAdapter(estimator, hold_out_ratio=0.2)
pu_estimator.fit(X_t,y_t)
y_pred = pu_estimator.predict(X_t)
diff=y_t[y_pred!=y_t]
print n,list(diff).count(-1),list(diff).count(1),estimator.oob_score_
genes=df[y_pred!=y_t]
genepos=genes[genes['CANCER']==0]
geneneg=genes[genes['CANCER']==1]


with open("pred_pos.txt",'w') as f:
	for element in list(genepos['Gene.Name']):
		f.write(element+'\n')
with open("pred_neg.txt",'w') as f:
	for element in list(geneneg['Gene.Name']):
		f.write(element+'\n')

