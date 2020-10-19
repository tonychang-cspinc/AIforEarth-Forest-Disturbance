import pandas as pd
import numpy as np
cn=48202851010497

treecsv=pd.read_csv('csv_data/CA_TREE.csv', usecols=('DIA', 'DRYBIO_AG','STATUSCD','PLT_CN','SUBP','INVYR', 'TREE'))
condcsv=pd.read_csv('csv_data/CA_COND.csv', usecols=('LIVE_CANOPY_CVR_PCT','CONDPROP_UNADJ','INVYR','PLT_CN','COND_STATUS_CD' ))
plotcsv=pd.read_csv('csv_data/CA_PLOT.csv', usecols=('PLOT_STATUS_CD','INVYR','CN', 'LAT','LON'))
subset=treecsv[treecsv['PLT_CN']==cn]
condsubset=condcsv[condcsv['PLT_CN']==cn]
plotsubset=plotcsv[plotcsv['CN']==cn]

print(subset)
print(condsubset)
print(plotsubset)
#lat, lon = plotcsv['LAT'].tolist(), plotcsv['LON'].tolist()
#latlon=list(zip(lat,lon))
#print(latlon)
#newdf=pd.DataFrame({'LAT':plotcsv['LAT'],'LON':plotcsv['LON']})

#print(len(latlon),len(set(latlon)))
print(str(len(subset)) + ' total trees.')
