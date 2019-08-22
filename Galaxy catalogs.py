#!/usr/bin/env python
# coding: utf-8

# In[128]:


import numpy as np
import pandas as pd
get_ipython().run_line_magic('config', 'InlineBackend.rc = {}')
import matplotlib
#matplotlib.rc_file("../../templates/matplotlibrc")
import astropy.coordinates as coord
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from   astropy.table  import Table
import astropy.units as u
import math
GLADE_dir           = "/mnt/d/LT_algorithm/GLADE_2.3.txt"
GLADE               = pd.read_csv(GLADE_dir, delimiter=' ',header = None)
GLADE.columns       = ['PGC','name', 'HyperLEDA name','2MASS name','SDSS-DR12 name',
                       'flag1','ra','dec','dist','dist_err','z','Bmag', 'Bmag_err','BMAG',
                       'Jmag','Jmag_err','Hmag','Hmag_err','Kmag','Kmag_err','flag2','flag3']  
GLADE            = GLADE.loc[:,['name','ra','dec', 'BMAG']]
#--------


# In[130]:


#------------------------------------------------------------------------
GWGC_dir            = "/mnt/d/LT_algorithm/GWGC_catalog.txt"
GWGC                = pd.read_csv(GWGC_dir, sep='\t') # galaxy list from NED
GWGC                = GWGC.replace(r'^\s+$', np.nan, regex=True)
GWGC.columns       = ['name', 'ra', 'dec', 'TT', 'Bmag', 'e_Bmag', 'a', 'e_a', 'b',
                    'e_b', 'b/a', 'e_b/a', 'PA', 'BMAG', 'e_BMAG', 'dist', 'e_Dist']
for col_name in GWGC.columns[1:]:  # changing the datatypes of the columns form object to float 
     GWGC.loc[:,col_name] = GWGC.loc[:,col_name].astype(float)
GWGC             = GWGC.loc[:,['name','ra','dec','BMAG']]
#-------------------------------------------------------------------------------


# In[148]:


def BMAG(CLU_NED):
    Dist =  CLU_NED.loc[:,'DISTMPC']
    Bmag  = CLU_NED.loc[:,'MAGB']
    BMAG = []
    for i in np.arange(len(CLU_NED)):
        mag    = 5 + Bmag[i]- 5*math.log10(Dist[i])   
        BMAG.append(mag)
    CLU_NED['BMAG'] = BMAG
    return CLU_NED    


CLU_NED_dir         = "/mnt/d/LT_algorithm/CLU_NEDonly.fits"
data                = Table.read(CLU_NED_dir, format = 'fits')
CLU_NED             = data.to_pandas()
for col_name in CLU_NED.columns:                                               # some column of the df has uncessary character (b'..') 
    value = CLU_NED.loc[0,col_name]    
    value = str(value)
    if value[0] == 'b':           
        CLU_NED[col_name] = CLU_NED[col_name].astype(str)                       # first all the values are converted from object to string to split.
        new               = CLU_NED[col_name].str.split("'", expand = True)     # the string is split
        CLU_NED[col_name] = new[1] 


no_BMAG = (CLU_NED.loc[:,'MAGB']!= 'NaN')           # removing galaxies that do not have absolute Blue mag
CLU_NED  = CLU_NED[no_BMAG]
CLU_NED = BMAG(CLU_NED)
CLU_NED = CLU_NED.loc[:,['NAME','RA','DEC','BMAG']] 


# In[ ]:





# In[118]:



no_BMAG = (GLADE.loc[:,'BMAG']!= 'NaN' )&  (GLADE.loc[:,'BMAG'] > -20.47)            # removing galaxies that do not have absolute Blue mag
GLADE  = GLADE[no_BMAG]
data= Table.from_pandas(GLADE)
ra = coord.Angle(data['ra'].filled(np.nan)*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(data['dec'].filled(np.nan)*u.degree)
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide" )
ax.scatter(ra.radian, dec.radian, s = 0.0005,color='red')
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())


# In[137]:


no_BMAG = (CLU_NED.loc[:,'BMAG'] > -20.47)            # removing galaxies that do not have absolute Blue mag
CLU_NED  = CLU_NED[no_BMAG]
data2 = Table.from_pandas(CLU_NED)
ra = coord.Angle(data2['RA']*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(data2['DEC']*u.degree)
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide" )
ax.scatter(ra.radian, dec.radian, s = 0.0005,color='red')
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())


# In[138]:


no_BMAG = (GWGC.loc[:,'BMAG'] > -20.47)            # removing galaxies that do not have absolute Blue mag
GWGC  = GWGC[no_BMAG]
data3 = Table.from_pandas(GWGC)
ra = coord.Angle(data3['ra']*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(data3['dec']*u.degree)
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide" )
ax.scatter(ra.radian, dec.radian, s = 0.0005,color='red')
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())


# In[90]:


fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide")
ax.hexbin(ra.radian, dec.radian, cmap=plt.get_cmap('Spectral_r',50))
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_major_formatter(plt.NullFormatter())
fig.colorbar(ax)


# In[ ]:




