#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#scheduling the targets
import subprocess
import pandas as pd 
from   PyAstronomy  import pyasl

LT_targets_dir = 
LT_targets = 
RA  = []
DEC = []
for i in np.arange(len(LT_targets)):  # ra and dec are converted from degrees to hms and dms respectively for sceduling purpose
        ra      = LT_targets.loc[i,'ra'] # ra of all galaxies
        dec     = LT_targets.loc[i,'dec']  # dec of all galaxies
        sexa    = pyasl.coordsDegToSexa(ra,dec, fmt=('%02d:%02d:%06.2f,', '%s%02d:%02d:%06.2f'))   # ra, dec in sexagimal format are stored as one string; comma separated               
        s       = sexa.split(',')  # splitting the string into two; ra  and dec separately obtained
        RA.append(s[0])  # list of all ra in sexagesimal
        DEC.append(s[1])   # list of all dec in sexagesimal
LT_targets.loc[:,'RA']  = RA # storing the ra list as a column in LT_targets
LT_targets.loc[:,'DEC'] = DEC  # storing the dec list as a column in LT_targets
sceduling_targets       = LT_targets.loc[:,['Name','RA', 'DEC']]  # creating a new dataframe only with information required for sceduling

for filter_choice in ['"G"']:    
        bin_dir      = "/home/extprajk/local_dev/rtml"
        docs_dir     = "/home/extprajk/local_dev/rtml/scripts/rtml_documents"
        project      = "gwem"
        username     = "copperwheat_chris"
        contactname  = "Copperwheat.Chris"
        password     = "WU2DGFmK"
        raoff        = "0.0"
        decoff       = "0.0"
        startdate    = "2016-01-01T00:00:00"
        enddate      = "2017-12-31T00:00:00"
        exposure     = "60"
        nexpose      = "1"
        filt         = filter_choice
        binning      = "2 2" 
        now_date     = str(datetime.time(datetime.now()))
        for i in np.arange(len(sceduling_targets)):
            target   = '"' + str(sceduling_targets.loc[i,'Name']).strip()+ '"'
            ra       = '"' + str(sceduling_targets.loc[i,'RA']).strip()  + '"'
            dec      = '"' + str(sceduling_targets.loc[i,'DEC']).strip() + '"'
            string   = bin_dir + "/create_rtml -rtml_version 3.1a -request -uid rtml://rtml-ioo-" + now_date + " -iauri file:/dev/null -project " + project + " -contact -contact_name " + contactname + " -contact_user " + username + " -observation -target -target_name " + target + " -ra " + ra + " -dec " + dec + " -ra_offset " + raoff + " -dec_offset " + decoff + " -start_date " + startdate + " -end_date " + enddate + " -device -device_name IO:O -device_type camera -device_spectral_region optical -device_filter " + filt + " -binning " + binning + " -exposure " + exposure + " s " + nexpose + " > " + docs_dir +"/ioo_rtml_document_" + str(i) + ".rtml"
            subprocess.call(string, shell = True)
            string2  = bin_dir + "/nodetester.tcsh -url http://161.72.57.3:8080/org_estar_nodeagent/services/NodeAgent -username " + username + " -password " + password + " -cookie -axis -urn  urn:/node_agent -method handle_rtml -rtml \"" + docs_dir +"/ioo_rtml_document_" + str(i) + ".rtml\""
            print (string2,'\n')
            subprocess.call(string2, shell = True)

