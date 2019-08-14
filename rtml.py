#!/usr/bin/env python

import subprocess
from datetime import datetime

bin_dir = "/home/cmc/local_dev/rtml"
docs_dir = "/home/cmc/local_dev/rtml/scripts/rtml_documents"

project = "gwem"
username = "copperwheat_chris"
contactname = "Copperwheat.Chris"
password = "WU2DGFmK"

target = "MyStar"
ra = "13:48:08.60"
dec = "+02:27:55.20"
raoff = "0.0"
decoff = "0.0"
startdate = "2016-01-01T00:00:00"
enddate = "2017-12-31T00:00:00"
exposure = "60"
nexpose = "2"
filt = "G"
binning = "2 2"

now_date = str(datetime.time(datetime.now()))

string = bin_dir + "/create_rtml -rtml_version 3.1a -request -uid rtml://rtml-ioo-" + now_date + " -iauri file:/dev/null -project " + project + " -contact -contact_name " + contactname + " -contact_user " + username + " -observation -target -target_name " + target + " -ra " + ra + " -dec " + dec + " -ra_offset " + raoff + " -dec_offset " + decoff + " -start_date " + startdate + " -end_date " + enddate + " -device -device_name IO:O -device_type camera -device_spectral_region optical -device_filter " + filt + " -binning " + binning + " -exposure " + exposure + " s " + nexpose + " > " + docs_dir +"/ioo_rtml_document_" + now_date + ".rtml"


subprocess.call(string, shell=True)


string2 = bin_dir + "/nodetester.tcsh -url http://161.72.57.3:8080/org_estar_nodeagent/services/NodeAgent -username " + username + " -password " + password + " -cookie -axis -urn  urn:/node_agent -method handle_rtml -rtml \"" + docs_dir +"/ioo_rtml_document_" + now_date + ".rtml\""

print string2


subprocess.call(string2, shell=True)


