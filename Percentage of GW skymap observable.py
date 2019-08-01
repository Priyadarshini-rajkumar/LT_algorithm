#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import astropy.coordinates
from astropy.time import Time
import warnings  
warnings.filterwarnings('ignore')
import astropy.units as u
import healpy as hp
import numpy as np



bayestar_dir            = "/data1/extprajk/GW170817/LALInference_v2.fits.gz"
prob                    = hp.read_map(bayestar_dir, field = [0], verbose = False)   # obtaining the probability for each pixel if it is in or outside the error region
npix                    = len(prob)
nside                   = hp.npix2nside(npix)
time                    = Time("2018-06-26 23:00")
#or
#time                    = time.now()
LT                      = astropy.coordinates.EarthLocation( lat = 28.7583332*u.deg, lon = - 17.8799*u.deg, height=2327*u.m) # information about the site   
frame                   = astropy.coordinates.AltAz(obstime = time, location = LT)
theta, phi              = hp.pix2ang(nside, np.arange(npix))
radecs                  = astropy.coordinates.SkyCoord(ra = phi*u.rad, dec = (0.5*np.pi - theta)*u.rad)
altaz                   = radecs.transform_to(frame) 
sun_altaz               = astropy.coordinates.get_sun(time).transform_to(altaz)    
observing_probability   = round(prob[(sun_altaz.alt <= -18*u.deg) & (altaz.secz <= 2.5)].sum()*100)  # least 18 degrees below the horizon and that the airmass (secant of zenith angle approximation) is at most 2.5.
hp.mollview(prob)
print ('Percentage of GW skymap observable', observing_probability,'%')
    
    
 


# In[ ]:




