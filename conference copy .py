#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math 
import pandas                as pd
import numpy                 as np
import astropy.units         as u
import warnings  
warnings.filterwarnings('ignore')
import healpy                as hp 
import numpy                 as np 
import matplotlib.pyplot     as plt
import astropy.coordinates
from astropy.time            import Time
from   astropy               import cosmology
from   datetime              import datetime
from scipy.stats             import norm
from scipy                   import integrate
from astropy.utils.data      import download_file
from astropy.cosmology       import default_cosmology, FlatLambdaCDM
from   astropy.table         import Table
from   astropy.time          import Time
from   astropy.coordinates   import SkyCoord,EarthLocation
from   astroplan             import (Observer, FixedTarget, AltitudeConstraint,AirmassConstraint,
                                     AtNightConstraint, MoonSeparationConstraint, is_observable,is_always_observable, ObservingBlock)
from   astroplan.constraints import TimeConstraint
from ligo.skymap.postprocess import find_greedy_credible_levels
from   astroplan.scheduling  import (Transitioner, Schedule, PriorityScheduler)
from   astroplan.plots       import (plot_airmass, plot_schedule_airmass)
#======================================================================================================
start_time          = Time("2019-06-29 21:00") # start time of LT observation 
end_time            = Time("2019-06-30 5:00")  # end time of LT observation
LT                  = Observer.at_site('lapalma')  # information about the site
LT_constraints      = [AltitudeConstraint(20 * u.deg), AtNightConstraint.twilight_nautical(),
                           MoonSeparationConstraint(5 * u.deg), AirmassConstraint(None)]     # constraits for observability
Dist                = 40  # estimated distance of the merger event; units in Mega parsec
Dist_err            = 20  # distance error in Mega parsecs
Event_ID            = 'GW170817'
observation_time    = 2 * 60  #mins
seeing              = 5.0   # units in arcsec
required_SNR        = 50.0  # used to calculate an estimate for exposure time
sky_brightness      = 2     
number_of_exposures = 2   
no_of_filters       = 1
no_of_exposures     = [1]  # no of exposures in each filter
GLADE_dir           = "/mnt/d/LT_algorithm/GLADE_2.3.txt"
bayestar_dir        = "/mnt/d/LT_algorithm/S190408an_bayestar.fits.gz,0"     

#=======================================================================================================
def gaussian(x, Dist, Dist_err):
    return np.exp(-np.power(x - Dist, 2.) / (2 * np.power(Dist_err, 2.)))

def distance_filter(Dist, Dist_err, data):
    no_dist         = data.loc[:,'dist']!='NaN'   # removing galaxies that do not have a distance
    data            = data[no_dist]               # altered dataframe with non-zero values for dist
    dist_max        = Dist + Dist_err      
    dist_min        = Dist - Dist_err
    dist_limit      = (data.loc[:,'dist'] < dist_max) & (data.loc[:,'dist'] > dist_min) # removing galaxies which donot satisfy GraceDB distance approx
    host_candidates = data[dist_limit]       
    score_list      = []                             # following steps reorder all list of galaxies according to the distance probability by GraceDB
    scoring_factor  = Dist_err / 10.0
    for i in host_candidates.loc[:,"dist"]:                 
            i       = Dist_err - abs(i-Dist)
            score   = i / scoring_factor
            score_list.append(score)
    host_candidates.loc[:,'dist_score'] = score_list
    return  host_candidates

def maximum_altitude(df):  # function to calculate maximum_altitude of all galaxies 
    lat           = math.radians(28.7624)  # latitude of the LT_observatory 
    max_Alt       = []
    for i in df.dec:
        i         = math.radians(i)  # dec in radians
        sin_alt   = math.sin(i)*math.sin(lat)+math.cos(i)*math.cos(lat)   # relation between latitude of LT, dec and alt
        alt       = math.asin(sin_alt)   # alt in radians
        alt       = math.degrees(alt) # alt in degrees
        max_Alt.append(alt)
    df.loc[:,"max_alt"] = max_Alt
    return df

def credibility(LT_targets, nside, credible_levels):
    credible_list          = []
    probability_levels     = []
    for i in LT_targets.index:
        ra                 = LT_targets.loc[i,['ra']]
        ra                 = int(ra) # functions deg2rad and ang2pix need interger values of phi
        dec                = LT_targets.loc[i,'dec']
        theta              = 0.5 * np.pi - np.deg2rad(dec)
        phi                = np.deg2rad(ra)
        ipix               = hp.ang2pix(nside, theta, phi)
        prob               = round(credible_levels[ipix],2)   # probabilty_levels of each galaxy 
        probability_levels.append(prob)
        credible           = prob <= 0.9
        credible_list.append(credible)
   
    return probability_levels, credible_list

def scoring(df):                                                                        
    score_list     = []
    max_value      = np.nanmax(df)
    min_value      = np.nanmin(df)
    diff           = max_value - min_value
    scoring_factor = diff / 10.0  
    for i in df:
        i         = i - min_value
        score     = i / scoring_factor
        score     = round(score,1)
        score_list.append(score)
    scale          = [10]*len(score_list)
    score_list     = np.array(scale) - np.array(score_list)
    return score_list

def final_score(df):   # function to calucated a weighted score to rank the galaxies
    BMAG_weight        =  50          
    prob_region_weight =  80
    Dist_weight        =  80
    total              = BMAG_weight + prob_region_weight + Dist_weight
    final_score_list   = (df.BMAG_score * BMAG_weight) +  (df.prob_region_score * prob_region_weight) + (df.dist_score * Dist_weight) 
    final_score_list   = final_score_list / total
    return final_score_list
 
def BNS_magnitude(Dist, Dist_err):
    Dist_range   = np.arange(Dist-Dist_err, Dist + Dist_err, 0.1) 
    Dist_range   = Dist_range * 1e6 # parces#
    absolute_mag  = -15.8   # expected magnitude of the optical counterpart of BNS mergers
    mag_range    = []
    for dist in Dist_range:
        mag    = absolute_mag+ 5*math.log10(dist)-5      # apparent magnitude to be used to approximate for exposure time
        mag_range.append(mag)
    mag          = round(np.average(mag_range),2)
    return mag    


def exposure_time(mag,snr,seeing,sky_tonight):
    zp             = 25.14     # SDSS_g; unit in mags
    skybr          = 21.7   # SDSS_g; unit in mags
    skyoff         = 1.0   # SDSS_g; correction for moon brightness
    pixscale       = 0.15    # IO:O
    darkcurrent    = 0    # IO:O
    readnoise      = 10   # IO:O
    binn           = 2
    skymag         = skybr - (skyoff * sky_tonight)   # stepping through sky brightnesses (in mags)
    areaofdisk     = np.pi * seeing * seeing  
    numberofpixels = areaofdisk / (pixscale * pixscale * binn* binn)
    starphotons    = 10**((zp-mag)/2.5) 
    skyphotons     = 10**((zp-skymag)/2.5)                              
    a              = starphotons * starphotons
    b              = -snr * snr * (starphotons + skyphotons + darkcurrent)  # solving the quadratic equation obtained by re-arranging the SNR estimation eq.                     
    c              = -snr * snr * numberofpixels * readnoise * readnoise                              
    t_e1           = (-b + math.sqrt(b*b - 4*a*c))/(2*a)   # solution 1 
    t_e2           = (-b - math.sqrt(b*b - 4*a*c))/(2*a) # solution 2
    exposure_time  = max(t_e1,t_e2)   # maximum of the two solutions
    if exposure_time < 30:
        exposure_time = 30
    return exposure_time

def total_observation_time(exposure_time, n, no_of_filters, no_of_exposures): # n =  no. of targets
    readout_time        = 18.5   # for 2X2 on-chip binning 
    aquisition_time     = 60     # time taken to slew the telescope on target  
    filter_change_time  = 15
    all_obs = 0    # to
    for i in np.arange(no_of_filters):
            one_obs    = (aquisition_time + filter_change_time + no_of_exposures[i] * (exposure_time + readout_time)) # observation time for one object
            total_obs  = n * one_obs
            all_obs    = all_obs + total_obs  # adding each filter observation time to find total time
    one_obs = round(one_obs/60)
    total_obs = round(all_obs/60)   # from seconds to minutes         
    return one_obs, total_obs

def schechter():
    binwidth    = 0.2
    H0          = 70
    h_50        = H0 / 50
    phi_star    = 0.0002 * (h_50 ** 3) #Mpc-3
    L_star      = 10.38
    alpha       = -1.25 
    L           = np.linspace(6, 11, 1e6)
    schechter   = (1/binwidth)*phi_star * ((10**L / 10**L_star)** (alpha))* np.exp (- 10**L / 10**L_star) * ((10**(L+0.1))+(10**(L-0.1)))/10**L_star
    plt.loglog(10**L,schechter,'r')
    plt.xlabel('Lumninosity')
    plt.title('Schecter Lumninosity function')
    plt.grid(True)
    plt.ylabel('Number density/ Mpc3')
    plt.show()
    objects_per_Mpc3 = integrate.simps(schechter, L)
    return objects_per_Mpc3


#=====================================================================================================
GLADE               = pd.read_csv(GLADE_dir, delimiter=' ',header = None)
GLADE.columns       = ['PGC','name', 'HyperLEDA name','2MASS name','SDSS-DR12 name',
                       'flag1','ra','dec','dist','dist_err','z','Bmag', 'Bmag_err','BMAG',
                       'Jmag','Jmag_err','Hmag','Hmag_err','Kmag','Kmag_err','flag2','flag3']  
GLADE            = GLADE.loc[:,['name','ra','dec','dist','BMAG']]
#====================================================================================================
prob, distmu, distsigma, distnorm = hp.read_map(bayestar_dir, field = [0, 1, 2, 3], verbose = False)
credible_levels       = find_greedy_credible_levels(prob) # obtaining the probability for each pixel if it is in or outside the error region
npix                   = len(prob)
nside                  = hp.npix2nside(npix)

time2                    = Time("2018-06-26")
LT2                      = astropy.coordinates.EarthLocation( lat = 28.7583332*u.deg, lon = - 17.8799*u.deg, height=2327*u.m) # information about the site   
frame                   = astropy.coordinates.AltAz(obstime = time2, location = LT2)
theta, phi              = hp.pix2ang(nside, np.arange(npix))
radecs                  = astropy.coordinates.SkyCoord(ra = phi*u.rad, dec = (0.5*np.pi - theta)*u.rad)
altaz                   = radecs.transform_to(frame) 
sun_altaz               = astropy.coordinates.get_sun(time2).transform_to(altaz)    
observing_probability   = round(prob[(sun_altaz.alt <= -18*u.deg) & (altaz.secz <= 2.5)].sum()*100)  # least 18 degrees below the horizon and that the airmass (secant of zenith angle approximation) is at most 2.5.
hp.mollview(prob)
print ('Percentage of GW skymap observable', observing_probability,'%')
plt.show()
#=====================================================================================================
#distance_probability plotting
r     = np.linspace(0, 200)
dp_dr = gaussian(r, Dist, Dist_err)
plt.plot(r, dp_dr)
plt.xlabel('distance (Mpc)')
plt.ylabel('prob Mpc$ˆ{-1}$')
plt.show()
#======================================================================================================
host_candidates                            = distance_filter(Dist, Dist_err, GLADE)
probability_levels, credible_list          = credibility(host_candidates, nside,credible_levels)
host_candidates.loc[:,'probabilty_region'] = probability_levels
host_candidates                            = host_candidates[credible_list]
df                                         = host_candidates.loc[:,["name","ra","dec"]] 
targets_astropy                            = Table.from_pandas(df) # converting the dataframe into an astropy table
time_range                              = Time([start_time, end_time])   # time of observation
targets                                 = [FixedTarget(coord = SkyCoord(ra = ra*u.deg, dec=dec*u.deg), name = name) 
                                            for name, ra, dec in targets_astropy]      
LT                                      = Observer.at_site('lapalma')  # information about the site
observable                              = is_observable(LT_constraints,LT,targets, time_range=time_range)  # is observable atleast for a short duration for the given time range
LT_targets                              = host_candidates[observable]
no_BMAG                                 = LT_targets.loc[:,'BMAG']   != 'NaN'                               # removing galaxies that do not have absolute Blue mag
LT_targets                              = LT_targets[no_BMAG]
BMAG_score                              = scoring(LT_targets.BMAG)   # a list of score values for targets based on Blue Magnitude
LT_targets.loc[:,"BMAG_score"]          = BMAG_score  # adding Bmag_score as a column to the datframe
Prob_region_score                       = scoring(LT_targets.probabilty_region)
LT_targets.loc[:, 'prob_region_score']  = Prob_region_score
final_score_list                        = final_score(LT_targets.loc[:, ["BMAG_score","prob_region_score",'dist_score']])
LT_targets.loc[:, "Final_score"]        = final_score_list
LT_targets                              = LT_targets.sort_values('Final_score', ascending = False)
BNS_mag                                 = BNS_magnitude(Dist, Dist_err)
print ('Expected BNS merger apparent magnitude: ', BNS_mag)
exposure_time                           = exposure_time(BNS_mag, required_SNR, seeing, sky_brightness)
print ('Exposure time:',exposure_time,'s')
n                         = len(LT_targets)
one_filter_obs, total_obs = total_observation_time(exposure_time, n, no_of_filters, no_of_exposures)
one_target_obs            = one_filter_obs * no_of_filters
print ('Time to complete observation of one target:', one_target_obs,'m')
print  ('Total approximated time of observation:', total_obs,'m')
#======================================================================================================
objects_per_Mpc3 = schechter()
print ('No of galaxies per Mpc^3:', objects_per_Mpc3)
GW_sq_deg       = np.sum(credible_levels <= 0.9) * hp.nside2pixarea(nside, degrees=True)  #compute the 90% credible area by counting the number of pixels inside the 90% credible region and multiplying by the area per pixel
print ('Area of event', GW_sq_deg)
whole_sphere_deg= 41253
Dist_range      = np.linspace(Dist - Dist_err, Dist + Dist_err, 1e3)
area_range      = 4 * np.pi * (GW_sq_deg/whole_sphere_deg) * (Dist_range **2)
GW_vol_Mpc3      = integrate.simps(area_range, Dist_range)
print ('GW Probability volume:', GW_vol_Mpc3,'Mpc3')
no_of_obj       = int(GW_vol_Mpc3*objects_per_Mpc3)
print ('Luminous Galaxies in GW Probability Volume:', no_of_obj)
#======================================================================================================

prob         = [(i) for i in dp_dr if i>1e-5]
id_1          = [i for i in np.arange(len(dp_dr)) if prob[0] == dp_dr[i]]
id_2        = [i for i in np.arange(len(dp_dr)) if prob[-1] == dp_dr[i]]
s =r[id_1[0]]
e = r[id_2[0]]
d = e-s
print (s,e,d)
Dist         = np.arange(s,e,d/len(prob))  # Distance distribution
RCF_d        = 0.8 - 0.00147 * Dist
plt.plot(Dist[0:len(RCF_d)], RCF_d,'r')
plt.title('Catalog Completeness')
plt.xlabel('Distance (Mpc)')
plt.grid()
plt.ylabel('Redshift Correctness Factor (RCF)')                # redshift correctness factor
plt.show()

def weighted_RCF(RCF_d,prob): # crossmatching 
    sum_of_weights = sum(prob)
    s = 0             
    for i in np.arange(len(RCF_d)):
        s            = RCF_d[i] * prob[i] + s  # sum of Number X weighting factor
    weighted_average = round(s / sum_of_weights,4)  #weighted average
    return weighted_average
weighted_RCF_value  = weighted_RCF(RCF_d, prob)  # completeness of the catalog
print('Catalog Completeness:', weighted_RCF_value*100,'%')
Bmag                = LT_targets[0:100]['BMAG'].sum()
total_Bmag          = host_candidates['BMAG'].sum()            
probability         = Bmag / total_Bmag                   #completeness of the observation
completeness        = round(weighted_RCF_value * probability * 100,2)
print ('Completeness of the observation: ',completeness,'%')
print ('DIS', Dist)


# In[ ]:




