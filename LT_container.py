import math 
import astroplan
import astropy 
import pandas as pd
import numpy as np
import healpy as hp  
import matplotlib.pyplot as plt
from ligo.skymap.postprocess import find_greedy_credible_levels

#=======================================================================================================
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

def weighted_RCF(RCF_d,prob): # crossmatching 
    sum_of_weights = sum(prob)
    s = 0             
    for i in np.arange(len(RCF_d)):
        s            = RCF_d[i] * prob[i] + s  # sum of Number X weighting factor
    weighted_average = round(s / sum_of_weights,4)  #weighted average
    return weighted_average
