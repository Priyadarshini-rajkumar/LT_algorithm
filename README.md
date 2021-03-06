**LIVERPPOOL TELESCOPE OPTIMIZATION ALGORITHM**

This algorithm can be used to generate ranked lists of candidate galaxies within LIGO/Virgo localisation regions. 

The algorithm performs the following operations :
* Calculates the percentage of GW Skymap obsevable from a particular location (Check GW_skymap_observability.py)
* Loads galaxy catalogues/ GW Bayestar Maps
* Filters the galaxy catalogs based on the GW error map provided by LIGO/ Virgo GW detectors 
* Produces a ranked list of candidates by combining GW error-map probability, observability, distance, and absolute magnitude data 
* Provides a "probability of success" estimate to help make a go/no-go decision on triggering the plan
* Creates an observing schedule based on the time available for observation.
* Sends the observing plan to LT’s fully autonomous Database to begin EM follow-up.


