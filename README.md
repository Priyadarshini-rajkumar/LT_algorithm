This algorithm can be used to generate ranked lists of candidate galaxies within LIGO/Virgo localisation regions. 
The algorithm performs the following operations :
1.  Calculates the percentage of GW Skymap obsevable from a particular location  
2.  Load galaxy catalogues/ GW Bayestar Maps
3.	Filter the galaxy catalogs based on the GW error map provided by LIGO/ Virgo GW detectors 
4.	Produce a ranked list of candidates by combining GW error-map probability, observability, distance, and absolute magnitude data 
5.	Provide a "probability of success" estimate to help make a go/no-go decision on triggering the plan
6.	Create an observing schedule based on the time available for observation.
7.	Send the observing plan to LT’s fully autonomous Database to begin EM follow-up.


