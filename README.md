# Introduced species Checker

Use WRIMS + MarineRegions to determine if a given organism+location is invasive or not. The steps to checking are:
 - Provide sample location and AphiaID of detected organism
 - MarineRegion GIDs are retrieved from WRIMS. These are associated with the distribution of the AphiaID
 - For each MRGID a geom is pulled from the web, or from the internal cache, and inserted into a dataframe
 - The sample location is compared to the geoms to find whether the sample location is in, or near, a known region that hold introduced species
 - The human readable summary, as well as the dataframe (minus geom) is returned to the user

One point: up until now the species have been called "invasive". There is a technical difference between "Alien", "Introduced", and "Invasive" organisms. The returns from this system are known to be "introduced" and may, or may not, be "invasive". 


## Getting running in Docker
A docker-compose file and .env file are provided. Together they deploy a container that contains a REST API with Swagger documentation. Currently there are only a few REST endpoints and they'll look like this when running on a local system with the default .env file: 

 - localhost:8090 - The root of the API. Returns a simple message
 - localhost:8090/docs - The swagger documentation
 - localhost:8090/check - Takes an aphia_id, lon and lat and returns the summary and dataframe
 - localhost:8090/clear_cache - Clears the internal geom cache of the Python class. 

## Human Readable Summary

The json return for a request like http://localhost:8090/check/132818?lon=2.5&lat=51.5 would be:

```
{
  "summary": {
    "aphia_id": "132818",
    "sample location [WKT]": "POINT (2.5 51.5)",
    "sample location within aphia distibution": true,
    "sample location within <buffer> of aphia distribution": true,
    "buffer [deg]": 5,
    "species known to be introduced at sample location": true,
    "nearest introduced location": ["North Sea"],
    "distance [km] to nearest introduced location": 0,
    "distance [deg] to nearest introduced location": 0,
    "nearest introduced MRGID": [21912],
    "AphiaDistribution URL": "http://www.marinespecies.org/rest/AphiaDistributionsByAphiaID/132818"
  },
  "details": {
    "locality": {
      "1": "North Sea",
      "2": "Virginian"
    },
    "locationID": {
      "1": "http://marineregions.org/mrgid/21912",
      "2": "http://marineregions.org/mrgid/21853"
    },
    "typeStatus": {
      "1": "",
      "2": "holotype"
    },
    "establishmentMeans": {
      "1": "Introduced",
      "2": ""
    },
    "qualityStatus": {
      "1": "checked",
      "2": "checked"
    },
    "MRGID": {
      "1": 21912,
      "2": 21853
    },
    "contains_sample_point": {
      "1": true,
      "2": false
    },
    "distance_to_distribution": {
      "1": 0,
      "2": 71.25285240155577
    }
  }
}

  }
```
and the values mean:
```
summary": {
    "aphia_id": Aphia ID used in query,
    "sample location [WKT]": WKT of the sample point,
    "sample location within aphia distibution": Sample location is within any geom of the known organism distibution
    "sample location within <buffer> of aphia distribution": Sample location is within <buffer> degrees of any known organism distribution
    "buffer [deg]": The number of degrees to use for distance checks
    "species known to be introduced at sample location": Organism is known to be "introduced" to a marine region that the sample location is within. 
    "nearest introduced location": A list of names of the closest marine region where the species in "introduced". Multiple names indicate that the sample region is within multiple MarineRegions
    "distance [km] to nearest introduced location": An estimate of the shortest distance, in kilometers, from the sample location to the known introduced geom. 
    "distance [deg] to nearest introduced location": An estimate of the shortest distance in deg.
    "nearest introduced MRGID": A list of MRGIDs for the nearest introduced marine region/s.
    "AphiaDistribution URL": A URL to the distribution of the aphia_id from WORMS/WRIMS
  }
```

The "details" variables are:
```
"details": {
    "locality": Names of the known MarineRegion distributions
    "locationID": URL's for the known MarineRegion distributions
    "typeStatus": MarineRegion types
    "establishmentMeans": Organism location status. Could be Null, Introduced, Native, Native - Endemic, Native - Non-endemic, Origin uncertain, Origin unknown. The vast majority of returns will be Null or Introduced
    "qualityStatus": Status flag
    "MRGID": MarineRegion ID's
    "contains_sample_point": Does the MRGID contain the Sample Location? 
    "distance_to_distribution": Distance [deg] from Sample Location to geom
  }
}
```


## Some handy variables:
The following are useful for testing the system:

Aphia ID's

132818 - Sponge that is known to be introduced to few location around Europe

107451 - Crab that is known to be introduced to many locations around the world

126436 - Cod, not known to be introduced anywhere


## Python Code 

To check whether a species, identified using a WORMS Aphia ID, is known to be native or "introduced" in a specific location use the following snippet of code:

```python
from invasive_checker import invasive_checker

lon = 2.5
lat = 51.5
aphia_id = 132818

aphia_checker = invasive_checker.Aphia_Checker()
human_dict, pandas_dataframe = aphia_checker.check_aphia(lon, lat, aphia_id)
```
which should return a summary of the data found for the input data
```
print(human_dict)

{'aphia_id': 126436,
 'sample location [WKT]': 'POINT (2.5 51.5)',
 'sample location within aphia distibution': True,
 'sample location within <buffer> of aphia distribution': True,
 'buffer [deg]': 5,
 'species known to be introduced at sample location': False,
 'nearest introduced location': array([], dtype=object),
 'distance [km] to nearest introduced location': None,
 'distance [deg] to nearest introduced location': "No known 'introduced' locations",
 'nearest introduced MRGID': array([], dtype=int64),
 'AphiaDistribution URL': 'http://www.marinespecies.org/rest/AphiaDistributionsByAphiaID/126436'}
```
as well as a dataframe for the known distribution of the organism near to the sampling site:

```
pandas_dataframe.head()

locality	locationID	typeStatus	establishmentMeans	qualityStatus	MRGID	contains_sample_point	distance_to_distribution
0	European waters (ERMS scope)	http://marineregions.org/mrgid/7130	None	None	unreviewed	7130	True	0.000000
2	Belgian Exclusive Economic Zone	http://marineregions.org/mrgid/3293	None	None	unreviewed	3293	True	0.000000
3	Westerschelde	http://marineregions.org/mrgid/4752	None	None	unreviewed	4752	False	1.051264
6	British Isles	http://marineregions.org/mrgid/3140	None	None	unreviewed	3140	False	1.057823
8	Dutch Exclusive Economic Zone	http://marineregions.org/mrgid/5668	None	None	unreviewed	5668	False	0.342922

```
