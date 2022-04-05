# Introduced species Checker

Use WRIMS + MarineRegions to determine if a given organism+location is invasive or not. 

## Steps 

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
