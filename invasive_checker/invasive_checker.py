import requests
import pandas as pd
import geopandas as gpd
import pyproj
import re 
import rdflib
import os
import json
import functools
import logging
import shapely
import math
import warnings
from shapely.errors import ShapelyDeprecationWarning
from shapely.geometry import Point, Polygon
from shapely.ops import nearest_points
import datetime

from rdflib import Graph 
from rdflib.namespace import RDF 

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

log = logging.getLogger('invasive_checker') 

def derive_status(this_aphia_df):
    '''
    Provide some human readable results...
    There might be some issues with this:
      - What happens when a sample isn't enclosed by a marineregion?
      - What happens if the only marineregion is large? Something like "Atlantic Ocean"
      - How do we say that a sample is unrecorded at a marineregion? How do we find "missing" marineregions?
    '''
    log.debug(f'Preparing results...') 
    if this_aphia_df is None:
        status = ['Unrecorded']
        mrgids = ['None']
    else:
            
        status = list(this_aphia_df.establishmentMeans.values)
        mrgids =list(this_aphia_df.MRGID.values)

    derived_status = {  'Status':status,
                        'Within':mrgids,
                        } 
        
    return derived_status 


def check_aphia(lon, lat, id, source='worms'):
    '''
    Check if the aphia is invasive, native or unknown. Return MRGID's that are within <buffer> degrees 
    of the sample location.
    
    CRS for everything is WGS84.
    
    aphia_df is a dataframe of MRGID's and status' that are associated with the aphia_id. Retrieved from WRIMS
    
    geom_df is a dataframe of MRGID geoms retrieved from MR. 
    geom_store is a dict of geom_df
    '''
    
    log.debug(f'Received request for {id} for location {lon}/{lat} ')
    
    # Get Aphia status, geoms associated with the aphia and whether the sample location is 
    # inside those geoms or not. 
    # =============
    if source == 'worms':
        aphia_id = id
        this_aphia_df = get_aphia_status(aphia_id)
    else:
        aphia_id = get_external_status(id, source)
        this_aphia_df = get_aphia_status(aphia_id)

    if this_aphia_df is None:
        status_dict = {'aphia_id': aphia_id,
                       'Error': 'No distribution found for this Aphia_ID'}
        return status_dict, None
 
     
    # Finds MRGIDs that intersect with the sample location 
    log.debug(f'  -Finding MarineRegions that intersect with sample...')
    sample_mr_response = get_mrgid_from_latlon(lat,lon)
    sample_mrgids = list(set([d.get('MRGID') for d in sample_mr_response]))

    # Find intersecting MRGIDs that are also in the WRIMS response
    this_aphia_df = this_aphia_df[this_aphia_df['MRGID'].isin(sample_mrgids)]
    #-----------------------

    # this_aphia_gdf['distance_to_distribution'] = this_aphia_gdf.distance(sample_point)
    this_aphia_df.establishmentMeans.replace('Alien','Introduced',inplace=True)
    this_aphia_df.establishmentMeans.fillna('Recorded',inplace=True)
    
    # =============
    # Apply some human logic to determine whether the above results are of interest or not
    # =============
    log.debug(f'  -Applying logic...')
    status_dict = derive_status(this_aphia_df)

    # =============
    # Prep for returning to user
    #   - Drop geom column since it can get very big 
    #   - drop non-useful or confusing columns
    # =============
    # this_aphia_gdf = this_aphia_gdf[this_aphia_gdf['distance_to_distribution'] < buffer]
    invasive_df = this_aphia_df.drop(['decimalLongitude', 
                                    'decimalLatitude', 
                                    'higherGeography',
                                    'higherGeographyID'],axis=1)
    
    # =============
    log.debug(f'  -Done:')
    log.debug(status_dict)
    return status_dict, invasive_df
   
def get_external_status( external_id, id_source):
    '''
    Get the APHIA ID from an externalID. See
    https://marinespecies.org/rest/AphiaRecordByExternalID/

    algaebase: Algaebase species ID
    bold: Barcode of Life Database (BOLD) TaxID
    dyntaxa: Dyntaxa ID
    fishbase: FishBase species ID
    iucn: IUCN Red List Identifier
    lsid: Life Science Identifier
    ncbi: NCBI Taxonomy ID (Genbank)
    tsn: ITIS Taxonomic Serial Number
    gisd: Global Invasive Species Database
    '''
    id_sources = {'algaebase': 'Algaebase species ID',
                    'bold': 'Barcode of Life Database (BOLD) TaxID',
                    'dyntaxa': 'Dyntaxa ID',
                    'fishbase': 'FishBase species ID',
                    'iucn': 'IUCN Red List Identifier',
                    'lsid': 'Life Science Identifier',
                    'ncbi': 'NCBI Taxonomy ID (Genbank)',
                    'tsn': 'ITIS Taxonomic Serial Number',
                    'gisd': 'Global Invasive Species Database'}

    if id_source not in id_sources:
        log.warning(f'Unknown external source: {id_source}.')
        return None 
    try:
        aphia_url = f'https://marinespecies.org/rest/AphiaRecordByExternalID/{external_id}?type={id_source}'
        aphia_return =  requester(aphia_url)
        if aphia_return is not None:
            aphia_id = aphia_return.json()['AphiaID']
            return aphia_id
        else:
            return None
    except Exception as err:
        log.warning(f'Error retrieving distribution for {id_source}: {external_id}')
        log.warning(err)
        return None 

def get_aphia_status(aphia_id):
    '''
    Get the MRGIDs and invasive status for the aphia_id specified. Return dataframe with MRGID's of 
    known distribution and the the native/alien status of the MRGID/APHIA pair.
    Good test values are aphiaID = 107451 (chinese mitten crab, invasive)
    '''
    wrms_distribution = f'http://www.marinespecies.org/rest/AphiaDistributionsByAphiaID/{aphia_id}'
    try:
        req_return = requester(wrms_distribution)
        if req_return is not None:
            wrms_dist = req_return.json()
        else:
            return None

        if len(wrms_dist) == 0 :
            log.warning(f'No distribution for Aphia {aphia_id}')
            return None
        else:
            wrms_dist_df = pd.DataFrame(wrms_dist) 
            # Filter the data based off of comments on confluence:
            # =============
            wrms_dist_df = wrms_dist_df.drop_duplicates()
            wrms_dist_df = wrms_dist_df[wrms_dist_df.recordStatus == 'valid']
            # =============
            wrms_dist_df[['root','MRGID']] = wrms_dist_df.locationID.str.rsplit('/',1,expand=True)
            wrms_dist_df['MRGID'] = wrms_dist_df['MRGID'].astype(int)
            wrms_dist_df = wrms_dist_df.drop(['root'],axis=1)
            return wrms_dist_df
    except Exception as err:
        log.warning(f'Error retrieving distribution for aphia {aphia_id}')
        log.warning(err)
        return None

@functools.cache
def requester(url):
    '''
    Do a safe request and return the result.
    '''
    log.debug('    -Doing URL request: {0}'.format(url))
    reply = requests.get(url)
    if reply.status_code == 200:
        try:
            r = reply
        except Exception as error:
            log.error(error)
            r = None
        return r
    elif reply.status_code == 204:
        log.warning('No Content for {0}...'.format(url))
        return None
    else: 
        # Something not right with the request...
        log.warning(reply)
        log.warning(reply.text)
        return None

def get_aphia_from_lineage(tax_string, sep = ';'):
    '''
    Given a taxon lineage string (Eukaryota;Chordata;Ascidiacea;Enterogona;Ascidiidae;Ascidiella;Ascidiella scabra)
    get the aphia_id for the lowest level.
    '''
    tax_lineage = tax_string.split(sep)
    log.debug('Checking taxon: {0}'.format(tax_lineage))
    req_return = None
    while req_return is None:
        try:
            req_return = get_aphia_from_taxname(tax_lineage[-1])
            tax_lineage.pop(-1)
        except IndexError:
            log.warning('Reached end of taxon lineage without success...')
            break
    return req_return

def get_aphia_from_taxname(taxa_name):
    '''
    Given a taxon name string, get the aphia_id/s that are associated with it. 

    https://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]=<some name>&marine_only=true
    '''
    taxamatch_url = f'https://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]={taxa_name}&marine_only=true'
    try:
        req_return = requester(taxamatch_url)

        if (req_return is None):
            log.warning(f'No AphiaID for taxname {taxa_name} found...')
            return None
        elif (req_return.status_code == 204):
            log.warning(f'No AphiaID for taxname {taxa_name} found...')
            return None
        elif req_return.status_code == 200:
            log.debug(f'Returns: {req_return}') 
            req_return = req_return.json()[0][0]
            return req_return
        else: 
            log.warning('Not sure how I got here...')
            log.warning(req_return)
            return None
    except Exception as err:
        log.warning(f'Error retrieving aphia_id for sci-name: {taxamatch_url}')
        log.warning(err)
        return None


def get_mrgid_from_latlon(lat,lon):
    '''
    Given a the location of a sample, find the Marineregions that intersect with it. 
    https://www.marineregions.org/rest/getGazetteerRecordsByLatLong.json/{lat}{lon}/?offset=0
    '''
    mr_url = f'https://www.marineregions.org/rest/getGazetteerRecordsByLatLong.json/{lat}/{lon}/?offset=0'
    try:
        req_return = requester(mr_url)
        if (req_return is None):
            log.warning(f'No MarineRegions for {lat}/{lon} found...')
            return None
        elif (req_return.status_code == 204):
            log.warning(f'No MarineRegions for {lat}/{lon} found...')
            return None
        elif req_return.status_code == 200:
            log.debug(f'Returns: {req_return}') 
            req_return = req_return.json() 
            return req_return
        else: 
            log.warning('Not sure how I got here...')
            log.warning(req_return)
            return None
    except Exception as err:
        log.warning(f'Error retrieving MRGIDs for location: {lat}/{lon}')
        log.warning(err)
        return None