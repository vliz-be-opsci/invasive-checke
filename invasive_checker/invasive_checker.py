import requests
import pandas as pd
import geopandas as gpd
import pyproj
import re 
import rdflib
import os
import json
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


class Aphia_Checker : 
    def __init__(self):
        log.debug('Inititialising...')
        self.geom_store = {} 
        self.cache_persist_period = int(os.getenv('CACHE_PERIOD',default=7))
    
    def check_aphia(self,lon,lat,id, buffer = 5, source='worms'):
        '''
        Check if the aphia is invasive, native or unknown. Return MRGID's that are within <buffer> degrees 
        of the sample location.
        
        CRS for everything is WGS84.
        
        aphia_df is a dataframe of MRGID's and status' that are associated with the aphia_id. Retrieved from WRIMS
        
        geom_df is a dataframe of MRGID geoms retrieved from MR. 
        geom_store is a dict of geom_df
        '''
        
        log.debug(f'Received request for {id} for location {lon}/{lat}:buffer={buffer}')
        sample_point = Point(lon,lat)

        # Get Aphia status, geoms associated with the aphia and whether the sample location is 
        # inside those geoms or not. 
        # =============
        if source == 'worms':
            aphia_id = id
            this_aphia_df, this_aphia_dist_url = self.get_aphia_status(aphia_id)
        else:
            aphia_id = self.get_external_status(id, source)
            this_aphia_df, this_aphia_dist_url = self.get_aphia_status(aphia_id)
        if this_aphia_df is None:
            human_readable_results = {'aphia_id': aphia_id,
                                  'Error': 'No distribution found for this Aphia_ID'}
            return human_readable_results, None
        
        log.debug(f'Fetching geom...')
        this_aphia_df['geom'] = this_aphia_df.apply(lambda x: self.get_rdf_geom(x['MRGID']), axis=1)
        log.debug(f'Doing GIS work...')
        this_aphia_gdf = gpd.GeoDataFrame(this_aphia_df, geometry='geom')        
        this_aphia_gdf['contains_sample_point'] = this_aphia_gdf.contains(sample_point)
        this_aphia_gdf['distance_to_distribution'] = this_aphia_gdf.distance(sample_point)
        
        # =============
        # Apply some human logic to determine whether the above results are of interest or not
        #   - replace "Alien" with "introduced"
        # =============
        this_aphia_gdf.establishmentMeans.replace('Alien','Introduced',inplace=True)
        introduced_rows = this_aphia_gdf[this_aphia_gdf['establishmentMeans'] == 'Introduced']
        sample_within_MR = this_aphia_gdf.contains_sample_point.any()
        sample_closeto_MR = (this_aphia_gdf.distance_to_distribution < buffer).any()
        MR_introduced = len(introduced_rows[(introduced_rows['contains_sample_point'] == True)].index) > 0
        MR_dist_introduced = min(introduced_rows.distance_to_distribution, default="No known 'introduced' locations")
        closest_MR = introduced_rows[introduced_rows.distance_to_distribution == introduced_rows.distance_to_distribution.min()]
        closest_distance = self.calc_distance(closest_MR, sample_point) 

        log.debug(f'Preparing results...')
        human_readable_results = {'aphia_id': aphia_id,
                                  'sample location [WKT]': sample_point.wkt,
                                  'sample location within aphia distibution': bool(sample_within_MR),
                                  'sample location within <buffer> of aphia distribution': bool(sample_closeto_MR),
                                  'buffer [deg]': buffer,  
                                  'species known to be introduced at sample location': bool(MR_introduced),
                                  'nearest introduced location': closest_MR.locality.values.tolist(),
                                  'distance [km] to nearest introduced location': closest_distance,
                                  'distance [deg] to nearest introduced location':MR_dist_introduced,
                                  'nearest introduced MRGID': closest_MR.MRGID.values.tolist(),
                                  'AphiaDistribution URL':this_aphia_dist_url}
        
        # =============
        # Prep for returning to user
        #   - Drop geom column since it can get very big 
        #   - drop non-useful or confusing columns
        # =============
        # this_aphia_gdf = this_aphia_gdf[this_aphia_gdf['distance_to_distribution'] < buffer]
        invasive_df = pd.DataFrame(this_aphia_gdf.drop(['geom'],axis=1))  
        invasive_df = invasive_df.drop(['decimalLongitude', 
                                        'decimalLatitude',
                                        'recordStatus',
                                        'higherGeography',
                                        'higherGeographyID'],axis=1)
        
        # =============
        log.debug(f'Done:')
        log.debug(human_readable_results)
        return human_readable_results, invasive_df
    
    def get_rdf_geom(self, mrgid):
        '''
        Take single MRGID and return a single geom that's a combination of all geom's associated with that MRGID. 
        Use some simple caching to speed things up. 
        '''
        if mrgid in self.geom_store:
            data_store = self.geom_store[mrgid]
            single_geom = data_store['geom']
            cache_datetime = data_store['store_date']
            today = datetime.date.today()
            if (today - cache_datetime).days > self.cache_persist_period:
                log.info(f'MRGID {mrgid} geom is old. Throwing it out...')
                data_store.pop[mrgid]
            return single_geom
        else:
            mr_geom_request = f'https://marineregions.org/rest/getGazetteerGeometries.jsonld/{mrgid}/'
            rdf = self.requester(mr_geom_request)
            g = Graph()
            g.parse(rdf.content, format='json-ld')
            wkt_pred = rdflib.term.URIRef('http://www.opengis.net/ont/geosparql#asWKT')

            geoms = []
            for s, p, o in g: 
                if p == wkt_pred:
                    try:
                        bad_wkt = o.n3()
                        xx = re.search('^.*\>\s(.*)\".*$', bad_wkt)
                        good_wkt = xx[1]
                        geoms.append(shapely.wkt.loads(good_wkt))
                    except:
                        pass

            num_geoms = len(geoms) 
            log.debug(f'  -combining {num_geoms} geoms for {mrgid}')
            single_geom = shapely.ops.unary_union(geoms)


            self.geom_store[mrgid] = {'geom':single_geom,
                                      'store_date': datetime.date.today()}
            return single_geom


    def get_external_status(self, external_id, id_source):
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
            log.warning(f'Unkown external source: {id_source}.')
            return None 
        try:
            aphia_url = f'https://marinespecies.org/rest/AphiaRecordByExternalID/{external_id}?type={id_source}'
            aphia_return =  self.requester(aphia_url).json()
            if aphia_return is not None:
                aphia_id = aphia_return['AphiaID']
                return aphia_id
            else:
                return None
        except Exception as err:
            log.warning(f'Error retriving distribution for {id_source}: {external_id}')
            log.warning(err)
            return None 

    def get_aphia_status(self, aphia_id):
        '''
        Get the MRGIDs and invasive status for the aphia_id specified. Return dataframe with MRGID's of 
        known distribution and the the native/alien status of the MRGID/APHIA pair
        '''
        wrms_distribution = f'http://www.marinespecies.org/rest/AphiaDistributionsByAphiaID/{aphia_id}'
        try:
            wrms_dist = self.requester(wrms_distribution).json()

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
                return wrms_dist_df, wrms_distribution
        except Exception as err:
            log.warning(f'Error retriving distribution for aphia {aphia_id}')
            log.warning(err)
            return None,None
    
    def calc_distance(self, df_row, point):
        '''
        Calculate the nearest distance in Km from the point to the polygon.
        Round it off to imply that it's not super accurate. 
        '''

        if len(df_row.index) > 0:  
            polygon = df_row['geom'].values[0]
            # The points are returned in the same order as the input geometries:
            p1, p2 = nearest_points(polygon, point)

            # Haversine formula
            dLat = math.radians(p2.y) - math.radians(p1.y)
            dLon = math.radians(p2.x) - math.radians(p1.x)
            a = math.sin(dLat/2) * math.sin(dLat/2) + math.cos(math.radians(p1.y)) * math.cos(math.radians(p2.y)) * math.sin(dLon/2) * math.sin(dLon/2)
            distance = 6371 * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
            distance = round(distance, -1)
            return distance
        else: 
            return None
    
    def requester(self, url):
        '''
        Do a safe request and return the result.
        '''
        reply = requests.get(url)
        if reply.status_code == 200:
            try:
                r = reply
            except Exception as error:
                log.error(error)
                raise Exception("No known distribution")
                r = []
            return r
        elif reply.status_code == 204:
            log.warning('No items found...')
            return None
        else: 
            # Something not right with the request...
            log.warning(reply)
            log.warning(reply.text)
            return []

    
    def clear_cache(self):
        '''
        Dump the stored dicts. Forces next call to be retrieved from web.
        '''
        cache_len = len(self.geom_store)
        xx = f'Clearing cache: dumping {cache_len} geoms...'
        log.info(xx)
        self.geom_store = {} 
        return {'message':xx}
    
    def get_aphia_from_taxname(self, taxa_name):
        '''
        Given a taxon name string, get the aphia_id/s that are associated with it. 

        https://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]=<some name>&marine_only=true
        '''
        taxamatch_url = f'https://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]={taxa_name}&marine_only=true'
        try:
            req_return = self.requester(taxamatch_url)

            if req_return.status_code == 204 :
                log.warning(f'No aphia for Aphia {taxa_name} found...')
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
            log.warning(f'Error retriving aphia_id for sci-name {taxa_name}')
            log.warning(err)
            return None
