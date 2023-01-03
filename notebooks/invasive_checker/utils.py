import logging
import pandas as pd
import functools
import requests
import time

log = logging.getLogger('utils') 

def get_sample_location_df(AccessionIDs, meta_df):
    '''
    Return a dataframe of AccessionID metadata. 
    '''
    row = {}
    for sample_id in AccessionIDs: 
        log.info('Getting metadata for {0}'.format(sample_id))
        meta_row = meta_df[meta_df.isin([sample_id]).any(axis=1)]
        meta_cols = meta_row[meta_row.isin([sample_id])].dropna(axis='columns')
        gene = meta_cols.columns.values[0]
        if gene.startswith('negativeControl_gene'):
            negative_control = True
        else:
            negative_control = False

        if negative_control:
            longitude = None
            latitude = None
            material_id = None
        else: 
            log.debug(meta_row)
            if meta_row.shape[0] > 1:
                log.debug('   -Possible duplicate gene ID...')
                longitude = None
                latitude = None
                material_id = None
            elif meta_row.shape[0] == 0:
                log.warning('   -Empty row?...')
                longitude = None
                latitude = None
                material_id = None
            else:
                longitude = meta_row.loc[meta_row.index].longitude.values[0]
                latitude = meta_row.loc[meta_row.index].latitude.values[0]
                material_id = meta_row.loc[meta_row.index].MaterialSample_ID.values[0] 

        row[sample_id] = {'AccessionNumber':sample_id,
                           'isNegativeControlGene':negative_control,
                           'GeneType':gene,
                           'sampleLatitude':latitude,
                           'sampleLongitude':longitude,
                           'MaterialSample_ID': material_id,
                           }
    my_df = pd.DataFrame(row).T
    return my_df

@functools.cache
def requester(url):
    '''
    Do a safe request and return the result.
    ''' 
    reply = requests.get(url)
    if reply.status_code == 200:
        try:
            r = reply
        except Exception as error:
            print(error)
            raise Exception("No known distribution")
            r = []
        return r
    
    elif reply.status_code == 429:
        print('Going too fast!')
        time.sleep(2)
        reply = requests.get(url)
        try:
            r = reply
        except Exception as error:
            print(error)
            raise Exception("No known distribution")
            r = []
        return r
    
    elif reply.status_code == 204:
        print('No items found...')
        return None
    
    else: 
        # Something not right with the request...
        print(reply)
        print(reply.text)
        return []
    
def get_mrgids(lon, lat):
        '''
        Get the MRGIDs for regions containing lon/lat
        
        ''' 
        mr_url = f'https://www.marineregions.org/rest/getGazetteerRecordsByLatLong.json/{lat}/{lon}/'

        try:
            mrgids = requester(mr_url).json()

            if len(mrgids) == 0 :
                print(f'No distribution for lat/lon {lon}, {lat}')
                return None
            else:
                wrms_dist_df = pd.DataFrame(mrgids) 
                # =============
                wrms_dist_df = wrms_dist_df.drop_duplicates() 
                # ============= 
                
                return  wrms_dist_df
        except Exception as err:
            print(f'Error retrieving distribution for lat/lon {lon}, {lat}')
            print(err)
            return [] 
            
def get_aphia_status(aphia_id):
        '''
        Get the MRGIDs and invasive status for the aphia_id specified. Return dataframe with MRGID's of 
        known distribution and the the native/alien status of the MRGID/APHIA pair
        '''
        wrms_distribution = f'http://www.marinespecies.org/rest/AphiaDistributionsByAphiaID/{aphia_id}'
        try:
            req_return = requester(wrms_distribution)
            if req_return is not None:
                wrms_dist = req_return.json()
            else:
                return None,None

            if len(wrms_dist) == 0 :
                print(f'No distribution for Aphia {aphia_id}')
                return None,None
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
            print(f'Error retrieving distribution for aphia {aphia_id}')
            print(err)
            return None,None
