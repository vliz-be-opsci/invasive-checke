#--- Python libs ---
import os
import json
import logging 
import argparse
import traceback 
#--- Pip libs ---
# import numpy as np
import pandas as pd 
import functools
#--- Custom libs ---
from invasive_checker import invasive_checker

'''
This APP takes an input csv file and checks each row of the 
PEMA styled csv (+ additional info) for invasive/non-invasive 
organisms from WRIMS. 

Steps:
    - Read input CSV to pandas + READ semantic file
    - Check all required columns are available
    - Convert sci-name/NCBI column to Aphia-ID (need config on this)
    - Check lon/lat colum + Aphia-ID column for invasiveness
    - write bunch of files to do with invasiveness
        - output semantic file/config/PAPERTRAIL too
'''
log = logging.getLogger('main')  

def process_input_row(row, aphia_checker):
    '''
    Takes pandas/xarray input row and expands it with additional data
    from the invasive_checker lib.

    get aphia_id
    get locations
    get invasiveness
    '''
    # log.info('Processing row {0}'.format(row.Name) )
    sciname = row.get('classification')
    lon = row.get('lon',-64)
    lat = row.get('lat',-45)

    if lat == 0 and lon == 0:
        log.warning('Sample from Null Island!')

    aphia_json = aphia_checker.get_aphia_from_lineage(sciname) 
    aphia_id = aphia_json['AphiaID']
    results, invasive_df = aphia_checker.check_aphia(lon, lat, aphia_id,source='worms')
    if (results is not None):
        item_dict = results
    else: 
        item_dict = {}
        
    if (invasive_df is not None):
        invasive_dict = invasive_df.to_dict(orient='records')
    else: 
        invasive_dict = {}
# derived_status = {  'Status':status,
#                             'Within':mrgid,
#                             'Nearest_Native_MRGID': None,
#                             'Nearest_Introduced_MRGID': None,
#                             'Nearest_Prev_Recorded_MRGID': None,
#                             }
    # to_rvlab
    row['WRIMS Status'] = results.get('Status')
    row['Sample within MRGID'] = results.get('Within')
    row['Details'] = json.dumps(invasive_dict, indent=2)

    return row
 
def make_otu_unique(df):
    '''
    Make OTU's unique for the rvlab
    '''

def do_work(input_file, output_folder, cfg):
    '''
    The meat and potatoes

    Typical Worms input columns:
        - OTU (operational taxonomic unit)
        - n * Sample ID's
        - Sci Classification: NCBI Lineage
        - aphia_id (?Not used?)
        - location_id (Not used)
    Desired extra columns:
        - lat: SRID:4326/WGS84/GPS
        - lon: SRID:4326/WGS84/GPS

    Output Files: 
        - to_rvlab.tsv 
        - classification.csv

    '''
    meta_file = cfg.get('METAFILE')
    log.info('Input file: {0}'.format(input_file))
    log.info('Output folder: {0}'.format(output_folder))
    log.info('Input metadata file: {0}'.format(json.dumps(meta_file, indent=2)))
    log.info('Extra config: {0}'.format(json.dumps(cfg, indent=2)))

    aphia_checker = invasive_checker.Aphia_Checker()

    worms_df = pd.read_csv(input_file, sep = cfg.get('SEP'))
    worms_df = worms_df.apply(lambda row: process_input_row(row, aphia_checker), axis=1)

    # Write CSV containing all data:
    filepath = os.path.join(output_folder, 'classification.csv')
    log.info('Writing full classification file to {0}'.format( filepath))
    worms_df.to_csv(filepath)

    # Write one CSV per status type (Introduced, Native, Holotype etc)
    # status_types = worms_df['WRIMS Status'].unique() 
    # for status in status_types:
    #     filepath = os.path.join(output_folder, 'classification_' + str(status) + '.csv')
    #     log.info('Writing {0} classification file to {1}'.format(status, filepath))
    #     status_df = worms_df[worms_df['WRIMS Status'] == status]
    #     status_df.to_csv(filepath)

    # # Write to_rvlab TSV file 
    # worms_df.groupby(df['OTU'])\
    #   .cumcount()\
    #   .astype(str)\
    #   .str.replace('0', '')\
    #   .values

    # # Write metadata file
    # meta_df  = pd.read_csv(meta_file)
 

def get_config():
    '''
    Get the script config from the environment files
    '''
    cfg = { 'LLEVEL':os.getenv('LLEVEL', 'INFO'),
            'ID_SOURCE':os.getenv('ID_SOURCE', 'sciname'),
            'METAFILE':os.getenv('METAFILE', 'INFO'),
            'SEP':os.getenv('SEP', '\t'),}

    return cfg

def main(args):
    '''
    Setup logging, and args, then "do_work" via a schedular
    '''
    loglevel = os.getenv('LLEVEL', default='INFO')
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
                        level=getattr(logging, loglevel))    
    try: 
        cfg = get_config() 
        do_work(args.input_file, args.output_folder, cfg)
    except (KeyboardInterrupt, SystemExit):
        log.warning('Exiting script...')
        pass

    log.info('Script Ended...') 

if __name__ == "__main__":
    '''
    This takes the command line args and passes them to the 'main' function
    '''
    PARSER = argparse.ArgumentParser(
        description='Invasive Checker')
    PARSER.add_argument(
        '-l', '--loglevel', default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help="Set log level for service (%s)" % 'INFO')
    PARSER.add_argument(
        '-i', '--input_file', 
        help="Path to input csv file.")
    PARSER.add_argument(
        '-o', '--output_folder', default='/tmp/',
        help="Path to folder to write output files.")
    ARGS = PARSER.parse_args()
    try:
        main(ARGS)
    except KeyboardInterrupt:
        log.warning('Keyboard Interrupt. Exiting...')
        os._exit(0)
    except Exception as error:
        log.error('Other exception. Exiting with code 1...')
        log.error(traceback.format_exc())
        log.error(error)
        os._exit(1)
