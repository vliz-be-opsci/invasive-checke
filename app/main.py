#--- Python libs ---
import os
import json
import logging 
import argparse
import traceback 
from hashlib import md5
#--- Pip libs ---
# import numpy as np
import pandas as pd 
#--- Custom libs ---
from invasive_checker import invasive_checker, utils 

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

def process_input_row(row):
    '''
    Takes pandas/xarray input row and expands it with additional data
    from the invasive_checker lib.

    get aphia_id
    get locations
    get invasiveness
    '''
    log.info('Processing row {0}'.format(row.name) )
    log.debug('row: {0}'.format(row))
    sciname = row.get('classification')
    lon = row.get('sampleLongitude') or row.get('Longitude') or row.get('longitude') or row.get('lon') or 0
    lat = row.get('sampleLatitude') or row.get('Latitude') or row.get('latitude') or row.get('lat') or 0

    if lat == 0 and lon == 0:
        log.warning('Sample from Null Island! Lat=Lon=0')

    aphia_json = invasive_checker.get_aphia_from_lineage(sciname) 
    if aphia_json is not None:
        aphia_id = aphia_json.get('AphiaID')
        sci_name = aphia_json.get('scientificname')
        sci_rank = aphia_json.get('rank')
    else:
        aphia_id = 'No Match'
        sci_name = 'No Match'
        sci_rank = 'No Match'

    if row['isNegativeControlGene']:
        log.warning('Sample is Negative Control: not bothering with invasiveness check...')
    elif aphia_id is not None:
        results, invasive_df = invasive_checker.check_aphia(lon, lat, aphia_id,source='worms')
        if (results is not None): 
            wrims_status = results.get('Status',['Unrecorded'])
            within = results.get('Within',['None'])
        else: 
            wrims_status = ['Unrecorded']
            within = ['None']
            
        if (invasive_df is not None):
            invasive_dict = invasive_df.to_dict(orient='records')
        else: 
            invasive_dict = {}

        # to_rvlab
        row['Aphia_ID'] = aphia_id
        row['Worms SciName'] = sci_name
        row['Worms SciName Rank'] = sci_rank
        row['WRIMS Status at Sample Location'] = wrims_status
        row['MarineRegions with known occurrence at Sample Location'] = within

        # row['Details'] = json.dumps(invasive_dict, indent=2)
    return row

def do_work(input_file, output_folder, meta_file, cfg):
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
    log.info('  -Preparing data...')

    with open(f"{output_folder}/invasive_checker.log", "w", encoding="utf8") as f:
        with open(input_file, "rb") as fi:
            f.write(f"{md5(fi.read()).hexdigest()}\t(md5sum of {input_file})\n")
        with open(meta_file, "rb") as fm:
            f.write(f"{md5(fm.read()).hexdigest()}\t(md5sum of {meta_file})\n")


    worms_df = pd.read_csv(input_file, sep = cfg.get('SEP'))
    worms_df = clean_up_dataframes(worms_df, cfg)

    worms_df_unpivot = pd.melt(worms_df, id_vars=['classification','OTU'],var_name='AccessionID', value_name='Count')
    meta_df = pd.read_csv(meta_file) 
    sample_df = utils.get_sample_location_df(worms_df_unpivot.AccessionID.unique(),meta_df)
    worms_df_unpivot = pd.merge(worms_df_unpivot,sample_df,how="left",left_on='AccessionID',right_on='AccessionNumber')
    
    os.makedirs('/mnt/tests/output/', exist_ok=True)

    worms_df_unpivot.to_csv('/mnt/tests/output/unpivot.csv',index=False)
    log.info('  -Looping through rows...') 
    wrims_df = worms_df_unpivot.apply(lambda row: process_input_row(row), axis=1)
    wrims_df.to_csv('/mnt/tests/output/wrims_df.csv',index=False)

    # Clean up table
    wrims_df = wrims_df.drop('AccessionNumber',axis=1)

    # Write input + aphia file
    aphia_df = pd.merge(worms_df,wrims_df[['classification','Aphia_ID']],how="left",on='classification')
    aphia_df = aphia_df.drop_duplicates()
    aphia_filepath = os.path.join(output_folder, cfg.get('WORMS_OUTPUT_FILE'))
    
    if aphia_df['OTU'].is_unique:
        log.info('Writing worms.csv classification file to {0}'.format( aphia_filepath))
        aphia_df.to_csv(aphia_filepath,sep = cfg.get('SEP'),index=False)
    else:
        log.warning('Non-unique OTU column. Will add suffix to duplicate OTUs...') 

        aphia_df['duplicate_otu'] = aphia_df['OTU'].duplicated()
        aphia_df['duplicate_count'] = aphia_df.groupby(aphia_df.OTU).cumcount().values
        aphia_df['OTU'] =  aphia_df.apply(lambda row: row['OTU'] + '_' + str(row['duplicate_count']) if row['duplicate_otu'] else row['OTU'],axis=1)
        aphia_df.drop(['duplicate_otu','duplicate_count'],axis=1,inplace=True)
        log.info('Writing worms.csv classification file to {0}'.format( aphia_filepath))
        aphia_df.to_csv(aphia_filepath,sep = cfg.get('SEP'),index=False)

    # Write CSV containing classification data:
    filepath = os.path.join(output_folder, cfg.get('CLASS_OUTPUT_FILE'))
    log.info('Writing full classification file to {0}'.format( filepath))
    wrims_df = wrims_df[wrims_df['Count'] > 0]
    wrims_df.to_csv(filepath,index=False)

# take a dataframe and change several column names to match column names defined in the cfg file
def clean_up_dataframes(df, cfg):

    otu_col_name = cfg.get('OTU_COL_NAME')
    class_col_name = cfg.get('CLASS_COL_NAME')
    log.info(df.columns)
 
    # df['classification'] = df[class_col_name]
    # df['OTU'] = df[otu_col_name]
    
    df.rename(columns={class_col_name:'classification',
                           otu_col_name: 'OTU'}, inplace=True)
    df = df.drop(['location_id','aphia_id'], axis=1, errors='ignore')
    log.info(df.columns)
    log.info('Done preparing data...')
    
    return df

def get_config():
    '''
    Get the script config from the environment files
    '''
    cfg = { 'LLEVEL':os.getenv('LLEVEL', 'INFO'),
            'ID_SOURCE':os.getenv('ID_SOURCE', 'sciname'),
            # 'METAFILE':os.getenv('METAFILE', 'INFO'),
            'SEP':os.getenv('SEP', '\t'),
            'OTU_COL_NAME':os.getenv('OTU_COL_NAME', 'OTU'),
            'CLASS_COL_NAME':os.getenv('CLASS_COL_NAME', 'classification'),
            'CLASS_OUTPUT_FILE':os.getenv('CLASS_OUTPUT_FILE', 'classification.csv'),
            'WORMS_OUTPUT_FILE':os.getenv('WORMS_OUTPUT_FILE', 'worms.csv'),}
    return cfg

def main(args):
    '''
    Setup logging, and args, then "do_work" via a scheduler
    '''
    loglevel = args.loglevel
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
                        level=getattr(logging, loglevel))    
    try: 
        cfg = get_config() 
        log.info('Input file: {0}'.format(args.input_file))
        log.info('Output folder: {0}'.format(args.output_folder))
        log.info('Input metadata File: {0}'.format(args.meta_file))
        log.info('Extra config: {0}'.format(json.dumps(cfg, indent=2)))
        do_work(args.input_file, args.output_folder, args.meta_file, cfg)
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
        '-m', '--meta_file', 
        help="Path to sample metadata csv file.")
    PARSER.add_argument(
        '-o', '--output_folder', default='/mnt/',
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
