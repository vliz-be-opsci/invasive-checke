from fastapi import FastAPI, Query
import os
import json
import logging 
import numpy as np
from invasive_checker import invasive_checker

class CustomJSONizer(json.JSONEncoder):
    def default(self, obj):
        return super().encode(bool(obj)) \
            if isinstance(obj, np.bool_) \
            else super().default(obj)

log = logging.getLogger('main') 
loglevel = os.getenv('LLEVEL', default='INFO')
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
                    level=getattr(logging, loglevel))
log.setLevel(getattr(logging,loglevel))

app = FastAPI()
aphia_checker = invasive_checker.Aphia_Checker() 

@app.get("/")
async def root():
    return {"message": "Aphia Checker. Provide Sample Lon/Lat and AphiaID of organism detected there"}

@app.get("/check/{aphia_id}")
async def read_user_item(aphia_id: str, 
                        source: Query('worms', enum=['worms', 
                                                    'algaebase',
                                                    'bold',
                                                    'dyntaxa',
                                                    'fishbase',
                                                    'iucn',
                                                    'lsid',
                                                    'ncbi',
                                                    'tsn',
                                                    'gisd']),
                         lon: float, 
                         lat: float):
    aa, bb = aphia_checker.check_aphia(lon, lat, aphia_id,source=source)
    item_dict = {'summary':aa,
            'details':bb.fillna('').to_dict(orient='records')}
    return item_dict 

@app.get("/check/from_sciname/{sciname}")
async def from_sciname(sciname: str, lon: float, lat: float):
    aphia_json = aphia_checker.get_aphia_from_taxname(sciname) 
    aphia_id = aphia_json['AphiaID']
    aa, bb = aphia_checker.check_aphia(lon, lat, aphia_id,source='worms')
    item_dict = {'summary':aa,
            'details':bb.fillna('').to_dict(orient='records')}
    return item_dict

@app.post("/clear_cache")
async def clear_geom():
    clear_result = aphia_checker.clear_cache()
    return clear_result