import logging

import pandas as pd

log = logging.getLogger("utils")


def get_sample_location_df(AccessionIDs, meta_df):
    """
    Return a dataframe of AccessionID metadata.
    """
    row = {}
    for sample_id in AccessionIDs:
        log.info("Getting metadata for {0}".format(sample_id))
        meta_row = meta_df[meta_df.isin([sample_id]).any(axis=1)]
        meta_cols = meta_row[meta_row.isin([sample_id])].dropna(axis="columns")
        gene = meta_cols.columns.values[0]
        if gene.startswith("negativeControl_gene"):
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
                log.debug("   -Possible duplicate gene ID...")
                longitude = None
                latitude = None
                material_id = None
            elif meta_row.shape[0] == 0:
                log.warning("   -Empty row?...")
                longitude = None
                latitude = None
                material_id = None
            else:
                longitude = meta_row.loc[meta_row.index].longitude.values[0]
                latitude = meta_row.loc[meta_row.index].latitude.values[0]
                material_id = meta_row.loc[meta_row.index].Sample_ID.values[0]

        row[sample_id] = {
            "AccessionNumber": sample_id,
            "isNegativeControlGene": negative_control,
            "GeneType": gene,
            "latitude": latitude,
            "longitude": longitude,
            "Sample_ID": material_id,
        }
    my_df = pd.DataFrame(row).T
    return my_df
