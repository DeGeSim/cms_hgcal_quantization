from pathlib import Path

import awkward as ak
import pandas as pd
import uproot

from hgcal_quant import conf

pickle_lut_path = Path(conf.path.geo_lut).with_suffix(".pd")
if not pickle_lut_path.is_file():
    with uproot.open(conf.path.geo_lut) as rf:
        geo_lut = rf["analyzer/tree;1"].arrays(library="ak")
    geo_lut = ak.to_dataframe(geo_lut)
    geo_lut.set_index("globalid", inplace=True)
    # Fix the layers being the same for all subdetector
    # assign negative values for z<0
    factor = 1 - 2 * (geo_lut.z < 0)
    # get the layerid within the subdetector
    layerid = geo_lut.layerid
    # calculate the needed shift
    detectors = sorted(geo_lut.detectorid.value_counts().index.values)
    # check that there is only the 8 detectors
    assert detectors == [8, 9, 10]
    # shift all layers from other subdetectors by the layers in the EE
    shift = layerid[geo_lut.detectorid == 8].max() * (geo_lut.detectorid != 8)
    geo_lut["layer"] = factor * (layerid + shift)
    geo_lut.to_pickle(pickle_lut_path)


geo_lut = pd.read_pickle(pickle_lut_path)
