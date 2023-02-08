from pathlib import Path
from typing import Optional

import awkward as ak
import uproot
from tqdm import tqdm
from glob import glob

from hgcal_quant import conf

# event branches
branches = ["runNumber", "eventNumber", "luminosityNumber"]

# # gen particle
# branches += [
#     "genPart_n",
#     "genPart_E",
#     "genPart_px",
#     "genPart_py",
#     "genPart_pz",
#     "genPart_pT",
#     "genPart_eta",
#     "genPart_phi",
# ]
# # simhits
# branches += [
#     "simHit_n",
#     "simHit_E",
#     "simHit_x",
#     "simHit_y",
#     "simHit_z",
#     "simHit_eta",
#     "simHit_phi",
#     "simHit_layer",
#     "simHit_zside",
#     "simHit_detector",
#     "simHit_detId",
# ]

# rechits
branches += [
    "recHit_n",
    "recHit_E",
    "recHit_x",
    "recHit_y",
    "recHit_z",
    "recHit_eta",
    "recHit_phi",
    "recHit_layer",
    "recHit_zside",
    "recHit_detector",
    "recHit_detId",
]


rootprefix = "treeMaker/tree"


def readpath(
    fn: Path,
) -> ak.highlevel.Array:
    with uproot.open(fn) as rfile:
        roottree = rfile[rootprefix]
        array = roottree.arrays(branches, library="ak")
        # Todo
        # Some cells don't have correct eta value
        # Alternative: a) don't use eta b) recompute eta
        nan_mask = ~ak.is_none(ak.nan_to_none(ak.sum(array["recHit_eta"], 1)))
        z_mask = ak.all(array["recHit_z"] > 200, 1)
        return array[nan_mask & z_mask]


filelist = glob(str(Path(conf.path.dataset).expanduser()))
ds = ak.concatenate([readpath(Path(fn)) for fn in tqdm(filelist)])
event0 = ds[0]
