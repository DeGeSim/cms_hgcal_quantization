from hgcal_quant.geo import geo_lut
import numpy as np

for celltype in geo_lut.celltype.unique():
    refcell = geo_lut[
        (geo_lut.detectorid == 8) & (geo_lut.celltype == celltype)
    ].iloc[0]
    pos = np.array([refcell[v] for v in ["x", "y", "z"]])
    for ineighbor in range(6):
        neighbor = geo_lut.loc[refcell[f"n{ineighbor}"]]
        npos = np.array([neighbor[v] for v in ["x", "y", "z"]])
        dist = np.sqrt(((pos - npos) ** 2).sum())
        print(celltype, dist)

celltype_to_dist = {0: 0.8, 1: 1.2, 2: 1.2}
