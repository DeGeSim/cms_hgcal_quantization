# %%
%load_ext autoreload
%autoreload 2
from hgcal_quant.geo import geo_lut
import numpy as np

from util import get_neighbors, postr, rms


celltype_to_dist = {0: 0.8, 1: 1.2, 2: 1.2}
# default value for missing cells
geo_lut.loc[0] = np.ones(len(geo_lut.columns)) * -9999

refcell = geo_lut[(geo_lut.detectorid == 8) & (geo_lut.celltype == 0)].iloc[0]


neighbors = get_neighbors(refcell)
nneighbors = [
    geo_lut.loc[detid]
    for detid in {e.name for ne in neighbors for e in get_neighbors(ne)}
    if detid != refcell.name and detid not in {e.name for e in neighbors}
]

names = ["center"] + [f"n{i+1}" for i in range(len(neighbors))]
nnames = [f"n{i+1}" for i in range(len(nneighbors))]

delta_x = np.abs(min([refcell.x - n.x for n in neighbors]))
delta_y = np.abs(min([refcell.y - n.y for n in neighbors]))

# %%


pos = postr(refcell)

# %%


def delta_benno(nvec, nmask):
    r = np.sqrt(rms(nvec * nmask).reshape(-1, 1) * np.random.rand(len(nvec), 1) / 2)
    phi = np.random.rand(*r.shape) * 2 * np.pi
    dx = r * np.sin(phi)
    dy = r * np.cos(phi)
    delta = np.hstack([dx, dy])
    return delta


def delta_uniform_avg(nvec, nmask):
    # offcenter correction
    nvec = nvec - (nvec * nmask).sum(1).reshape(-1, 1, 2).repeat(12, 1) / 12
    # scale = np.random.beta(0.5, 0.5, (len(nvec), 12, 1))
    scale = np.random.rand(len(nvec), 12, 1) ** 2
    scale *= 1 / (scale * nmask[..., :1]).sum(1).reshape(-1, 1, 1).repeat(12, 1)
    scale = (scale / (scale * nmask[..., :1]).sum()) * nvec.shape[0]
    scale = scale ** (1 / 2)
    scale = scale.repeat(2, -1)
    delta = (nvec * scale * nmask).sum(1)
    return delta


def smearcell(cells):
    # get the xy pos of the neighbors
    xy = geo_lut.loc[cells][["x", "y"]].to_numpy(dtype="float32")
    # detids of the neihgbors
    xynbghidx = geo_lut.loc[cells][[f"n{i}" for i in range(12)]]
    # save which ones are real to mask them out later
    nmask = (xynbghidx != 0).values.reshape(-1, 12, 1).repeat(2, -1)

    # get the xy values of the neigbors
    nvec = geo_lut[["x", "y"]].loc[xynbghidx.values.reshape(-1)].values.reshape(
        -1, 12, 2
    ) - xy.reshape(-1, 1, 2).repeat(12, 1)

    # delta = delta_uniform_avg(nvec, nmask)
    delta = delta_benno(nvec, nmask)

    print(delta.mean(0))

    xy += delta
    layer = (geo_lut.loc[cells]["layer"].values + np.random.rand(len(xy))).reshape(
        -1, 1
    )
    return np.hstack([xy, layer])


cells = [refcell] + neighbors + nneighbors
ncells = len(cells)
scatter = (
    smearcell([e.name for e in cells] * 10000)[..., :2]
    .reshape(-1, ncells, 2)
    .swapaxes(0, 1)
)


# %%


from plotting import all_plots

all_plots(
    dir="plots/disks",
    scatter=scatter,
    names=names,
    nnames=nnames,
    refcell=refcell,
    neighbors=neighbors,
    nneighbors=nneighbors,
)
# %%
