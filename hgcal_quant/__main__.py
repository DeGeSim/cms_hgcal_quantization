# %%
# from hgcal_quant.events import event0
from hgcal_quant.geo import geo_lut
import numpy as np

# %%
import matplotlib.pyplot as plt
import matplotlib as mpl

# default value for missing cells
geo_lut.loc[0] = np.ones(len(geo_lut.columns)) * -9999

refcell = geo_lut[(geo_lut.detectorid == 8) & (geo_lut.celltype == 0)].iloc[0]
neighbors = [geo_lut.loc[refcell[f"n{i}"]] for i in range(6)]

delta_x = np.abs(min([refcell.x - n.x for n in neighbors]))
delta_y = np.abs(min([refcell.y - n.y for n in neighbors]))

# %%


def postr(cell):
    pos = np.array([cell[v] for v in ["x", "y"]])
    # pos[0] = pos[0] / delta_x
    # pos[1] = pos[1] / delta_y  # + pos[0]
    return pos


pos = postr(refcell)

# %%
def smearcell(cells):
    # get the xy pos of the neighbors
    xy = geo_lut.loc[cells][["x", "y"]].to_numpy(dtype="float32")
    # detids of the neihgbors
    xynbghidx = geo_lut.loc[cells][[f"n{i}" for i in range(12)]]
    # save which ones are real to mask them out later
    xynb_mask = (xynbghidx != 0).values.reshape(-1, 12, 1).repeat(2, -1)

    # get the xy values of the neigbors
    # geo_lut.x.get(xynbghidx.values.reshape(-1),0)
    xyneigbors_delta = geo_lut[["x", "y"]].loc[
        xynbghidx.values.reshape(-1)
    ].values.reshape(-1, 12, 2) - xy.reshape(-1, 1, 2).repeat(12, 1)

    scale = (
        np.random.beta(2, 2, (len(xyneigbors_delta), 12, 1)).repeat(2, -1)
    ) ** (0.5)
    scale = scale / (scale * xynb_mask).sum() * len(cells)
    # scale = np.zeros((len(xyneigbors_delta), 1))
    # scale[0, 0] = 1.0
    xy += (xyneigbors_delta * scale * xynb_mask).sum(1)
    layer = (geo_lut.loc[cells]["layer"].values + np.random.rand(len(xy))).reshape(
        -1, 1
    )
    return np.hstack([xy, layer])


scatter = smearcell([refcell.name] * 10000)[..., :2]

# %%
fig, ax = plt.subplots(1, 1)

h = ax.hist2d(
    scatter[..., 0],
    scatter[..., 1],
    bins=[np.linspace(-42.2, -40.5, 100), np.linspace(-5.8, -3.9, 100)],
)[3]
fig.colorbar(h, ax=ax)

# plt.scatter([0], [0])
for ineighbor, neighbor in enumerate([refcell] + neighbors):
    npos = postr(neighbor)
    # dist = np.round(pos - npos, 2)

    c = [mpl.cm.Set1.colors[ineighbor]]
    c = "white" if ineighbor == 0 else "red"

    ax.scatter(npos[0], npos[1], c=c, s=100)
    ax.annotate(str(npos), npos, c=c)
    # scatter = dist.reshape(1, 2).repeat(100, axis=0) + np.random.uniform(
    #     size=(100, 2)
    # )
    # ax.scatter(scatter[:, 0], scatter[:, 1], c=c)
    # print(ineighbor, dist, c)


celltype_to_dist = {0: 0.8, 1: 1.2, 2: 1.2}
# %%
plt.hist(scatter[:, 0], bins=100)
xpos = [e.x for e in neighbors]
xmin = min(xpos)
xmax = max(xpos)
plt.axvline(xmin + (xmax - xmin) * 1 / 4, c="red")
plt.axvline(xmin + (xmax - xmin) * 3 / 4, c="red")
# %%
