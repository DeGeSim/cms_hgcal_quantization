# %%
# from hgcal_quant.events import event0
from hgcal_quant.geo import geo_lut
import numpy as np

# %%
import matplotlib.pyplot as plt
import matplotlib as mpl

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
def smearcell(cell):
    xy = cell[["x", "y"]].to_numpy(dtype="float32")
    xynbghidx = [e for e in [cell[f"n{i}"] for i in range(12)] if e != 0]
    xyneigbors_delta = geo_lut[["x", "y"]].loc[xynbghidx].values - xy
    scale = (np.random.beta(2, 2, (len(xyneigbors_delta), 1))) ** 4
    scale /= scale.sum()
    # scale = np.zeros((len(xyneigbors_delta), 1))
    # scale[0, 0] = 1.0
    xy += (xyneigbors_delta * scale.repeat(2, 1)).sum(0)
    layer = cell["layer"] + np.random.rand()
    return np.hstack([xy, layer])


scatter = np.stack([smearcell(refcell)[:2] for _ in range(10_000)])

# %%
fig, ax = plt.subplots(1, 1)

h = ax.hist2d(
    scatter[:, 0],
    scatter[:, 1],
    bins=[np.linspace(-42.2, -40.5, 100), np.linspace(-5.8, -3.9, 100)],
)
fig.colorbar(h[3], ax=ax)

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
plt.hist(scatter[:, 0], bins=10)
# %%
