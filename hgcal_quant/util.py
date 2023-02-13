from hgcal_quant.geo import geo_lut
import numpy as np


def get_neighbors(cell):
    nl = [geo_lut.loc[cell[f"n{i}"]] for i in range(12)]
    return list(filter(lambda x: x.name != 0, nl))


def rms(a):
    return np.sqrt(np.power(a, 2).sum(2)).mean(1)


def postr(cell):
    pos = np.array([cell[v] for v in ["x", "y"]])
    # pos[0] = pos[0] / delta_x
    # pos[1] = pos[1] / delta_y  # + pos[0]
    return pos
