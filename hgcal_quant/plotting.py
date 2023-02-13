import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from util import postr
import os

plt.rcParams["savefig.bbox"] = "tight"
# plt.rcParams["backend"] = "Agg"
plt.rcParams["figure.dpi"] = 150
np.set_printoptions(formatter={"float_kind": "{:.3g}".format})


def all_plots(*, dir, **kwargs):
    os.makedirs(dir, exist_ok=True)
    for f, name in [
        (scatter, "scatter"),
        (hist2d, "2dhist"),
        (yhist, "yhist"),
        (xhist, "xhist"),
    ]:
        fig: mpl.figure.Figure = f(**kwargs)
        fig.tight_layout()
        fig.savefig(f"{dir}/{name}.png")
        fig.clf()


def scatter(*, scatter, names, **kwargs):
    fig, ax = plt.subplots(1, 1)
    for name, arr in list(zip(names, scatter))[::-1]:
        ax.scatter(arr[..., 0], arr[..., 1], alpha=0.1, label=name)
    ax.legend()
    return fig


def hist2d(*, scatter, names, nnames, refcell, neighbors, nneighbors, **kwargs):
    fig, ax = plt.subplots(1, 1)
    h = ax.hist2d(
        scatter[..., 0].reshape(-1),
        scatter[..., 1].reshape(-1),
        bins=100,  # [np.linspace(-42.2, -40.5, 100), np.linspace(-5.8, -3.9, 100)],
        # norm=mpl.colors.LogNorm(),
    )[3]
    fig.colorbar(h, ax=ax)
    for ineighbor, neighbor in enumerate([refcell] + neighbors + nneighbors):
        npos = postr(neighbor)
        # c = [mpl.cm.Set1.colors[ineighbor]]
        # c = "brown" if ineighbor == 0 else "red"
        ax.scatter(
            npos[0],
            npos[1],
            s=100,
            label=(names + nnames)[ineighbor],
            edgecolors="black",
        )
    ax.legend(bbox_to_anchor=(1.5, 1.0))
    return fig


def xhist(*, scatter, names, neighbors, **kwargs):
    fig, ax = plt.subplots(1, 1)
    bins = 100  # np.linspace(-43.2, -39.8, 100)

    arrs = [e[..., 0] for e in scatter] + [scatter[..., 0].reshape(-1)[::3]]
    # for name, arr in zip(
    #     names,arrs
    # ):
    ax.hist(
        arrs,
        bins=bins,
        histtype="step",
        # stacked=True,
        label=names + ["sum/3"],
    )
    xmin, xmax = [f([e.x for e in neighbors]) for f in [min, max]]

    ax.axvline(xmin, c="green")
    ax.axvline(xmax, c="green")
    ax.axvline(xmin + (xmax - xmin) * 1 / 4, c="red")
    ax.axvline(xmin + (xmax - xmin) * 3 / 4, c="red")
    ax.set_title("x Histogramm")
    ax.legend(prop={"size": 10})
    fig.tight_layout()
    return fig


def yhist(*, scatter, names, **kwargs):
    fig, ax = plt.subplots(1, 1)
    arrs = [e[..., 1] for e in scatter] + [scatter[..., 1].reshape(-1)[::3]]
    ax.hist(
        arrs,
        bins=100,
        histtype="step",
        label=names + ["sum/3"],
    )
    # ymin, ymax = [f([e.y for e in neighbors]) for f in [min, max]]

    # ax.axvline(ymin, c="green")
    # ax.axvline(ymax, c="green")
    # ax.axvline(ymin + (ymax - ymin) * 1 / 4, c="red")
    # ax.axvline(ymin + (ymax - ymin) * 3 / 4, c="red")
    ax.set_title("y Histogramm")
    ax.legend(prop={"size": 10})
    fig.tight_layout()
    return fig
