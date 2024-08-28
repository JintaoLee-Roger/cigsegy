# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.
"""
Plotting functions for visualizing survey coordinates, trace header keys, and other SEG-Y related data.

This module provides a collection of plotting functions tailored for geophysical data visualization. 
These functions are designed to help users graphically represent various aspects of SEG-Y files, 
such as survey coordinates and trace header keys. 

**The module may require additional dependencies for full functionality.**
"""

import numpy as np
from typing import Tuple
from cigse.cpp._CXX_SEGY import Pysegy
from cigse import ExceptionWrapper, tools
from cigse import SegyNP

try:
    import matplotlib.pyplot as plt
except BaseException as E:
    plt = ExceptionWrapper(E, "run `pip install matplotlib` to install the dependency") # yapf: disable


try:
    import cigvis
except BaseException as E:
    cigvis = ExceptionWrapper(E, "run `pip install cigvis` to install the dependency") # yapf: disable


def plot_region(fname: str,
                mode: str = 'line',
                save: str = None,
                **kwargs) -> None:
    """
    plot the region map (x and y axis are inline and crossline)

    Parameters
    ----------
    segy : str or Pysegy
        input segy file
    loc : list
        contains 4 values, as [iline, xline, istep, xstep]
    save : str
        save to a png image
    """
    if mode not in ['line', 'cdpxy', 'xy']:
        raise RuntimeError(f"mode must be 'line' or 'cdpxy'/'xy', mode = {mode}") # yapf: disable

    geom = tools.get_lineInfo(fname, mode='geom', **kwargs)
    if mode == 'line':
        x, y = geom[:, 0], geom[:, 1]
    else:
        x, y = geom[:, 2], geom[:, 3]

    plt.fill(x, y, color=(0.9, 0.9, 0.9))
    plt.plot(x, y)
    plt.gca().invert_yaxis()
    # plt.gca().xaxis.set_ticks_position('top')

    plt.grid(True, linestyle='--')
    if mode == 'line':
        xlabel = f"Inline Number"
        ylabel = f"Crossline Number"
    else:
        xlabel = f"X"
        ylabel = f"Y"

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Region')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches='tight', pad_inches=0.0)
    plt.show()


def plot_trace_keys(fname: str,
                    keyloc: int,
                    beg: int = 0,
                    end: int = 1000,
                    save: str = None) -> None:
    """
    plot the values (at keyloc in each trace) of the traces 
    range from beg to end .
    """
    if isinstance(fname, Pysegy):
        segy = fname
    else:
        segy = Pysegy(str(fname))
    assert beg >= 0 and end > beg + 1 and end <= segy.ntrace, "invalid beg and end"
    keys = segy.get_trace_keys([keyloc], [4], beg, end).squeeze()
    x = np.arange(beg, end)
    plt.plot(x, keys)
    plt.xlabel("Trace Number")
    plt.ylabel(f"Key values (at {keyloc})")
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches='tight', pad_inches=0.0)
    plt.show()


def plot_trace_ix(fname: str,
                  iline: int,
                  xline: int,
                  beg: int = 0,
                  end: int = 1000,
                  figsize: Tuple = None,
                  save: str = None) -> None:
    """
    plot inline and crossline number of the traces range from
    beg to end 
    """
    if isinstance(fname, Pysegy):
        segy = fname
    else:
        segy = Pysegy(str(fname))
    assert beg >= 0 and end > beg + 1, "invalid beg and end"
    assert end < segy.ntrace, "end > trace_count"
    keys = segy.get_trace_keys([iline, xline], [4, 4], beg, end).squeeze()
    iv = keys[:, 0]
    xv = keys[:, 1]
    x = np.arange(beg, end)

    if figsize is None:
        figsize = (9, 4)
    plt.figure(figsize=figsize)
    plt.subplot(1, 2, 1)
    plt.plot(x, iv)
    plt.xlabel('Trace Number')
    plt.ylabel(f"Inline (at {iline})")
    plt.grid(True, linestyle='--')
    plt.subplot(1, 2, 2)
    plt.plot(x, xv)
    plt.xlabel('Trace Number')
    plt.ylabel(f"Crossline (at {xline})")
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches='tight', pad_inches=0.0)
    plt.show()


def plot_trace_ixo(fname: str,
                   iline: int,
                   xline: int,
                   offset: int = 37,
                   beg: int = 0,
                   end: int = 1000,
                   figsize: Tuple = None,
                   save: str = None) -> None:
    """
    plot inline, crossline and offset number of the traces range from
    beg to end 
    """
    if isinstance(fname, Pysegy):
        segy = fname
    else:
        segy = Pysegy(str(fname))
    assert beg >= 0 and end > beg + 1, "invalid beg and end"
    assert end < segy.ntrace, "end > trace_count"
    keys = segy.get_trace_keys([iline, xline, offset], [4, 4, 4], beg, end)
    iv = keys[:, 0]
    xv = keys[:, 1]
    ov = keys[:, 2]
    x = np.arange(beg, end)

    if figsize is None:
        figsize = (12, 4)
    plt.figure(figsize=figsize)
    plt.subplot(1, 3, 1)
    plt.plot(x, iv)
    plt.xlabel('Trace Number')
    plt.ylabel(f"Inline (at {iline})")
    plt.grid(True, linestyle='--')
    plt.subplot(1, 3, 2)
    plt.plot(x, xv)
    plt.xlabel('Trace Number')
    plt.ylabel(f"Crossline (at {xline})")
    plt.grid(True, linestyle='--')
    plt.subplot(1, 3, 3)
    plt.plot(x, ov)
    plt.xlabel('Trace Number')
    plt.ylabel(f"offset (at {offset})")
    plt.grid(True, linestyle='--')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches='tight', pad_inches=0.0)
    plt.show()


def plot3d(fname: str):
    d = SegyNP(fname)
    nodes = cigvis.create_slices(d)
    cigvis.plot3D(nodes)
