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
from cigsegy.cpp._CXX_SEGY import Pysegy
from cigsegy import ExceptionWrapper, tools
from cigsegy.interp import arbitray_line

try:
    import matplotlib.pyplot as plt
except BaseException as E:
    plt = ExceptionWrapper(E, "run `pip install matplotlib` to install the dependency") # yapf: disable


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

    k1, k2 = [x[0], y[0]], [x[-2], y[-2]]
    s = (x.max() - x.min()) / 50
    plt.arrow(k1[0], k1[1], k2[0]-k1[0], k2[1]-k1[1], head_width=s, head_length=s, ec='red', fc='red', zorder=4)

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


def extract_arbitrary_line_by_view(data: np.ndarray,
                                   bmap: str = 'data',
                                   draw_arb: bool = True,
                                   *,
                                   line: bool = True,
                                   idx: int = 50,
                                   cline='#F3AA3C'):
    """
    extract arbitrary line from seismic data by clicking

    Parameters
    ----------
    - data: np.ndarray 
        3D seismic data
    - bmap: str
        background map, 'data' or 'blank'
    - line : bool
        whether to draw the broken line 
    - idx: int 
        the slice index of the seismic data if bmap is 'data'
    - cline: str
        color of the line

    Returns
    -------
    - out: np.ndarray
        extracted arbitrary line
    - p: np.ndarray
        extracted arbitrary line path
    - coords: np.ndarray
        the coordinates by clicking
    """
    fig, ax = plt.subplots()
    if bmap == 'data':
        img = data[:, :, idx]
        ax.imshow(img.T, cmap='gray', aspect='auto')
    else:
        img = np.zeros((data.shape[1], data.shape[0], 4), dtype=np.uint8)
        img = img + int(0.9 * 255)
        ax.imshow(img, aspect='auto')
        ax.grid(True)

    ax.set_title('Obtain a path by clicking, press "u" to undo, "enter" to finish') # yapf: disable
    ax.set_xlabel('Inline')
    ax.set_ylabel('Crossline')
    coords = []
    lines = []
    points = []
    out = None
    p = None
    indices = None

    def _draw_arbline():
        if (len(coords) > 0):
            ax.clear()
            nonlocal out, p, indices
            out, p, indices = arbitray_line(data, coords)
            ax.grid(False)
            ax.imshow(out.T, cmap='gray', aspect='auto')
            if line:
                for i in range(len(indices)):
                    ax.plot([indices[i]] * out.shape[1],
                            range(out.shape[1]),
                            'w--',
                            lw=0.5)
                    if i > 0 and i < len(indices) - 1:
                        ax.text(indices[i] - out.shape[0] / 50,
                                30,
                                str(indices[i]),
                                color='w',
                                fontsize=10)
            ax.set_title("Arbitrary Line")
            ax.set_ylabel('Time')
            fig.canvas.draw()
            # plt.close()

    def _click_event(event):
        if event.inaxes:
            coords.append((round(event.xdata, 2), round(event.ydata, 2)))
            p = ax.plot(event.xdata, event.ydata, 'ro')[0]
            points.append(p)
            if len(coords) > 1:
                x_coords, y_coords = zip(*coords[-2:])
                line, = ax.plot(x_coords, y_coords, c=cline)
                lines.append(line)
            fig.canvas.draw()

    def _undo_last(event):
        if event.key == 'u':
            if coords:
                coords.pop()
                p = points.pop()
                p.remove()
                if len(lines) > 0:
                    l = lines.pop()
                    l.remove()
                fig.canvas.draw()
        if event.key in ('enter', 'return', 'escape'):
            if draw_arb:
                fig.canvas.mpl_disconnect(cid_click)
                fig.canvas.mpl_disconnect(cid_key)
                _draw_arbline()
            else:
                plt.close(fig)

    cid_click = fig.canvas.mpl_connect('button_press_event', _click_event)
    cid_key = fig.canvas.mpl_connect('key_press_event', _undo_last)

    plt.show()

    return out, p, indices, np.array(coords)
