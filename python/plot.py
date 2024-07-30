import numpy as np
from typing import List, Tuple, Dict, Union
from pathlib import Path
from .cigsegy import Pysegy # type: ignore
from . import utils
from . import tools
import matplotlib.pyplot as plt


def plot_region(segy: Union[str, Pysegy],
                mode: str = 'line',
                loc: list = None,
                cdpxy_loc: list = None,
                save: str = None) -> None:
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
    if mode not in ['line', 'cdpxy']:
        raise RuntimeError(f"mode must be 'line' or 'cdpxy', mode = {mode}")

    assert loc is None or len(loc) == 4
    assert cdpxy_loc is None or len(cdpxy_loc) == 2
    cdpx = 73 if cdpxy_loc is None else cdpxy_loc[0]
    cdpy = 77 if cdpxy_loc is None else cdpxy_loc[1]

    if isinstance(segy, Pysegy):
        try:
            segy.scan()
        except:
            loc = utils.guess(segy)[0]
            segy.setInlineLocation(loc[0])
            segy.setCrosslineLocation(loc[1])
            segy.setSteps(loc[2], loc[3])
            segy.setXLocation(loc[4])
            segy.setYLocation(loc)[5]
            segy.scan()
        lineinfo = segy.get_lineInfo()
    elif isinstance(segy, (str, Path)):
        if loc is None:
            loc = utils.guess(str(segy))[0]
        segy = Pysegy(str(segy))
        segy.setInlineLocation(loc[0])
        segy.setCrosslineLocation(loc[1])
        segy.setSteps(loc[2], loc[3])
        segy.scan()
        lineinfo = segy.get_lineInfo()
    else:
        raise RuntimeError("Invalid type of `segy`")

    if mode == 'line':
        x = np.concatenate((lineinfo[:, 0], lineinfo[::-1, 0]))
        y = np.concatenate((lineinfo[:, 1], lineinfo[::-1, 2]))
        x = np.append(x, x[0])
        y = np.append(y, y[0])
    else:
        ni = lineinfo.shape[0]
        N = ni * 2 + 1
        x = np.zeros(N, dtype=int)
        y = np.zeros(N, dtype=int)
        for i in range(ni):
            x[i] = utils.get_trace_keys(segy, cdpx, lineinfo[i, 3])
            y[i] = utils.get_trace_keys(segy, cdpy, lineinfo[i, 3])
            x[N - i - 2] = utils.get_trace_keys(segy, cdpx, lineinfo[i, 4])
            y[N - i - 2] = utils.get_trace_keys(segy, cdpy, lineinfo[i, 4])
        x[-1] = x[0]
        y[-1] = y[0]

    istep = x[1] - x[0]
    xstep = (lineinfo[0, 2] - lineinfo[0, 1]) // (lineinfo[0, 5] - 1)

    plt.fill(x, y, color=(0.9, 0.9, 0.9))
    plt.plot(x, y)
    plt.gca().invert_yaxis()
    # plt.gca().xaxis.set_ticks_position('top')

    plt.grid(True, linestyle='--')
    if mode == 'line':
        xlabel = f"Inline Number/interval={istep}"
        ylabel = f"Crossline Number/interval={xstep}"
    else:
        xlabel = f"CDP X"
        ylabel = f"CDP Y"

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Region')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches='tight', pad_inches=0.0)
    plt.show()


def plot_trace_keys(segy: Union[str, Pysegy],
                    keyloc: int,
                    beg: int = 0,
                    end: int = 1000,
                    save: str = None) -> None:
    """
    plot the values (at keyloc in each trace) of the traces 
    range from beg to end .
    """
    assert beg >= 0 and end > beg, "invalid beg and end"
    assert end < tools.trace_count(segy), "end > trace_count"
    keys = utils.get_trace_keys(segy, keyloc, beg, end)
    x = np.arange(beg, end)
    plt.plot(x, keys)
    plt.xlabel("Trace Number")
    plt.ylabel(f"Key values (at {keyloc})")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches='tight', pad_inches=0.0)
    plt.show()


def plot_trace_ix(segy: Union[str, Pysegy],
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
    assert beg >= 0 and end > beg, "invalid beg and end"
    assert end < tools.trace_count(segy), "end > trace_count"
    keys = utils.get_trace_keys(segy, [iline, xline], beg, end)
    iv = keys[:, 0]
    xv = keys[:, 1]
    x = np.arange(beg, end)

    if figsize is not None:
        plt.figure(figsize=figsize)
    plt.subplot(1, 2, 1)
    plt.plot(x, iv)
    plt.xlabel('Trace Number')
    plt.ylabel(f"Inline (at {iline})")
    plt.subplot(1, 2, 2)
    plt.plot(x, xv)
    plt.xlabel('Trace Number')
    plt.ylabel(f"Crossline (at {xline})")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches='tight', pad_inches=0.0)
    plt.show()


def plot_trace_ixo(segy: str,
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
    assert beg >= 0 and end > beg, "invalid beg and end"
    assert end < tools.trace_count(segy), "end > trace_count"
    keys = utils.get_trace_keys(segy, [iline, xline, offset], beg, end)
    iv = keys[:, 0]
    xv = keys[:, 1]
    ov = keys[:, 2]
    x = np.arange(beg, end)

    if figsize is not None:
        plt.figure(figsize=figsize)
    plt.subplot(1, 3, 1)
    plt.plot(x, iv)
    plt.xlabel('Trace Number')
    plt.ylabel(f"Inline (at {iline})")
    plt.subplot(1, 3, 2)
    plt.plot(x, xv)
    plt.xlabel('Trace Number')
    plt.ylabel(f"Crossline (at {xline})")
    plt.subplot(1, 3, 3)
    plt.plot(x, ov)
    plt.xlabel('Trace Number')
    plt.ylabel(f"offset (at {offset})")
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=200, bbox_inches='tight', pad_inches=0.0)
    plt.show()
