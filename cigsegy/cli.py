# Copyright (c) 2026 Jintao Li.
# Zhejiang University (ZJU).
# All rights reserved.
"""Command line interface for cigsegy."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from . import SegyWriter
from .factories import collect, fromfile, metaInfo, textual_header, to_npy, tofile
from .tools import read_header


class _Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0

    try:
        return int(args.func(args) or 0)
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130
    except Exception as exc:
        if getattr(args, "verbose", False):
            raise
        print(f"cigsegy: error: {exc}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cigsegy",
        description="Inspect, read, and write SEG-Y files with cigsegy.",
        epilog="""Examples:
  cigsegy textual input.sgy
  cigsegy meta input.sgy
  cigsegy header input.sgy --trace 100
  cigsegy fromfile input.sgy volume.npy
  cigsegy to-npy input.sgy volume.npy
  cigsegy collect input.sgy traces.npy --beg 0 --end 1000
  cigsegy tofile input.sgy samples.dat
  cigsegy create template.sgy processed.npy processed.sgy --overwrite
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("--verbose", action="store_true", help="show Python traceback on error")

    subparsers = parser.add_subparsers(dest="command", metavar="command")
    _add_textual(subparsers)
    _add_meta(subparsers)
    _add_header(subparsers)
    _add_fromfile(subparsers)
    _add_to_npy(subparsers)
    _add_collect(subparsers)
    _add_tofile(subparsers)
    _add_create(subparsers)
    return parser


def _add_textual(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "textual",
        aliases=["textual-header", "text"],
        help="print the 3200-byte textual header",
        description="Print the SEG-Y textual header.",
        epilog="""Examples:
  cigsegy textual input.sgy
  cigsegy textual input.sgy --coding e --output textual.txt
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("segy", help="input SEG-Y file")
    parser.add_argument("--coding", choices=["u", "a", "e"], default="u", help="'u' guesses, 'a' is ASCII, 'e' is EBCDIC")
    parser.add_argument("--output", "-o", help="optional text output file")
    parser.set_defaults(func=_cmd_textual)


def _add_meta(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "meta",
        aliases=["info", "scan"],
        help="scan and print SEG-Y geometry metadata",
        description="Scan SEG-Y geometry metadata.",
        epilog="""Examples:
  cigsegy meta input.sgy
  cigsegy meta gather.sgy
  cigsegy meta input.sgy --iline 189 --xline 193
  cigsegy meta input.sgy --iline 189 --xline 193 --output meta.txt
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("segy", help="input SEG-Y file")
    _add_geometry_args(parser, include_xy=True)
    _add_dimensionality_args(parser)
    parser.add_argument("--no-scalar", action="store_true", help="do not apply coordinate/header scalars")
    parser.add_argument("--output", "-o", help="optional text output file")
    parser.set_defaults(func=_cmd_meta)


def _add_header(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "header",
        aliases=["hdr"],
        help="print binary or trace header",
        description="Print the binary header or one trace header.",
        epilog="""Examples:
  cigsegy header input.sgy --binary
  cigsegy header input.sgy --trace 100
  cigsegy header input.sgy --trace 100 --output trace100.json
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("segy", help="input SEG-Y file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--binary", "-b", action="store_true", help="print binary header")
    group.add_argument("--trace", "-t", type=int, help="print trace header by zero-based trace index")
    parser.add_argument("--output", "-o", help="optional JSON output file")
    parser.set_defaults(func=_cmd_header)


def _add_fromfile(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "fromfile",
        aliases=["tonpy"],
        help="read a regular SEG-Y volume and save it as .npy",
        description="Read a 3D or 4D SEG-Y volume into a NumPy .npy file.",
        epilog="""Examples:
  cigsegy fromfile input.sgy volume.npy
  cigsegy fromfile gather.sgy gathers.npy
  cigsegy fromfile input.sgy volume.npy --iline 189 --xline 193
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("segy", help="input SEG-Y file")
    parser.add_argument("output", help="output NumPy .npy file")
    _add_geometry_args(parser)
    _add_dimensionality_args(parser)
    parser.set_defaults(func=_cmd_fromfile)


def _add_collect(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "collect",
        help="collect traces and save them as .npy",
        description="Collect trace samples by file-order trace index and save them as a NumPy .npy file.",
        epilog="""Examples:
  cigsegy collect input.sgy traces.npy
  cigsegy collect input.sgy traces.npy --beg 0 --end 1000 --tbeg 100 --tend 800
  cigsegy collect input.sgy traces.npy --indices 0,10,20,30
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("segy", help="input SEG-Y file")
    parser.add_argument("output", help="output NumPy .npy file")
    parser.add_argument("--beg", type=int, default=-1, help="first trace index; <0 means all traces")
    parser.add_argument("--end", type=int, default=0, help="end trace index, exclusive; 0 means one trace at beg, <0 means file end")
    parser.add_argument("--tbeg", type=int, default=-1, help="first sample index; <0 means all samples")
    parser.add_argument("--tend", type=int, default=0, help="end sample index, exclusive; 0 means one sample at tbeg, <0 means trace end")
    parser.add_argument("--indices", help="comma-separated trace indices, or a .npy file containing 1D integer indices")
    parser.set_defaults(func=_cmd_collect)


def _add_to_npy(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "to-npy",
        aliases=["tonpy-stream"],
        help="stream SEG-Y samples directly to a .npy file",
        description="Write SEG-Y samples to a NumPy .npy file without loading the whole volume into memory.",
        epilog="""Examples:
  cigsegy to-npy input.sgy volume.npy
  cigsegy to-npy gather.sgy gathers.npy
  cigsegy to-npy input.sgy volume.npy --iline 189 --xline 193
  cigsegy to-npy input.sgy traces.npy --as2d
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("segy", help="input SEG-Y file")
    parser.add_argument("output", help="output NumPy .npy file")
    _add_geometry_args(parser)
    _add_dimensionality_args(parser)
    parser.add_argument("--as2d", action="store_true", help="ignore geometry and export traces as a 2D array")
    parser.set_defaults(func=_cmd_to_npy)


def _add_tofile(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "tofile",
        aliases=["toraw"],
        help="export SEG-Y samples to a raw binary file",
        description="Export SEG-Y samples to raw little-endian IEEE float32 data.",
        epilog="""Examples:
  cigsegy tofile input.sgy samples.dat
  cigsegy tofile input.sgy samples.dat --as2d
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("segy", help="input SEG-Y file")
    parser.add_argument("output", help="output raw binary file")
    _add_geometry_args(parser)
    _add_dimensionality_args(parser)
    parser.add_argument("--as2d", action="store_true", help="ignore geometry and export traces as a 2D array")
    parser.set_defaults(func=_cmd_tofile)


def _add_create(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "create",
        aliases=["write", "from-template"],
        help="write a SEG-Y from a template SEG-Y and NumPy/raw samples",
        description="Create a SEG-Y file by sharing geometry and headers from an existing SEG-Y.",
        epilog="""Examples:
  cigsegy create template.sgy processed.npy processed.sgy --overwrite
  cigsegy create template.sgy processed.dat processed.sgy --shape 100,200,300
  cigsegy create input_2ms.sgy output_1ms.npy output_1ms.sgy --sample-interval-us 1000 --non-strict
  cigsegy create input.sgy thin.npy thin.sgy --iline-slice ::2 --xline-slice ::2
""",
        formatter_class=_Formatter,
    )
    parser.add_argument("template", help="template SEG-Y file whose headers/geometry are reused")
    parser.add_argument("data", help="input .npy or raw float32 binary data")
    parser.add_argument("output", help="output SEG-Y file")
    parser.add_argument("--shape", help="data shape for raw binary input, e.g. 100,200,300 or 100x200x300")
    _add_geometry_args(parser)
    _add_dimensionality_args(parser)
    parser.add_argument("--keylocs", help="explicit keylocs as iline,xline,istep,xstep or iline,xline,offset,istep,xstep,ostep")
    parser.add_argument("--start", help="output start index for non-selected template writes, e.g. 0,0,0")
    parser.add_argument("--as2d", action="store_true", help="treat data as file-order trace/sample 2D data")
    parser.add_argument("--non-strict", action="store_true", help="allow changed time length/start/dt where supported")
    parser.add_argument("--start-time-from-zero", action="store_true", help="interpret time start from zero in non-strict template writes")
    parser.add_argument("--sample-interval-us", type=int, help="new sample interval in microseconds")
    parser.add_argument("--dt-new", type=int, default=0, help="legacy alias for --sample-interval-us; 0 keeps template dt")
    parser.add_argument("--textual", default="", help="optional textual header description or 3200-char header")
    parser.add_argument("--overwrite", action="store_true", help="overwrite output SEG-Y if it exists")
    parser.add_argument("--iline-slice", help="template iline selection, e.g. ::2 or 10:100:2")
    parser.add_argument("--xline-slice", help="template xline selection, e.g. ::2 or 10:100:2")
    parser.add_argument("--offset-slice", help="template offset selection for 4D data")
    parser.add_argument("--sample-slice", help="template sample selection, e.g. ::2")
    parser.add_argument("--trace-slice", help="file-order trace selection for 2D/template trace writes")
    parser.add_argument("--block-traces", type=int, default=16384, help="trace block size for selected template writes")
    parser.set_defaults(func=_cmd_create)


def _add_geometry_args(parser: argparse.ArgumentParser, *, include_xy: bool = False) -> None:
    parser.add_argument("--iline", type=int, help="optional inline byte-location override")
    parser.add_argument("--xline", type=int, help="optional crossline byte-location override")
    parser.add_argument("--offset", type=int, help="optional offset byte-location override")
    parser.add_argument("--istep", type=int, help="optional inline value-step override")
    parser.add_argument("--xstep", type=int, help="optional crossline value-step override")
    parser.add_argument("--ostep", type=int, help="optional offset value-step override")
    if include_xy:
        parser.add_argument("--xloc", type=int, help="X coordinate byte location in trace header")
        parser.add_argument("--yloc", type=int, help="Y coordinate byte location in trace header")


def _add_dimensionality_args(parser: argparse.ArgumentParser) -> None:
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--is4d", action="store_true", help="force 4D geometry")
    group.add_argument("--is3d", action="store_true", help="force 3D geometry")


def _cmd_textual(args: argparse.Namespace) -> int:
    text = textual_header(args.segy, coding=args.coding, printtext=args.output is None)
    if args.output is not None:
        _write_text(args.output, text)
    return 0


def _cmd_meta(args: argparse.Namespace) -> int:
    text = metaInfo(
        args.segy,
        iline=args.iline,
        xline=args.xline,
        offset=args.offset,
        istep=args.istep,
        xstep=args.xstep,
        ostep=args.ostep,
        xloc=args.xloc,
        yloc=args.yloc,
        is4d=_is4d_arg(args),
        apply_scalar=not args.no_scalar,
        printtext=args.output is None,
    )
    if args.output is not None:
        _write_text(args.output, text)
    return 0


def _cmd_header(args: argparse.Namespace) -> int:
    header_type = "th" if args.trace is not None else "bh"
    trace_index = 0 if args.trace is None else args.trace
    if args.output is None:
        read_header(args.segy, header_type, n=trace_index, printstr=True)
        return 0
    data = read_header(args.segy, header_type, n=trace_index, printstr=False)
    _write_json(args.output, data)
    return 0


def _cmd_fromfile(args: argparse.Namespace) -> int:
    data = fromfile(
        args.segy,
        iline=args.iline,
        xline=args.xline,
        offset=args.offset,
        istep=args.istep,
        xstep=args.xstep,
        ostep=args.ostep,
        is4d=_is4d_arg(args),
    )
    _save_npy(args.output, data)
    print(f"saved {args.output}: shape={data.shape}, dtype={data.dtype}")
    return 0


def _cmd_collect(args: argparse.Namespace) -> int:
    indices = None if args.indices is None else _indices_arg(args.indices)
    data = collect(
        args.segy,
        beg=args.beg,
        end=args.end,
        tbeg=args.tbeg,
        tend=args.tend,
        indices=indices,
    )
    _save_npy(args.output, data)
    print(f"saved {args.output}: shape={data.shape}, dtype={data.dtype}")
    return 0


def _cmd_to_npy(args: argparse.Namespace) -> int:
    shape = to_npy(
        args.segy,
        args.output,
        iline=args.iline,
        xline=args.xline,
        offset=args.offset,
        istep=args.istep,
        xstep=args.xstep,
        ostep=args.ostep,
        is4d=_is4d_arg(args),
        as2d=args.as2d,
    )
    print(f"saved {args.output}: shape={shape}, dtype=float32")
    return 0


def _cmd_tofile(args: argparse.Namespace) -> int:
    tofile(
        args.segy,
        args.output,
        iline=args.iline,
        xline=args.xline,
        offset=args.offset,
        istep=args.istep,
        xstep=args.xstep,
        ostep=args.ostep,
        is4d=_is4d_arg(args),
        as2d=args.as2d,
    )
    print(f"saved {args.output}")
    return 0


def _cmd_create(args: argparse.Namespace) -> int:
    shape = None if args.shape is None else _parse_int_tuple(args.shape, "shape")
    sample_interval = args.sample_interval_us
    if sample_interval is not None and args.dt_new:
        raise ValueError("--sample-interval-us and --dt-new cannot both be set")

    builder = SegyWriter.from_template(
        args.template,
        args.output,
        sample_interval_us=sample_interval,
        dt_new=args.dt_new,
        overwrite=args.overwrite,
    )
    _apply_create_options(builder, args)

    data: Any
    data_path = Path(args.data)
    if data_path.suffix.lower() == ".npy":
        data = np.load(data_path, mmap_mode="r")
        if shape is not None and tuple(shape) != tuple(data.shape):
            raise ValueError(f"--shape {shape} does not match npy shape {data.shape}")
    else:
        if shape is None:
            raise ValueError("--shape is required when input data is a raw binary file")
        data = str(data_path)

    with builder.open() as writer:
        writer.write(data, shape=shape)
    print(f"saved {args.output}")
    return 0


def _apply_create_options(builder: SegyWriter, args: argparse.Namespace) -> None:
    keylocs = _keylocs_arg(args)
    if keylocs is not None:
        if len(keylocs) == 4:
            builder.keylocs(iline=keylocs[0], xline=keylocs[1], istep=keylocs[2], xstep=keylocs[3])
        else:
            builder.keylocs(
                iline=keylocs[0],
                xline=keylocs[1],
                offset=keylocs[2],
                istep=keylocs[3],
                xstep=keylocs[4],
                ostep=keylocs[5],
            )

    if args.is4d:
        builder.as_4d()
    elif args.is3d:
        builder.as_3d()
    if args.as2d:
        builder.as2d()
    if args.non_strict:
        builder.strict(False)
    if args.start_time_from_zero:
        builder.start_time_from_zero(True)
    if args.start is not None:
        builder.start(_parse_int_tuple(args.start, "start"))
    if args.textual:
        builder.textual(args.textual)
    if args.block_traces:
        builder.block_traces(args.block_traces)

    selection = {
        "iline": _selection_arg(args.iline_slice, "iline-slice"),
        "xline": _selection_arg(args.xline_slice, "xline-slice"),
        "offset": _selection_arg(args.offset_slice, "offset-slice"),
        "sample": _selection_arg(args.sample_slice, "sample-slice"),
        "trace": _selection_arg(args.trace_slice, "trace-slice"),
    }
    if any(value is not None for value in selection.values()):
        builder.select(**selection)


def _keylocs_arg(args: argparse.Namespace) -> list[int] | None:
    if args.keylocs is not None:
        values = list(_parse_int_tuple(args.keylocs, "keylocs"))
        if len(values) not in (4, 6):
            raise ValueError("--keylocs must contain 4 or 6 integers")
        return values
    if args.iline is None and args.xline is None and args.offset is None:
        return None
    if args.iline is None or args.xline is None:
        raise ValueError("--iline and --xline must be provided together")
    if args.offset is None:
        return [args.iline, args.xline, args.istep or 1, args.xstep or 1]
    return [args.iline, args.xline, args.offset, args.istep or 1, args.xstep or 1, args.ostep or 1]


def _is4d_arg(args: argparse.Namespace) -> bool | None:
    if getattr(args, "is4d", False):
        return True
    if getattr(args, "is3d", False):
        return False
    return None


def _selection_arg(value: str | None, name: str) -> Any:
    if value is None:
        return None
    value = value.strip()
    if not value:
        raise ValueError(f"--{name} cannot be empty")
    if ":" in value:
        return _parse_slice(value, name)
    values = _parse_int_tuple(value, name)
    if len(values) == 1:
        return values[0]
    return np.asarray(values, dtype=np.int64)


def _parse_slice(value: str, name: str) -> slice:
    parts = value.split(":")
    if len(parts) > 3:
        raise ValueError(f"--{name} slice must be start:stop[:step]")
    parsed = []
    for part in parts:
        parsed.append(None if part == "" else int(part))
    while len(parsed) < 3:
        parsed.append(None)
    return slice(parsed[0], parsed[1], parsed[2])


def _parse_int_tuple(value: str | Iterable[int], name: str) -> tuple[int, ...]:
    if isinstance(value, str):
        text = value.replace("x", ",").replace("X", ",")
        parts = [part.strip() for part in text.split(",") if part.strip()]
        if not parts:
            raise ValueError(f"--{name} must contain at least one integer")
        return tuple(int(part) for part in parts)
    return tuple(int(v) for v in value)


def _indices_arg(value: str) -> np.ndarray:
    path = Path(value)
    if path.suffix.lower() == ".npy" and path.exists():
        out = np.load(path)
    else:
        out = np.asarray(_parse_int_tuple(value, "indices"), dtype=np.int64)
    if out.ndim != 1:
        raise ValueError("--indices must be a 1D integer array")
    if not np.issubdtype(out.dtype, np.integer):
        raise ValueError("--indices must contain integers")
    return np.ascontiguousarray(out, dtype=np.int32)


def _save_npy(path: str, data: np.ndarray) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as file:
        np.save(file, data)


def _write_text(path: str, text: str) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")


def _write_json(path: str, data: Any) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=_json_default), encoding="utf-8")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


if __name__ == "__main__":
    raise SystemExit(main())
