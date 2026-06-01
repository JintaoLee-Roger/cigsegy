# Copyright (c) 2026 Jintao Li.
# Zhejiang University (ZJU).
# All rights reserved.
"""Incremental SEG-Y writer API.

``SegyWriter`` is an additive API.  It does not change the behavior of
``create_by_sharing_header``; template writes delegate to that stable path.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from . import utils
from .cpp import _CXX_SEGY
from .createtool import assemble_metainfo, generate_textual
from .factories import create_by_sharing_header


TEXTUAL_HEADER_SIZE = 3200
BINARY_HEADER_SIZE = 400
TRACE_HEADER_SIZE = 240
STANDARD_FIRST_TRACE_OFFSET = TEXTUAL_HEADER_SIZE + BINARY_HEADER_SIZE
SUPPORTED_SAMPLE_BYTES = {
    1: 4,
    2: 4,
    3: 2,
    5: 4,
    8: 1,
    10: 4,
    11: 2,
    16: 1,
}


@dataclass
class _WriterConfig:
    mode: str
    out_segy: Path
    overwrite: bool = False
    header_segy: Path | None = None
    textual: bytes | str | list[str] | None = None
    binary: bytes | bytearray | np.ndarray | None = None
    extended_textual: bytes | bytearray | np.ndarray | None = None
    data_trailer: bytes | bytearray | np.ndarray | None = None
    trace_headers: np.ndarray | None = None
    shape: tuple[int, ...] | None = None
    sample_format: int | None = None
    sample_count: int | None = None
    sample_interval_us: int = 2000
    sample_interval_override: int | None = None
    start_time_ms: int = 0
    start: tuple[int, ...] | None = None
    interval: tuple[int, ...] | None = None
    iline_start: int = 1
    xline_start: int = 1
    offset_start: int = 1
    iline_step: int = 1
    xline_step: int = 1
    offset_step: int = 1
    x_start: int = 0
    y_start: int = 0
    x_step: int = 1
    y_step: int = 1
    geometry: dict[str, Any] | None = None
    keylocs_value: list[int] | None = None
    is4d: bool | None = None
    as2d: bool = False
    strict_value: bool = True
    start_time_from_zero_value: bool = False
    dt_new_value: int = 0
    block_traces_value: int = 16384
    selection: "_TemplateSelection | None" = None
    attrs: dict[str, Any] = field(default_factory=dict)


@dataclass
class _TemplateSelection:
    iline: Any = None
    xline: Any = None
    offset: Any = None
    sample: Any = None
    trace: Any = None
    sample_is_explicit: bool = False


@dataclass
class _AxisSelection:
    indices: np.ndarray
    explicit: bool = False

    @property
    def size(self) -> int:
        return int(self.indices.size)


class SegyWriter:
    """Builder for writing SEG-Y files.

    Examples
    --------
    Configure first, then open the writer:

    >>> b = SegyWriter.from_template("orig.sgy", "out.sgy")
    >>> b.keylocs(iline=189, xline=193).overwrite(True)
    >>> with b.open() as w:
    ...     w.write(data)

    Write from stored file-order headers:

    >>> b = SegyWriter.from_headers("out.sgy")
    >>> b.textual(textual).binary(binary).sample_format(1)
    >>> with b.open() as w:
    ...     for trace_headers, samples in blocks:
    ...         w.write_trace_block(trace_headers, samples)

    Create minimal headers from scratch:

    >>> b = SegyWriter.create("out.sgy", shape=(10, 20, 300))
    >>> b.sample_format(5).sample_interval_us(2000)
    >>> with b.open() as w:
    ...     w.write(data)
    """

    def __init__(self, config: _WriterConfig):
        self._config = config
        self._active: _OpenSegyWriter | None = None

    @classmethod
    def from_template(
        cls,
        header_segy: str | PathLike[str],
        out_segy: str | PathLike[str],
        *,
        sample_interval_us: int | None = None,
        dt_new: int = 0,
        overwrite: bool = False,
    ) -> "SegyWriter":
        dt_value = _resolve_dt_override(sample_interval_us, dt_new)
        return cls(
            _WriterConfig(
                mode="template",
                header_segy=Path(header_segy),
                out_segy=Path(out_segy),
                dt_new_value=dt_value,
                sample_interval_us=dt_value if dt_value > 0 else 2000,
                overwrite=overwrite,
            )
        )

    @classmethod
    def template(
        cls,
        header_segy: str | PathLike[str],
        out_segy: str | PathLike[str],
        *,
        sample_interval_us: int | None = None,
        dt_new: int = 0,
        overwrite: bool = False,
    ) -> "SegyWriter":
        return cls.from_template(
            header_segy,
            out_segy,
            sample_interval_us=sample_interval_us,
            dt_new=dt_new,
            overwrite=overwrite,
        )

    @classmethod
    def from_headers(
        cls,
        out_segy: str | PathLike[str],
        *,
        textual: bytes | bytearray | np.ndarray | str | None = None,
        binary: bytes | bytearray | np.ndarray | None = None,
        extended_textual: bytes | bytearray | np.ndarray | None = None,
        trace_headers: np.ndarray | None = None,
        sample_format: int | None = None,
        sample_count: int | None = None,
        sample_interval_us: int | None = None,
        overwrite: bool = False,
    ) -> "SegyWriter":
        return cls(
            _WriterConfig(
                mode="headers",
                out_segy=Path(out_segy),
                textual=textual,
                binary=binary,
                extended_textual=extended_textual,
                trace_headers=trace_headers,
                sample_format=sample_format,
                sample_count=sample_count,
                sample_interval_override=(
                    None if sample_interval_us is None else _sample_interval_value(sample_interval_us)
                ),
                overwrite=overwrite,
            )
        )

    @classmethod
    def headers(
        cls,
        out_segy: str | PathLike[str],
        **kwargs: Any,
    ) -> "SegyWriter":
        return cls.from_headers(out_segy, **kwargs)

    @classmethod
    def create(
        cls,
        out_segy: str | PathLike[str],
        *,
        shape: Iterable[int] | None = None,
        sample_format: int = 5,
        sample_interval_us: int = 2000,
        start_time_ms: int = 0,
        overwrite: bool = False,
    ) -> "SegyWriter":
        shape_tuple = None if shape is None else _as_int_tuple(shape, "shape")
        return cls(
            _WriterConfig(
                mode="create",
                out_segy=Path(out_segy),
                shape=shape_tuple,
                sample_format=int(sample_format),
                sample_interval_us=int(sample_interval_us),
                start_time_ms=int(start_time_ms),
                overwrite=overwrite,
            )
        )

    def as2d(self, value: bool = True) -> "SegyWriter":
        self._config.as2d = bool(value)
        return self

    def as_2d(self, value: bool = True) -> "SegyWriter":
        return self.as2d(value)

    def as_3d(self) -> "SegyWriter":
        self._config.as2d = False
        self._config.is4d = False
        return self

    def as_4d(self) -> "SegyWriter":
        self._config.as2d = False
        self._config.is4d = True
        return self

    def binary(self, value: bytes | bytearray | np.ndarray) -> "SegyWriter":
        self._config.binary = value
        return self

    def block_traces(self, value: int) -> "SegyWriter":
        self._config.block_traces_value = max(1, int(value))
        return self

    def data_trailer(self, value: bytes | bytearray | np.ndarray) -> "SegyWriter":
        self._config.data_trailer = value
        return self

    def extended_textual(self, value: bytes | bytearray | np.ndarray) -> "SegyWriter":
        self._config.extended_textual = value
        return self

    def geometry(self, **kwargs: Any) -> "SegyWriter":
        geom = dict(self._config.geometry or {})
        geom.update(kwargs)
        self._config.geometry = geom
        return self

    def grid(
        self,
        *,
        iline_start: int | None = None,
        xline_start: int | None = None,
        offset_start: int | None = None,
        iline_step: int | None = None,
        xline_step: int | None = None,
        offset_step: int | None = None,
    ) -> "SegyWriter":
        if iline_start is not None:
            self._config.iline_start = int(iline_start)
        if xline_start is not None:
            self._config.xline_start = int(xline_start)
        if offset_start is not None:
            self._config.offset_start = int(offset_start)
        if iline_step is not None:
            self._config.iline_step = int(iline_step)
        if xline_step is not None:
            self._config.xline_step = int(xline_step)
        if offset_step is not None:
            self._config.offset_step = int(offset_step)
        return self

    def keylocs(
        self,
        *,
        iline: int | None = None,
        xline: int | None = None,
        offset: int | None = None,
        istep: int = 1,
        xstep: int = 1,
        ostep: int = 1,
    ) -> "SegyWriter":
        if iline is None or xline is None:
            self._config.keylocs_value = None
        elif offset is None:
            self._config.keylocs_value = [int(iline), int(xline), int(istep), int(xstep)]
        else:
            self._config.keylocs_value = [
                int(iline),
                int(xline),
                int(offset),
                int(istep),
                int(xstep),
                int(ostep),
            ]
        return self

    def origin(
        self,
        *,
        x_start: int | None = None,
        y_start: int | None = None,
        x_step: int | None = None,
        y_step: int | None = None,
    ) -> "SegyWriter":
        if x_start is not None:
            self._config.x_start = int(x_start)
        if y_start is not None:
            self._config.y_start = int(y_start)
        if x_step is not None:
            self._config.x_step = int(x_step)
        if y_step is not None:
            self._config.y_step = int(y_step)
        return self

    def overwrite(self, value: bool = True) -> "SegyWriter":
        self._config.overwrite = bool(value)
        return self

    def sample_count(self, value: int) -> "SegyWriter":
        self._config.sample_count = int(value)
        return self

    def sample_format(self, value: int) -> "SegyWriter":
        self._config.sample_format = int(value)
        return self

    def sample_interval_us(self, value: int) -> "SegyWriter":
        value = _sample_interval_value(value)
        self._config.sample_interval_us = value
        if self._config.mode == "template":
            self._config.dt_new_value = value
        elif self._config.mode == "headers":
            self._config.sample_interval_override = value
        return self

    def select(
        self,
        *,
        iline: Any = None,
        xline: Any = None,
        offset: Any = None,
        sample: Any = None,
        trace: Any = None,
    ) -> "SegyWriter":
        if self._config.mode != "template":
            raise ValueError("select() is only available for from_template mode")
        if trace is not None and any(axis is not None for axis in (iline, xline, offset)):
            raise ValueError("trace selection cannot be combined with iline/xline/offset selection")
        self._config.selection = _TemplateSelection(
            iline=iline,
            xline=xline,
            offset=offset,
            sample=sample,
            trace=trace,
            sample_is_explicit=sample is not None,
        )
        return self

    def spatial_stride(self, *, iline: int = 1, xline: int = 1, offset: int = 1) -> "SegyWriter":
        if iline <= 0 or xline <= 0 or offset <= 0:
            raise ValueError("spatial stride values must be positive")
        selection = self._config.selection or _TemplateSelection()
        selection.iline = slice(None, None, int(iline))
        selection.xline = slice(None, None, int(xline))
        selection.offset = slice(None, None, int(offset))
        self._config.selection = selection
        return self

    def dt_new(self, value: int) -> "SegyWriter":
        value = _dt_new_value(value)
        self._config.dt_new_value = value
        if self._config.mode == "template" and value > 0:
            self._config.sample_interval_us = value
        elif self._config.mode == "headers":
            self._config.sample_interval_override = None if value == 0 else value
        elif self._config.mode == "create" and value > 0:
            self._config.sample_interval_us = value
        return self

    def shape(self, *value: int | Iterable[int]) -> "SegyWriter":
        if len(value) == 1 and not isinstance(value[0], (int, np.integer)):
            self._config.shape = _as_int_tuple(value[0], "shape")
        else:
            self._config.shape = _as_int_tuple(value, "shape")
        return self

    def start(self, *value: int | Iterable[int]) -> "SegyWriter":
        if len(value) == 1 and not isinstance(value[0], (int, np.integer)):
            self._config.start = _as_int_tuple(value[0], "start")
        else:
            self._config.start = _as_int_tuple(value, "start")
        return self

    def start_time_from_zero(self, value: bool = True) -> "SegyWriter":
        self._config.start_time_from_zero_value = bool(value)
        return self

    def start_time_ms(self, value: int) -> "SegyWriter":
        self._config.start_time_ms = int(value)
        return self

    def strict(self, value: bool = True) -> "SegyWriter":
        self._config.strict_value = bool(value)
        return self

    def textual(self, value: bytes | bytearray | np.ndarray | str | list[str]) -> "SegyWriter":
        self._config.textual = value
        return self

    def trace_headers(self, value: np.ndarray) -> "SegyWriter":
        self._config.trace_headers = _trace_headers_array(value)
        return self

    def open(self) -> "_OpenSegyWriter":
        if self._active is not None and not self._active.closed:
            raise RuntimeError("this SegyWriter builder already has an active writer")
        self._active = _OpenSegyWriter(self._config)
        return self._active

    def __enter__(self) -> "_OpenSegyWriter":
        return self.open().__enter__()

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> bool:
        if self._active is None:
            return False
        return self._active.__exit__(exc_type, exc, tb)


class _OpenSegyWriter:
    def __init__(self, config: _WriterConfig):
        self._config = config
        self._backend: Any | None = None
        self._closed = False
        self._finalized = False
        self._trace_count = 0
        self._sample_count = config.sample_count
        self._sample_format = config.sample_format
        self._textual: bytes | None = None
        self._binary: bytearray | None = None
        self._extended_textual: bytes = _bytes_from_optional(config.extended_textual)
        self._trace_headers = None if config.trace_headers is None else _trace_headers_array(config.trace_headers)
        if config.mode in {"headers", "create"}:
            self._open_file()

    @property
    def closed(self) -> bool:
        return self._closed

    @property
    def trace_count(self) -> int:
        if self._backend is not None:
            return int(self._backend.trace_count)
        return self._trace_count

    def close(self) -> None:
        if self._closed:
            return
        if self._backend is not None:
            self._backend.close()
            self._backend = None
        self._closed = True

    def finalize(self, data_trailer: bytes | bytearray | np.ndarray | None = None) -> None:
        if self._finalized:
            return
        if self._config.mode == "template":
            self._finalized = True
            self._closed = True
            return
        self._require_open_file()
        trailer = data_trailer if data_trailer is not None else self._config.data_trailer
        self._backend.finalize(_uint8_array(_bytes_from_optional(trailer)))
        self._trace_count = int(self._backend.trace_count)
        self._finalized = True
        self.close()

    def write(self, data: Any, *, start: Iterable[int] | None = None, shape: Iterable[int] | None = None) -> None:
        if self._config.mode == "template":
            self._write_template(data, start=start, shape=shape)
            return

        array = np.asarray(data)
        if array.ndim < 2:
            raise ValueError("data must have at least 2 dimensions: (..., sample)")
        if self._sample_count is None:
            self._sample_count = int(array.shape[-1])
        if int(array.shape[-1]) != int(self._sample_count):
            raise ValueError(f"data sample count {array.shape[-1]} does not match {self._sample_count}")
        samples = array.reshape((-1, int(self._sample_count)))

        if self._config.mode == "headers":
            if self._trace_headers is None:
                raise ValueError("write(data) in from_headers mode requires builder.trace_headers(...); use write_trace_block otherwise")
            if self._trace_headers.shape[0] != samples.shape[0]:
                raise ValueError(
                    f"trace_headers rows ({self._trace_headers.shape[0]}) must match trace count ({samples.shape[0]})"
                )
            self.write_trace_block(self._trace_headers, samples)
            return

        headers = self._make_trace_headers_for_block(array.shape, start=start)
        self.write_trace_block(headers, samples)

    def write_block(self, data: Any, *, start: Iterable[int] | None = None) -> None:
        if self._config.mode == "template":
            raise NotImplementedError(
                "template write_block is not implemented yet; use write(...) or from_headers(...).write_trace_block(...)"
            )
        self.write(data, start=start)

    def write_file(
        self,
        src: str | PathLike[str],
        *,
        shape: Iterable[int] | None = None,
        dtype: np.dtype[Any] | str = "float32",
        order: str = "C",
        start: Iterable[int] | None = None,
    ) -> None:
        path = Path(src)
        if self._config.mode == "template":
            self._write_template(str(path), start=start, shape=shape)
            return
        if path.suffix.lower() == ".npy":
            data = np.load(path, mmap_mode="r")
        else:
            if shape is None:
                raise ValueError("shape is required when writing a raw binary file")
            data = np.memmap(path, dtype=np.dtype(dtype), mode="r", shape=tuple(shape), order=order)
        self.write(data, start=start)

    def write_raw_trace_block(self, trace_headers: Any, sample_bytes: Any) -> None:
        self._require_open_file()
        headers = _trace_headers_array(trace_headers)
        headers = _trace_headers_with_sample_interval(headers, self._config.sample_interval_override)
        raw = np.ascontiguousarray(np.asarray(sample_bytes, dtype=np.uint8))
        self._backend.write_raw_trace_block(headers, raw)
        self._trace_count = int(self._backend.trace_count)
        self._sample_count = int(self._backend.sample_count)
        self._sample_format = int(self._backend.sample_format)

    def write_trace_block(self, trace_headers: Any, samples: Any) -> None:
        self._require_open_file()
        headers = _trace_headers_array(trace_headers)
        headers = _trace_headers_with_sample_interval(headers, self._config.sample_interval_override)
        data = np.asarray(samples)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.ndim != 2:
            raise ValueError("samples must be a 2D array with shape (ntrace, nsample)")
        if headers.shape[0] != data.shape[0]:
            raise ValueError(f"trace_headers rows ({headers.shape[0]}) must match samples rows ({data.shape[0]})")
        if self._sample_count is None:
            self._sample_count = int(data.shape[1])
        if data.shape[1] != int(self._sample_count):
            raise ValueError(f"samples shape[1] ({data.shape[1]}) must match sample_count ({self._sample_count})")
        data = np.ascontiguousarray(data, dtype=np.float32)
        self._backend.write_trace_block(headers, data)
        self._trace_count = int(self._backend.trace_count)
        self._sample_count = int(self._backend.sample_count)
        self._sample_format = int(self._backend.sample_format)

    def __enter__(self) -> "_OpenSegyWriter":
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> bool:
        if exc_type is None:
            self.finalize()
        else:
            self.close()
        return False

    def _open_file(self) -> None:
        out = self._config.out_segy
        if out.exists():
            if not self._config.overwrite:
                raise FileExistsError(out)
        out.parent.mkdir(parents=True, exist_ok=True)

        if self._config.mode == "headers":
            self._textual = _strict_textual_bytes(self._config.textual)
            self._binary = bytearray(_strict_binary_bytes(self._config.binary))
            if self._config.sample_interval_override is not None:
                _set_i2(self._binary, 17, int(self._config.sample_interval_override))
            if self._sample_format is None:
                self._sample_format = _read_i2(self._binary, 25)
            else:
                _set_i2(self._binary, 25, int(self._sample_format))
            if self._sample_count is None:
                binary_sample_count = _read_i2(self._binary, 21)
                self._sample_count = binary_sample_count if binary_sample_count > 0 else None
            else:
                _set_i2(self._binary, 21, int(self._sample_count))
        else:
            self._textual, self._binary = _created_headers(self._config)
            self._sample_format = self._config.sample_format
            self._sample_count = self._config.shape[-1] if self._config.shape is not None else None

        _validate_extended_textual(self._extended_textual)
        self._backend = _CXX_SEGY.SegyBlockWriter(
            str(out),
            _uint8_array(self._textual),
            _uint8_array(bytes(self._binary)),
            _uint8_array(self._extended_textual),
            int(self._sample_format or 0),
            int(self._sample_count or 0),
            bool(self._config.overwrite),
        )

    def _write_template(self, data: Any, *, start: Iterable[int] | None, shape: Iterable[int] | None) -> None:
        if self._finalized:
            raise RuntimeError("writer is already finalized")
        if self._config.selection is not None:
            self._write_template_selection(data, start=start, shape=shape)
            return
        out = self._config.out_segy
        if out.exists():
            if not self._config.overwrite:
                raise FileExistsError(out)
            out.unlink()
        if self._config.header_segy is None:
            raise ValueError("template mode requires header_segy")

        start_list = _optional_int_list(start if start is not None else self._config.start)
        shape_tuple = None if shape is None else tuple(int(v) for v in shape)
        create_by_sharing_header(
            str(out),
            str(self._config.header_segy),
            data,
            shape=shape_tuple,
            start=start_list,
            keylocs=self._config.keylocs_value,
            is4d=self._config.is4d,
            as2d=self._config.as2d,
            textual=self._config.textual or "",
            strict=self._config.strict_value,
            start_time_from_zero=self._config.start_time_from_zero_value,
            dt_new=self._config.dt_new_value,
        )
        self._finalized = True
        self._closed = True

    def _write_template_selection(
        self,
        data: Any,
        *,
        start: Iterable[int] | None,
        shape: Iterable[int] | None,
    ) -> None:
        if start is not None:
            raise ValueError("start cannot be used together with select(); put the spatial/time origin in select()")
        if self._config.header_segy is None:
            raise ValueError("template mode requires header_segy")

        array = _array_from_data(data, shape=shape)
        if array.ndim < 2:
            raise ValueError("selected template writes require data with at least 2 dimensions")

        out = self._config.out_segy
        if out.exists():
            if not self._config.overwrite:
                raise FileExistsError(out)
            out.unlink()
        out.parent.mkdir(parents=True, exist_ok=True)

        template = _open_template_segy(self._config)
        try:
            source_shape = _template_shape(template, self._config)
            selection = self._config.selection
            if selection is None:
                raise RuntimeError("template selection is not configured")
            trace_indices, spatial_shape = _selected_trace_indices(template, source_shape, selection, self._config)
            sample_axis = _selected_sample_axis(source_shape[-1], selection)
            if sample_axis.explicit and array.shape[-1] != sample_axis.size:
                raise ValueError(
                    f"data sample count {array.shape[-1]} does not match selected sample count {sample_axis.size}"
                )
            expected_shape = tuple(spatial_shape) + (int(array.shape[-1]), )
            if tuple(array.shape) != expected_shape:
                raise ValueError(f"data shape {array.shape} does not match selected output shape {expected_shape}")

            textual, binary, extended_textual = _template_file_headers(self._config.header_segy, self._config.textual)
            sample_count = int(array.shape[-1])
            sample_interval = _template_output_sample_interval(template, sample_axis, self._config)
            start_time = _template_output_start_time(template, sample_axis)
            _set_i2(binary, 17, sample_interval)
            _set_i2(binary, 21, sample_count)
            _set_i2(binary, 305, len(extended_textual) // TEXTUAL_HEADER_SIZE)
            _set_u8(binary, 321, STANDARD_FIRST_TRACE_OFFSET + len(extended_textual))

            sample_format = _read_i2(binary, 25)
            samples = np.asarray(array).reshape((-1, sample_count))
            self._backend = _CXX_SEGY.SegyBlockWriter(
                str(out),
                _uint8_array(textual),
                _uint8_array(bytes(binary)),
                _uint8_array(extended_textual),
                int(sample_format),
                sample_count,
                bool(self._config.overwrite),
            )
            block_traces = int(self._config.block_traces_value)
            for beg in range(0, trace_indices.size, block_traces):
                end = min(beg + block_traces, trace_indices.size)
                headers = _template_trace_headers(
                    template,
                    trace_indices[beg:end],
                    sample_count=sample_count,
                    sample_interval_us=sample_interval,
                    start_time_ms=start_time,
                )
                self._backend.write_trace_block(headers, np.ascontiguousarray(samples[beg:end], dtype=np.float32))
            self._backend.finalize(_uint8_array(_bytes_from_optional(self._config.data_trailer)))
            self._trace_count = int(self._backend.trace_count)
            self._sample_count = int(self._backend.sample_count)
            self._sample_format = int(self._backend.sample_format)
            self._finalized = True
        finally:
            template.close()
            self.close()

    def _make_trace_headers_for_block(self, shape: tuple[int, ...], *, start: Iterable[int] | None) -> np.ndarray:
        if self._config.mode != "create":
            raise RuntimeError("trace header generation is only available in create mode")
        if len(shape) < 2 or len(shape) > 4:
            raise ValueError(f"create mode supports 2D, 3D, or 4D arrays, got shape={shape}")
        full_shape = self._config.shape or tuple(int(v) for v in shape)
        block_start = tuple(0 for _ in shape) if start is None else _as_int_tuple(start, "start")
        if len(block_start) != len(shape):
            raise ValueError("start length must match data ndim")
        return _generate_trace_headers(self._config, shape, full_shape, block_start, self._trace_count)

    def _require_open_file(self) -> None:
        if self._closed or self._backend is None:
            raise RuntimeError("writer is closed")
        if self._finalized:
            raise RuntimeError("writer is already finalized")


def _as_int_tuple(value: Iterable[int], name: str) -> tuple[int, ...]:
    try:
        out = tuple(int(v) for v in value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an iterable of integers") from exc
    if not out:
        raise ValueError(f"{name} must not be empty")
    if any(v < 0 for v in out):
        raise ValueError(f"{name} values must be non-negative")
    return out


def _array_from_data(data: Any, *, shape: Iterable[int] | None) -> np.ndarray:
    if isinstance(data, (str, PathLike)):
        path = Path(data)
        if path.suffix.lower() == ".npy":
            return np.load(path, mmap_mode="r")
        if shape is None:
            raise ValueError("shape is required when writing a raw binary file with select()")
        return np.memmap(path, dtype=np.float32, mode="r", shape=tuple(int(v) for v in shape), order="C")
    return np.asarray(data)


def _open_template_segy(config: _WriterConfig) -> Any:
    if config.header_segy is None:
        raise ValueError("template mode requires header_segy")
    segy = _CXX_SEGY.Pysegy(str(config.header_segy))
    if config.as2d:
        return segy

    if config.keylocs_value is None:
        keylocs = utils.guess(segy)
        iline, xline, offset, istep, xstep, ostep, xloc, yloc, guessed_is4d = keylocs
        ndim = 4 if (config.is4d if config.is4d is not None else guessed_is4d) else 3
    elif len(config.keylocs_value) == 4:
        if config.is4d:
            segy.close()
            raise ValueError("4D template selection requires offset key location in keylocs()")
        iline, xline, istep, xstep = config.keylocs_value
        offset, ostep, xloc, yloc, ndim = 37, 1, 181, 185, 3
    elif len(config.keylocs_value) == 6:
        iline, xline, offset, istep, xstep, ostep = config.keylocs_value
        xloc, yloc = 181, 185
        ndim = 4 if config.is4d is not False else 3
    else:
        segy.close()
        raise ValueError("keylocs must contain 4 or 6 values")

    segy.setLocations(int(iline), int(xline), int(offset))
    segy.setSteps(int(istep), int(xstep), int(ostep))
    segy.setXYLocations(int(xloc), int(yloc))
    segy.set_segy_type(int(ndim))
    segy.scan()
    return segy


def _template_shape(segy: Any, config: _WriterConfig) -> tuple[int, ...]:
    if config.as2d:
        return int(segy.ntrace), int(segy.nt)
    return tuple(int(v) for v in segy.shape)


def _selected_trace_indices(
    segy: Any,
    source_shape: tuple[int, ...],
    selection: _TemplateSelection,
    config: _WriterConfig,
) -> tuple[np.ndarray, tuple[int, ...]]:
    if config.as2d or len(source_shape) == 2 or selection.trace is not None:
        if any(axis is not None for axis in (selection.iline, selection.xline, selection.offset)):
            raise ValueError("2D/trace template selection cannot use iline/xline/offset axes")
        trace_axis = _normalize_axis_selection(selection.trace, source_shape[0], "trace")
        if trace_axis.size == 0:
            raise ValueError("trace selection is empty")
        return trace_axis.indices.astype(np.int32, copy=False), (trace_axis.size, )

    spatial_shape = source_shape[:-1]
    if len(spatial_shape) not in (2, 3):
        raise ValueError(f"template selection supports 2D, 3D, or 4D template shapes, got {source_shape}")
    expected_ntrace = int(np.prod(spatial_shape))
    if expected_ntrace != int(segy.ntrace):
        raise RuntimeError(
            "template selection requires a regular template where product(shape[:-1]) equals ntrace; "
            f"got product={expected_ntrace}, ntrace={segy.ntrace}"
        )

    axes = [
        _normalize_axis_selection(selection.iline, spatial_shape[0], "iline"),
        _normalize_axis_selection(selection.xline, spatial_shape[1], "xline"),
    ]
    if len(spatial_shape) == 3:
        axes.append(_normalize_axis_selection(selection.offset, spatial_shape[2], "offset"))
    elif selection.offset is not None:
        raise ValueError("offset selection requires a 4D template")
    if any(axis.size == 0 for axis in axes):
        raise ValueError("spatial selection is empty")

    grids = np.meshgrid(*[axis.indices for axis in axes], indexing="ij")
    flat_coords = [grid.reshape(-1) for grid in grids]
    trace_indices = np.ravel_multi_index(flat_coords, spatial_shape).astype(np.int32)
    _validate_selected_template_headers(segy, trace_indices, flat_coords)
    return trace_indices, tuple(axis.size for axis in axes)


def _normalize_axis_selection(value: Any, axis_size: int, name: str) -> _AxisSelection:
    if value is None:
        return _AxisSelection(np.arange(axis_size, dtype=np.int64), explicit=False)
    if isinstance(value, (int, np.integer)):
        index = int(value)
        if index < 0:
            index += axis_size
        if index < 0 or index >= axis_size:
            raise IndexError(f"{name} index {value} is out of bounds for size {axis_size}")
        return _AxisSelection(np.array([index], dtype=np.int64), explicit=True)
    if isinstance(value, slice):
        start, stop, step = value.indices(axis_size)
        if step <= 0:
            raise ValueError(f"{name} selection requires a positive step")
        return _AxisSelection(np.arange(start, stop, step, dtype=np.int64), explicit=True)
    arr = np.asarray(value)
    if arr.ndim == 0:
        return _normalize_axis_selection(arr.item(), axis_size, name)
    if arr.ndim != 1:
        raise IndexError(f"{name} selection must be a scalar, slice, or 1D array")
    if arr.dtype == np.bool_:
        if arr.size != axis_size:
            raise IndexError(f"{name} boolean selection size {arr.size} does not match axis size {axis_size}")
        arr = np.flatnonzero(arr)
    elif not np.issubdtype(arr.dtype, np.integer):
        raise IndexError(f"{name} selection array must be integer or boolean")
    arr = arr.astype(np.int64, copy=False)
    arr = np.where(arr < 0, arr + axis_size, arr)
    if arr.size and ((arr < 0).any() or (arr >= axis_size).any()):
        raise IndexError(f"{name} selection contains indices outside [0, {axis_size})")
    return _AxisSelection(np.ascontiguousarray(arr, dtype=np.int64), explicit=True)


def _selected_sample_axis(source_nt: int, selection: _TemplateSelection) -> _AxisSelection:
    sample_axis = _normalize_axis_selection(selection.sample, source_nt, "sample")
    sample_axis.explicit = selection.sample_is_explicit
    if sample_axis.explicit:
        if sample_axis.size == 0:
            raise ValueError("sample selection is empty")
        _axis_positive_constant_step(sample_axis.indices, "sample")
    return sample_axis


def _axis_positive_constant_step(indices: np.ndarray, name: str) -> int:
    if indices.size <= 1:
        return 1
    steps = np.diff(indices)
    if not np.all(steps == steps[0]) or int(steps[0]) <= 0:
        raise ValueError(f"{name} selection must have a positive constant step")
    return int(steps[0])


def _validate_selected_template_headers(segy: Any, trace_indices: np.ndarray, flat_coords: list[np.ndarray]) -> None:
    keylocs = segy.get_keylocs()
    meta = segy.get_metainfo()
    keys = [int(keylocs["iline"]), int(keylocs["xline"])]
    expected = [
        int(meta["start_iline"]) + flat_coords[0].astype(np.int64) * int(keylocs["istep"]),
        int(meta["start_xline"]) + flat_coords[1].astype(np.int64) * int(keylocs["xstep"]),
    ]
    if len(flat_coords) == 3:
        keys.append(int(keylocs["offset"]))
        expected.append(int(meta["start_offset"]) + flat_coords[2].astype(np.int64) * int(keylocs["ostep"]))
    values = segy.get_trace_keys(keys, [4] * len(keys), trace_indices.astype(np.int32, copy=False))
    values = np.asarray(values, dtype=np.int64)
    for axis, exp in enumerate(expected):
        if not np.array_equal(values[:, axis], exp):
            raise RuntimeError(
                "template selection cannot safely map logical axes to trace headers; "
                "the selected trace header keys do not match the expected geometry"
            )


def _template_file_headers(
    template_path: Path,
    textual_override: bytes | bytearray | np.ndarray | str | list[str] | None,
) -> tuple[bytes, bytearray, bytes]:
    textual = _template_textual_bytes(template_path, textual_override)
    binary = bytearray(_read_file_range(template_path, TEXTUAL_HEADER_SIZE, BINARY_HEADER_SIZE))
    extended_textual = _read_template_extended_textual(template_path, binary)
    return textual, binary, extended_textual


def _template_textual_bytes(
    template_path: Path,
    textual_override: bytes | bytearray | np.ndarray | str | list[str] | None,
) -> bytes:
    if textual_override is None or (isinstance(textual_override, str) and textual_override == ""):
        return _read_file_range(template_path, 0, TEXTUAL_HEADER_SIZE)
    if isinstance(textual_override, str):
        return textual_override.encode("ascii", errors="replace").ljust(TEXTUAL_HEADER_SIZE, b" ")[:TEXTUAL_HEADER_SIZE]
    if isinstance(textual_override, list):
        lines = [str(line).ljust(80)[:80] for line in textual_override[:40]]
        lines.extend([" " * 80] * (40 - len(lines)))
        return "".join(lines).encode("ascii", errors="replace")
    return _strict_textual_bytes(textual_override)


def _read_template_extended_textual(template_path: Path, binary: bytearray) -> bytes:
    ext_count = max(0, _read_i2(binary, 305))
    if ext_count > 0:
        return _read_file_range(template_path, STANDARD_FIRST_TRACE_OFFSET, ext_count * TEXTUAL_HEADER_SIZE)
    first_trace_offset = _read_u8(binary, 321)
    if first_trace_offset > STANDARD_FIRST_TRACE_OFFSET and (first_trace_offset - STANDARD_FIRST_TRACE_OFFSET) % TEXTUAL_HEADER_SIZE == 0:
        return _read_file_range(template_path, STANDARD_FIRST_TRACE_OFFSET, first_trace_offset - STANDARD_FIRST_TRACE_OFFSET)
    return b""


def _read_file_range(path: Path, offset: int, size: int) -> bytes:
    with path.open("rb") as fp:
        fp.seek(int(offset))
        data = fp.read(int(size))
    if len(data) != size:
        raise ValueError(f"failed to read {size} bytes from {path} at offset {offset}")
    return data


def _template_output_sample_interval(segy: Any, sample_axis: _AxisSelection, config: _WriterConfig) -> int:
    if config.dt_new_value > 0:
        return _sample_interval_value(config.dt_new_value)
    source_dt = _source_sample_interval_us(segy)
    if sample_axis.explicit:
        source_dt *= _axis_positive_constant_step(sample_axis.indices, "sample")
    return _sample_interval_value(source_dt)


def _source_sample_interval_us(segy: Any) -> int:
    dt = int(segy.bkeyi2(17))
    if dt <= 0 and int(segy.ntrace) > 0:
        dt = int(segy.keyi2(0, 117))
    return _sample_interval_value(dt)


def _template_output_start_time(segy: Any, sample_axis: _AxisSelection) -> int:
    meta = segy.get_metainfo()
    start_time = int(meta.get("start_time", 0))
    if sample_axis.explicit and sample_axis.size:
        start_time += int(sample_axis.indices[0]) * _source_sample_interval_us(segy) // 1000
    return start_time


def _template_trace_headers(
    segy: Any,
    trace_indices: np.ndarray,
    *,
    sample_count: int,
    sample_interval_us: int,
    start_time_ms: int,
) -> np.ndarray:
    headers = np.empty((int(trace_indices.size), TRACE_HEADER_SIZE), dtype=np.uint8)
    for row, trace_index in enumerate(trace_indices):
        headers[row] = segy.get_trace_header(int(trace_index))
    _set_i2_columns(headers, 115, sample_count)
    _set_i2_columns(headers, 117, sample_interval_us)
    _set_i2_columns(headers, 105, start_time_ms)
    _set_i2_columns(headers, 109, start_time_ms)
    return np.ascontiguousarray(headers, dtype=np.uint8)


def _set_i2_columns(headers: np.ndarray, loc: int, value: int) -> None:
    start = int(loc) - 1
    raw = int(value).to_bytes(2, "big", signed=True)
    headers[:, start : start + 2] = np.frombuffer(raw, dtype=np.uint8)


def _sample_interval_value(value: int) -> int:
    out = int(value)
    if out <= 0:
        raise ValueError("sample_interval_us must be positive")
    if out > 32767:
        raise ValueError("sample_interval_us must fit in the SEG-Y 2-byte header field")
    return out


def _dt_new_value(value: int) -> int:
    out = int(value)
    if out < 0:
        raise ValueError("dt_new must be non-negative")
    if out > 32767:
        raise ValueError("dt_new must fit in the SEG-Y 2-byte header field")
    return out


def _resolve_dt_override(sample_interval_us: int | None, dt_new: int) -> int:
    dt_value = _dt_new_value(dt_new)
    if sample_interval_us is None:
        return dt_value
    sample_interval = _sample_interval_value(sample_interval_us)
    if dt_value > 0 and dt_value != sample_interval:
        raise ValueError("sample_interval_us and dt_new cannot specify different values")
    return sample_interval


def _bytes_from_optional(value: bytes | bytearray | np.ndarray | None) -> bytes:
    if value is None:
        return b""
    if isinstance(value, bytes):
        return value
    if isinstance(value, bytearray):
        return bytes(value)
    array = np.asarray(value, dtype=np.uint8)
    return np.ascontiguousarray(array).tobytes()


def _uint8_array(value: bytes | bytearray | np.ndarray) -> np.ndarray:
    if isinstance(value, np.ndarray):
        return np.ascontiguousarray(value, dtype=np.uint8)
    return np.frombuffer(bytes(value), dtype=np.uint8)


def _created_headers(config: _WriterConfig) -> tuple[bytes, bytearray]:
    if config.shape is None:
        raise ValueError("create mode requires shape")
    if len(config.shape) < 2 or len(config.shape) > 4:
        raise ValueError(f"create mode supports 2D, 3D, or 4D shapes, got {config.shape}")
    if config.sample_format not in SUPPORTED_SAMPLE_BYTES:
        raise ValueError(f"unsupported sample_format: {config.sample_format}")

    textual = _created_textual(config)
    binary = bytearray(BINARY_HEADER_SIZE)
    _set_i4(binary, 5, 1000)
    _set_i2(binary, 13, 1)
    _set_i2(binary, 17, config.sample_interval_us)
    _set_i2(binary, 19, config.sample_interval_us)
    _set_i2(binary, 21, config.shape[-1])
    _set_i2(binary, 23, config.shape[-1])
    _set_i2(binary, 25, int(config.sample_format))
    _set_i2(binary, 29, 4 if len(config.shape) == 3 else 2 if len(config.shape) == 4 else 1)
    _set_i2(binary, 55, 1)
    _set_i2(binary, 303, 1)
    _set_i2(binary, 305, len(_bytes_from_optional(config.extended_textual)) // TEXTUAL_HEADER_SIZE)
    _set_u8(binary, 321, STANDARD_FIRST_TRACE_OFFSET + len(_bytes_from_optional(config.extended_textual)))
    return textual, binary


def _created_textual(config: _WriterConfig) -> bytes:
    if isinstance(config.textual, (bytes, bytearray, np.ndarray)):
        return _strict_textual_bytes(config.textual)
    shape = config.shape
    if shape is None:
        raise ValueError("create mode requires shape")
    if config.start is not None:
        start = list(config.start)
    elif len(shape) == 2:
        start = config.start_time_ms
    elif len(shape) == 3:
        start = [config.iline_start, config.xline_start, config.start_time_ms]
    else:
        start = [config.iline_start, config.xline_start, config.offset_start, config.start_time_ms]
    if len(shape) == 2:
        interval = config.sample_interval_us
    else:
        interval = [config.x_step, config.y_step, config.sample_interval_us]
    meta = assemble_metainfo(shape, dformat=int(config.sample_format), start=start, interval=interval)
    meta["istep"] = int(config.iline_step)
    meta["xstep"] = int(config.xline_step)
    meta["ostep"] = int(config.offset_step)
    return bytes(generate_textual(meta, config.textual))


def _generate_trace_headers(
    config: _WriterConfig,
    block_shape: tuple[int, ...],
    full_shape: tuple[int, ...],
    start: tuple[int, ...],
    trace_offset: int,
) -> np.ndarray:
    ntrace = int(np.prod(block_shape[:-1]))
    headers = np.zeros((ntrace, TRACE_HEADER_SIZE), dtype=np.uint8)
    if len(block_shape) == 2:
        indices = np.arange(start[0], start[0] + block_shape[0], dtype=np.int64)
        for row, trace_index in enumerate(indices):
            _fill_trace_header(config, headers[row], trace_offset + row, (int(trace_index),), full_shape)
        return headers

    if len(block_shape) == 3:
        row = 0
        for ii in range(start[0], start[0] + block_shape[0]):
            for ix in range(start[1], start[1] + block_shape[1]):
                _fill_trace_header(config, headers[row], trace_offset + row, (ii, ix), full_shape)
                row += 1
        return headers

    row = 0
    for ii in range(start[0], start[0] + block_shape[0]):
        for ix in range(start[1], start[1] + block_shape[1]):
            for io in range(start[2], start[2] + block_shape[2]):
                _fill_trace_header(config, headers[row], trace_offset + row, (ii, ix, io), full_shape)
                row += 1
    return headers


def _fill_trace_header(
    config: _WriterConfig,
    header: np.ndarray,
    trace_index: int,
    coord: tuple[int, ...],
    full_shape: tuple[int, ...],
) -> None:
    seq = int(trace_index) + 1
    _set_i4(header, 1, seq)
    _set_i4(header, 5, seq)
    _set_i2(header, 29, 1)
    _set_i2(header, 71, 1)
    _set_i2(header, 105, config.start_time_ms)
    _set_i2(header, 109, config.start_time_ms)
    _set_i2(header, 115, config.sample_count or (config.shape[-1] if config.shape else 0))
    _set_i2(header, 117, config.sample_interval_us)

    iline_index = coord[0]
    xline_index = coord[1] if len(coord) >= 2 else 0
    offset_index = coord[2] if len(coord) >= 3 else None
    iline = int(config.iline_start + iline_index * config.iline_step)
    xline = int(config.xline_start + xline_index * config.xline_step)
    x_coord = int(config.x_start + iline_index * config.x_step)
    y_coord = int(config.y_start + xline_index * config.y_step)
    offset = None if offset_index is None else int(config.offset_start + offset_index * config.offset_step)
    iline = _geometry_value(config, "iline", full_shape, coord, iline)
    xline = _geometry_value(config, "xline", full_shape, coord, xline)
    x_coord = _geometry_value(config, "x", full_shape, coord, x_coord)
    y_coord = _geometry_value(config, "y", full_shape, coord, y_coord)
    if offset_index is not None:
        offset = _geometry_value(config, "offset", full_shape, coord, offset)

    _set_i4(header, 9, iline)
    _set_i4(header, 17, xline)
    _set_i4(header, 21, xline)
    _set_i4(header, 73, x_coord)
    _set_i4(header, 77, y_coord)
    _set_i4(header, 181, x_coord)
    _set_i4(header, 185, y_coord)
    _set_i4(header, 189, iline)
    _set_i4(header, 193, xline)
    if offset is not None:
        _set_i4(header, 37, offset)


def _geometry_value(
    config: _WriterConfig,
    name: str,
    full_shape: tuple[int, ...],
    coord: tuple[int, ...],
    fallback: int,
) -> int:
    geom = config.geometry or {}
    if name not in geom:
        return int(fallback)
    value = geom[name]
    arr = np.asarray(value)
    if arr.ndim == 0:
        return int(arr)
    trace_shape = tuple(int(v) for v in full_shape[:-1])
    if tuple(arr.shape) == trace_shape:
        return int(arr[coord])
    if arr.ndim == 1:
        index = int(np.ravel_multi_index(coord, trace_shape))
        if index >= arr.size:
            raise ValueError(f"geometry[{name!r}] has too few values for coordinate {coord}")
        return int(arr[index])
    raise ValueError(
        f"geometry[{name!r}] must be scalar, 1D file-order values, or have shape {trace_shape}; got {arr.shape}"
    )


def _optional_int_list(value: Iterable[int] | None) -> list[int] | None:
    if value is None:
        return None
    return [int(v) for v in value]


def _strict_binary_bytes(value: bytes | bytearray | np.ndarray | None) -> bytes:
    if value is None:
        raise ValueError("from_headers mode requires binary header bytes")
    out = _bytes_from_optional(value)
    if len(out) != BINARY_HEADER_SIZE:
        raise ValueError(f"binary header must be exactly {BINARY_HEADER_SIZE} bytes, got {len(out)}")
    return out


def _strict_textual_bytes(value: bytes | bytearray | np.ndarray | str | list[str] | None) -> bytes:
    if value is None:
        raise ValueError("from_headers mode requires textual header bytes")
    if isinstance(value, str):
        out = value.encode("ascii")
    elif isinstance(value, list):
        text = "".join(str(line).ljust(80)[:80] for line in value)
        out = text.encode("ascii")
    else:
        out = _bytes_from_optional(value)
    if len(out) != TEXTUAL_HEADER_SIZE:
        raise ValueError(f"textual header must be exactly {TEXTUAL_HEADER_SIZE} bytes, got {len(out)}")
    return out


def _trace_headers_array(value: Any) -> np.ndarray:
    out = np.asarray(value, dtype=np.uint8)
    if out.ndim == 1:
        if out.size != TRACE_HEADER_SIZE:
            raise ValueError(f"single trace header must have {TRACE_HEADER_SIZE} bytes")
        out = out.reshape(1, TRACE_HEADER_SIZE)
    if out.ndim != 2 or out.shape[1] != TRACE_HEADER_SIZE:
        raise ValueError(f"trace_headers must have shape (ntrace, {TRACE_HEADER_SIZE})")
    return np.ascontiguousarray(out, dtype=np.uint8)


def _trace_headers_with_sample_interval(headers: np.ndarray, sample_interval_us: int | None) -> np.ndarray:
    if sample_interval_us is None:
        return headers
    out = np.array(headers, dtype=np.uint8, copy=True, order="C")
    raw = int(sample_interval_us).to_bytes(2, "big", signed=True)
    out[:, 116:118] = np.frombuffer(raw, dtype=np.uint8)
    return out


def _validate_extended_textual(value: bytes) -> None:
    if value and len(value) % TEXTUAL_HEADER_SIZE != 0:
        raise ValueError("extended textual headers must be a multiple of 3200 bytes")


def _read_i2(header: bytes | bytearray | np.ndarray | None, loc: int) -> int:
    if header is None:
        raise ValueError("header is not available")
    start = int(loc) - 1
    return int(np.frombuffer(bytes(header[start : start + 2]), dtype=">i2")[0])


def _read_u8(header: bytes | bytearray | np.ndarray | None, loc: int) -> int:
    if header is None:
        raise ValueError("header is not available")
    start = int(loc) - 1
    return int(np.frombuffer(bytes(header[start : start + 8]), dtype=">u8")[0])


def _set_i2(header: bytearray | np.ndarray, loc: int, value: int) -> None:
    start = int(loc) - 1
    raw = int(value).to_bytes(2, "big", signed=True)
    _set_raw(header, start, raw)


def _set_i4(header: bytearray | np.ndarray, loc: int, value: int) -> None:
    start = int(loc) - 1
    raw = int(value).to_bytes(4, "big", signed=True)
    _set_raw(header, start, raw)


def _set_u8(header: bytearray | np.ndarray, loc: int, value: int) -> None:
    start = int(loc) - 1
    raw = int(value).to_bytes(8, "big", signed=False)
    _set_raw(header, start, raw)


def _set_raw(header: bytearray | np.ndarray, start: int, raw: bytes) -> None:
    if isinstance(header, np.ndarray):
        header[start : start + len(raw)] = np.frombuffer(raw, dtype=np.uint8)
    else:
        header[start : start + len(raw)] = raw
