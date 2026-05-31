# cigsegy

`cigsegy` 是一个用于读取、检查、转换和写出 SEG-Y 地震数据的 Python/C++ 工具库。

最常见的流程应该很直接：

1. 用 `textual_header` 看文本道头
2. 用 `metaInfo` 扫描 geometry，优先让 `cigsegy` 自动推断字段位置
3. 用 `fromfile` 直接读成 NumPy 数组
4. 在 NumPy 里处理数据
5. 用 `SegyWriter` 写回 SEG-Y

`SegyNP` 仍然很重要，但它更适合大文件、随机切片、交互式可视化、或者不想一次性读入内存的场景。

## 安装

```bash
pip install cigsegy
```

本地开发：

```bash
pip install -e . --config-settings editable_mode=strict
```

## 查看头信息

先看 textual header。很多 SEG-Y 会在里面写明 inline、crossline、坐标等字段的位置。

```python
import cigsegy

cigsegy.textual_header("input.sgy")
```

然后扫 metadata。默认情况下，`cigsegy` 会尽量从 trace header 自动推断 inline、crossline、offset、step 和坐标字段的位置：

```python
cigsegy.metaInfo("input.sgy")
```

只有在自动推断失败、结果不符合预期，或者 SEG-Y 头字段不规范时，才需要显式传入 byte location：

```python
cigsegy.metaInfo("input.sgy", iline=189, xline=193)
```

需要更细的 binary header 或 trace header 字段时：

```python
binary = cigsegy.tools.read_header("input.sgy", type="bh", printstr=False)
trace0 = cigsegy.tools.read_header("input.sgy", type="th", n=0, printstr=False)
```

## 直接读成 NumPy

很多处理脚本最自然的方式就是直接把 SEG-Y 读成 NumPy 数组。

```python
data = cigsegy.fromfile("input.sgy")
print(data.shape)  # (n_inline, n_xline, n_sample)
```

4D/pre-stack 也先尝试自动推断：

```python
gathers = cigsegy.fromfile("gather.sgy")
```

如果推断结果不对，再传入已知字段：

```python
data = cigsegy.fromfile("input.sgy", iline=189, xline=193)
gathers = cigsegy.fromfile("gather.sgy", iline=189, xline=193, offset=37)
```

2D trace collection：

```python
traces = cigsegy.collect("line.sgy")
```

如果数据太大，不想整体读入内存，可以直接流式写到硬盘：

```python
cigsegy.tofile("input.sgy", "samples.dat")  # raw float32，没有 shape 信息
cigsegy.to_npy("input.sgy", "samples.npy")  # .npy，可以后续 mmap 读取
```

## 处理后写回 SEG-Y

输出 trace 与原 SEG-Y 的 trace 一一对应时，使用 `SegyWriter.from_template`。
这就是新代码里替代 `create_by_sharing_header` 的主路径。

```python
processed = process(data)

b = cigsegy.SegyWriter.from_template("input.sgy", "processed.sgy")
b.overwrite(True)

with b.open() as w:
    w.write(processed)
```

连续子体：

```python
sub = data[100:300, 40:200, 0:800]

b = cigsegy.SegyWriter.from_template("input.sgy", "sub.sgy")
b.overwrite(True)

with b.open() as w:
    w.write(sub, start=(100, 40, 0))
```

时间超分，比如 `2ms -> 1ms`：

```python
b = cigsegy.SegyWriter.from_template("input_2ms.sgy", "output_1ms.sgy")
b.strict(False).sample_interval_us(1000).overwrite(True)

with b.open() as w:
    w.write(super_res_data)
```

空间抽稀，header 会从原始文件中真实存在的 trace 复制：

```python
thin = data[::2, ::3, :]

b = cigsegy.SegyWriter.from_template("input.sgy", "thin.sgy")
b.select(iline=slice(None, None, 2), xline=slice(None, None, 3))
b.overwrite(True)

with b.open() as w:
    w.write(thin)
```

空间抽稀 + 时间降采样：

```python
down = data[::2, ::3, ::4]

b = cigsegy.SegyWriter.from_template("input_1ms.sgy", "thin_4ms.sgy")
b.select(
    iline=slice(None, None, 2),
    xline=slice(None, None, 3),
    sample=slice(None, None, 4),
)
b.overwrite(True)

with b.open() as w:
    w.write(down)
```

## 使用 SegyNP 按需读取

`SegyNP` 不会一次性把整个文件读入内存，只有索引时才读取对应的数据。

```python
vol = cigsegy.SegyNP("input.sgy")

iline = vol[100, :, :]
xline = vol[:, 200, :]
time_slice = vol[:, :, 300]
cube = vol[100:140, 200:260, 300:700]
```

2D：

```python
traces = cigsegy.SegyNP("line.sgy", ndim=2)
trace = traces[100]
part = traces[1000:1200, :]
```

4D/pre-stack 也可以先让 `cigsegy` 自动推断：

```python
gathers = cigsegy.SegyNP("gather.sgy")
one_gather = gathers[20, 30, :, :]
```

## 使用已有 header block 写出

如果其他容器里已经保存了 textual header、binary header、trace headers，可以直接按 block 写。

```python
b = cigsegy.SegyWriter.from_headers("out.sgy")
b.textual(textual).binary(binary)
b.sample_format(1).sample_count(1001)
b.overwrite(True)

with b.open() as w:
    for trace_headers, samples in blocks:
        w.write_trace_block(trace_headers, samples)
```

如果 sample bytes 已经是 SEG-Y 编码后的原始字节，可以不重新编码：

```python
with b.open() as w:
    for trace_headers, sample_bytes in raw_blocks:
        w.write_raw_trace_block(trace_headers, sample_bytes)
```

## 从零创建规则 SEG-Y

```python
b = cigsegy.SegyWriter.create(
    "created.sgy",
    shape=(589, 762, 1001),
    sample_format=5,
    sample_interval_us=2000,
    overwrite=True,
)
b.grid(iline_start=1000, xline_start=2000)
b.origin(x_start=600000, y_start=4100000, x_step=25, y_step=25)

with b.open() as w:
    w.write(data)
```

4D：

```python
b = cigsegy.SegyWriter.create(
    "gather.sgy",
    shape=(120, 200, 48, 1500),
    sample_format=5,
    sample_interval_us=2000,
    overwrite=True,
)
b.as_4d().grid(iline_start=1000, xline_start=2000, offset_start=100)

with b.open() as w:
    w.write(gather_data)
```

## 旧的短函数

这些函数仍然可用：

- `cigsegy.fromfile`
- `cigsegy.collect`
- `cigsegy.tofile`
- `cigsegy.create_by_sharing_header`
- `cigsegy.create`

新写 SEG-Y 的代码建议优先使用 `SegyWriter`，因为它同时覆盖模板写出、已有 header block 写出、raw sample bytes 写出、以及从零生成规则 header。
