# cigsegy

<table>
  <tr>
    <td><a href="./README.rst">English</a></td>
    <td><b>中文</b></td>
  </tr>
</table>

## 中文版的README不是最新的, 只支持到1.1.5, 详情请参考英文版

**此文件最后更新版本为 1.1.5**

一个读写 `segy` 格式地震数据的 `python` 和 `c++` 工具。可以将 `segy` 格式文件读到内存或者直接转为二进制文件，也可以将一个 `numpy` 数据存储为`segy`格式的文件。

**特点**:

- 快，底层使用c++实现
- 可以在`python`中使用，使用了pybind11将c++包装为python可调用的库
- 可以处理**规则**和**不规则**的地震数据，比如工区不是一个矩形（有缺失）或数据间隔不为1
- 可以使用**已有的segy文件的道头**进行创建新的segy文件

你可以通过pip简单的安装 cigsegy
```bash
pip install cigsegy
```
你可以读取各种segy文件（无论是否规则）通过一个简单的函数:
```python
d = cigsegy.fromfile('f.segy', iline=, xline=, istep=, xstep=)
```
你可以在创建一个新的segy文件的时候，使用一个老的segy文件的道头，相当于对老的segy文件的数据进行替换:
```python
cigsegy.create_by_sharing_header('out.segy', 'header.segy', d, iline=, xline=, istep=, xstep=)
```

更多的用例见：[Cases](./Cases.md)

### v1.1.5 更新内容
- 在 jupyter notebook 中禁止显示进度条
- 允许 `numpy slice`， 即一个子数组传入 `create_by_sharing_header` 函数中，现在不需要在传入之前使用 `data.clone()`来使内存连续了。
- 添加了一个 `tools.read_header`来获取 binary/trace 的全部道头信息。
- 添加了一个 `tools.get_metaInfo()` 来获取元信息，以字典的形式存在

### 目录

- [安装](#Installation)
- [用法](#usage)
- [用到的第三方依赖](#ThirdPart)
- [两个在终端运行的可执行文件](#Executables)
- [局限性](#Limitations)
- [感谢](#Acknowledge)

<p id="Installation"></p>

### 安装

- 使用 pip 安装:


```
pip install cigsegy
```

- 本地安装
首先你需要安装依赖库 `fmt` 和 `pybind11`.
```bash
# linux
sudo apt-get install python3-pybind11 libfmt-dev

# mac
brew install pybind11 fmt
```
你也可以手动安装依赖哭 `fmt` 和 `pybind11` .
```bash
# 安装 fmt
mkdir thridPart && cd thridPart/
git clone https://github.com/fmtlib/fmt.git
# 现在 fmt 被安装到 /xxx/cigsegy/thridPart/fmt


# 使用 pypi 安装 pybind11
pip install pybind11
# 现在 pybind11 被安装到 /xxx/lib/python3.8/site-packages/pybind11/
```

如果你仅需要python的版本，你可以运行:
```bash
pip install -U pip
pip install .

# 如果 fmt 没有被包含在系统的路径里
pip install . --install-option="--fmt_root=/xxx/fmt"

# 如果你需要的是一个 wheel 文件
# pip wheel . --build-option="--fmt_root=/xxx/fmt"
```


如果你还需要使用C++编写，可以在终端运行的两个可执行文件:
```bash
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=your-install-dir
make
make install

# 如果 fmt 和 pybind11 没有在系统的搜索路径里
cmake .. -DCMAKE_INSTALL_PREFIX=your-install-dir \
  -Dpybind11_ROOT=/xxx/lib/python3.8/site-packages/pybind11/ \
  -Dfmt_ROOT=/xxx/cigsegy/thridPart/fmt/build/fmt/

# 如果你想要其他python版本的结果, 
# 添加 -DPYTHON_EXECUTABL 和 -DPYTHON_LIBRARIES 选项指定python的信息
cmake .. -DCMAKE_INSTALL_PREFIX=your-install-dir -DPYTHON_EXECUTABLE=/xxx/bin/python -DPYTHON_LIBRARIES=/xxx/lib/
```


<p id="usage"></p>

### 用法
> 更多例子: [Cases](Cases.md)

#### 读


```python
import cigsegy

# 显示segy文件的前3200个字节的文本道头信息
cigsegy.textual_header('fx.segy')

# 显示一个segy文件的一些关键信息，比如shape, 数据格式等
cigsegy.metaInfo('fx.segy', iline=189, xline=193, istep=1, xstep=1)

# 如果不知道iline等信息，可以使用
cigsegy.metaInfo('fx.segy', use_guess=True)

# 将fx.segy文件读入到内存中，d 为 numpy.ndarray 格式
d = cigsegy.fromfile('fx.segy')

# 读取fx.segy文件，并指定inline/crossline的存储位置和间隔
d = cigsegy.fromfile('fx.segy', iline=189, xline=193, istep=1, xstep=1)

# 不读入到内存，而是直接将其转存到二进制文件 fx.dat 中
# 在 fx.segy 非常大，超过内存大小的时候很有用
cigsegy.tofile('fx.segy', 'fx.dat', iline=189, xline=193, istep=1, xstep=1)

# 如果道头信息损坏，或者被人为隐藏了, 可以知道数据尺寸和格式，忽略掉道头信息
# 假设 shape: (inline, xline, dt) = (589, 762, 1001)
d = cigsegy.fromfile_ignore_header('fx.segy', 589, 762, 1001, format=5)
# 注意 format 5 意味着每个数据点为 4 字节的 IEEE 浮点数
# format 1 意味着每个数据点为 4 字节的 IBM 浮点数.

# 存储为一个二进制文件, i.e., 不读入到内存中
cigsegy.tofile_ignore_header('fx.segy', 'fx.dat', 589, 762, 1001, format=5)

# 仅仅获得数据作为一个2维数组
data = cigsegy.collect('fx.segy')
# data.shape = (ntraces, nt)

# 读取一个 2D 的 segy 数据
data = cigsegy.collect('fx2d.segy')

# 如果不知道 inline/crossline 的存储位置和间隔, 
# 可以使用 `guess` 函数去猜测
locs = cigsegy.tools.guess('fx.segy')

# 或者直接使用 `fromfile_by_guess`, `tofile_by_guess` 这两个函数
d = cigsegy.tools.fromfile_by_guess('fx.segy')
cigsegy.tools.tofile_by_guess('fx.segy', 'fx.dat')

# 通过 inline and crossline 画出工区的形状
cigsegy.tools.plot_region('fx.segy')

# 获取全部的道头信息, in dict
cigsegy.tools.read_header('fx.segy', type='bh') # binary header, print it
header = cigsegy.tools.read_header('fx.segy', type='bh', printstr=False) # in dict
header = cigsegy.tools.read_header('fx.segy', type='th', n=10, printstr=False) # 10-th trace header

# 获取元信息 in dict
meta = cigsegy.tools.get_metaInfo('fx.segy', iline=, xline=, ...)
```


#### 写

```python
import numpy as np
import cigsegy

# fx.segy

d = ... # np.ndarray, shape = (589, 762, 1001)

# 创建 "new_fx.segy", 使用数据 'd' 和来自于 fx.segy 的道头
cigsegy.create_by_sharing_header('new_fx.segy', 'fx.segy', d, 
            iline=189, xline=193, istep=1, xstep=1)

# new_fx.dat, shape = (589, 762, 1001), 
# 当 new_fx.dat 非常大时，直接内存映射
shape = (589, 762, 1001)
cigsegy.create_by_sharing_header('new_fx.segy', 'fx.segy', 'new_fx.dat', 
            shape, iline=189, xline=193, istep=1, xstep=1)

# sub_fx.dat, shape = (401, 535, 800), 相对于f3.segy整个数据
# 的offset是 (65, 73, 100), 你可以使用 `offset` 参数去创建一个小的segy文件，
# 但是使用大的（f3.segy）segy文件的道头
d = ...
shape = (401, 535, 800) # d.shape
cigsegy.create_by_sharing_header('new_fx.segy', 'fx.segy', d, 
            iline=189, xline=193, istep=1, xstep=1, offset=(65, 73, 100))
cigsegy.create_by_sharing_header('sub_fx.segy', 'fx.segy', 'sub_fx.dat', 
            shape, iline=189, xline=193, istep=1, xstep=1, offset=(65, 73, 100))

# 如果不知道 inline&crossline 的位置和间隔，可以使用
cigsegy.tools.create_by_sharing_header_guess('new_fx.segy', 'fx.segy', d)
cigsegy.tools.create_by_sharing_header_guess('new_fx.segy', 'fx.segy', 
            'new_fx.dat', shape)

# 简单创建, 不使用已有道头信息
cigsegy.create('new2_fx.segy', 'new_fx.dat', 589, 762, 1001, ...)
cigsegy.create('new2_fx.segy', d, 589, 762, 1001, ...)
```


所有的用于创建一个segy文件的函数都添加了`custom_info`参数用于自定义3200字节的文本道头的前12行信息,
比如：
```python
my_info = ['This is a test examples', 
           'For detail information, you can see: ', 
           'https://github.com/JintaoLee-Roger/cigsegy',
           'Thank you for using cigsegy']
cigsegy.create_by_sharing_header('new_fx.segy', 'fx.segy', d, 
            iline=189, xline=193, istep=1, xstep=1, offset=(65, 73, 100),
            custom_info=my_info)
cigsegy.textual_header('new_fx.segy')
# C01 This is a test examples
# C02 For detail information, you can see: 
# C03 https://github.com/JintaoLee-Roger/cigsegy
# C04 Thank you for using cigsegy
# C05
# ...
```


#### class `Pysegy`
```python
import cigsegy

# read
segy = cigsegy.Pysegy('fx.segy')
print(segy.textual_header())
segy.setInlineLocation(5)
segy.setCrosslineLocation(5)
d = segy.read()
d2 = segy.read(1, 5, 2, 10, 3, 9)
...
# 更多细节可以参考 cigsegy.pyi

# create
segy = cigsegy.Pysegy('fx.dat', 589, 762, 1001)
segy.set...
segy.create('new_fx.segy', d)
# m更多细节可以参考 cigsegy.pyi
```


<p id="ThirdPart"></p>

### 用到的第三方依赖

- `src/include/mio.hpp` is from [mandreyel/mio](https://github.com/mandreyel/mio): Cross-platform C++11 header-only library for memory mapped file IO

- `src/include/progressbar.hpp` is from [gipert/progressbar](https://github.com/gipert/progressbar): A very simple progress bar for C++ loops

- `tools/cxxopts.hpp` is from [jarro2783/cxxopts](https://github.com/jarro2783/cxxopts): Lightweight C++ command line option parser

- [pybind/pybind11](https://github.com/pybind/pybind11): Seamless operability between C++11 and Python

- [fmtlib/fmt](https://github.com/fmtlib/fmt): A modern formatting library


<p id="Executables"></p>


### 两个在终端运行的可执行文件
两个使用c++编写的位于 `your-install-dir/bin/` 的可执行文件: `SEGYRead` and `SEGYCreate`.

`SEGYRead`
```bash
./SEGYRead
# ./SEGYRead - a tool for segy file reading to binary file
# Usage:
#   ./SEGYRead [OPTION...] positional parameters

#   -o, --out arg               out binary file name
#   -f, --fills arg             the number to fill the miss trace, can be any 
#                               float or nan, or NAN
#   -z, --inline-loc arg        inline field in trace header, default is 189
#   -c, --crossline-loc arg     crossline field in trace header, default is 
#                               193
#       --istep arg             inline step
#       --xstep arg             crossline step
#       --xloc arg              X field in trace header, default is 73
#       --yloc arg              Y field in trace header, default is 77
#   -p, --print_textual_header  print 3200 bytes textual header
#   -m, --meta_info             print meta info
#       --ignore-header         reading segy by ignoring header and specify 
#                               shape
#   -d, --dimensions arg        the dimensions (x, y, z) or (nt, ncrossline, 
#                               ninline), use as '-d 128,128,256' (Required)

# Examples:
#   ./SEGYRead -p f3.segy             : show textual header
#   ./SEGYRead -m f3.segy             : show meta information
#   ./SEGYRead -o f3.dat f3.segy      : convert
#   ./SEGYRead -i f3.segy -o f3.dat   : convert
#   ./SEGYRead -o f3.dat -z 5 f3.segy : convert by specify inline field
#   ./SEGYRead -o f3.dat -z 5 --istep 2 f3.segy : convert by specify inline field and step
#   ./SEGYRead -o f3.dat -f nan f3.segy : convert and fill with nan
#   ./SEGYRead -o f3.dat --ignore-header -d 236,789,890 f3.segy : ignore header and specify shape
```

`SEGYCreate`
```bash
./SEGYCreate
# ./SEGYCreate - a tool for creating a segy file from a binary file
# Usage:
#   ./SEGYCreate [OPTION...] positional parameters

#   -o, --out arg            out segy file name (Required)
#   -d, --dimensions arg     the dimensions (x, y, z) or (nt, ncrossline, 
#                            ninline), use as '-d 128,128,256' (Required)
#   -z, --inline-loc arg     set inline field in trace header, default is 189
#   -c, --crossline-loc arg  set crossline field in trace header, default is 
#                            193
#       --dt arg             set sample interval, default is 4000
#   -f, --data_format arg    data format code, 4 bytes (1 for IBM, 5 for 
#                            IEEE) floating point, defualt is 5
#       --dx arg             set X interval
#       --dy arg             set Y interval
#       --min-inline arg     set start inline number
#       --min-crossline arg  set start crossline number
#       --start-time arg     set start time for each trace
#   -s, --share              create a segy file using a exist segy file 
#                            header
#       --header             the shared header segy file
#       --istep arg          inline step for the header segy file
#       --xstep arg          crossline step for the header segy file

# Examples:
#   ./SEGYCreate -i test.dat -o test.segy -d 128,128,256 : convert
#   ./SEGYCreate -o test.segy -d 128,128,256 test.dat : convert
#   ./SEGYCreate -o test.segy -d 128,128,256 -f 5 test.dat : specify data format
#   ./SEGYCreate -o test.segy -d 128,128,256 --dt 2000 test.dat : specify time interval
#   ./SEGYCreate -s -i test.dat -o test.segy --header header.segy -d 128,128,256 -z 189 -c 193 : using sharing header segy file
```

<p id="Limitations"></p>

### 局限性

- 只测试过叠后数据

<p id="Acknowledge"></p>

### 致谢

感谢 [jiarunyangustc](https://github.com/jiarunyangustc) 和 [shenghanlin](https://github.com/shenghanlin) 提供的bug和需求，并且帮助进行了多个数据的测试。
