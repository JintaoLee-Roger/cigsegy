# Copyright (c) 2023 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages


package_name = "segy"

description = "segy"

ext_modules = [
    Pybind11Extension(
        "segy._CXX_SEGYR",  # 模块名称
        ["segy/pysegyr.cpp"],  # 源文件路径
        # 可以在这里添加其他需要的编译器和链接器标志
        # extra_compile_args=['-O3', '-Wall'],
    ),
]

setup(
    name=package_name,
    author="Jintao Li",
    license='MIT',
    description=description,
    long_description_content_type="text/markdown",
    python_requires=">=3.7",
    packages=find_packages(),
    setup_requires=['pybind11'],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
