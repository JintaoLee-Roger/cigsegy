# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

from pathlib import Path
from setuptools import setup, find_packages
from distutils.ccompiler import get_default_compiler

try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext
except ImportError:
    raise RuntimeError(
        "pybind11 is required to build cigsegy. Install build dependencies first."
    )

cwd = Path(__file__).resolve().parent

package_name = "cigsegy"

version_file = cwd / 'VERSION.txt'
with open(version_file) as vf:
    version = vf.read().strip()

if not version:
    raise RuntimeError("Failed to parse version from VERSION")

if get_default_compiler() == 'msvc':
    extra_compile_args = ['/wd4244', '/wd4996', '/wd4819']
else:
    extra_compile_args = ["-std=c++14", "-O3", "-undefined dynamic_lookup"]

extra_compile_args += ["-DUSE_PYBIND11"]

ext_modules = [
    Pybind11Extension(
        "cigsegy.cpp._CXX_SEGY",
        [
            "cigsegy/cpp/segywrap.cpp",
            "cigsegy/cpp/segyrw.cpp",
            "cigsegy/cpp/segywriter.cpp",
        ],
        extra_compile_args=extra_compile_args,
    ),
]

setup(
    name=package_name,
    version=version,
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    author='Jintao Li',
    url='https://github.com/JintaoLee-Roger/cigsegy',
    license='MIT',
    install_requires=['numpy', 'tqdm'],
    python_requires=">=3.8",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: C++",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "License :: OSI Approved :: MIT License",
    ],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    entry_points={
        "console_scripts": [
            "cigsegy=cigsegy.cli:main",
        ],
    },
    packages=find_packages(exclude=['docs', 'python', 'tools', 'tests']),
    include_package_data=True,
    package_data={
        "cigsegy": ["*.pyi"],
        "cigsegy.cpp": ["*.pyi"],
    },
)
