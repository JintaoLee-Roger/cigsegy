# Copyright (c) 2024 Jintao Li.
# Computational and Interpretation Group (CIG),
# University of Science and Technology of China (USTC).
# All rights reserved.

import os, sys
import subprocess
from pathlib import Path
from setuptools import setup, find_packages
from distutils.ccompiler import get_default_compiler


def install_package(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])


try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext
except ImportError:
    install_package("pybind11")
    from pybind11.setup_helpers import Pybind11Extension, build_ext

cwd = Path(__file__).resolve().parent

package_name = "cigse"

version = '1.2.0'

if not version:
    raise RuntimeError("Failed to parse version from VERSION")

if get_default_compiler() == 'msvc':
    extra_compile_args = ['/wd4244', '/wd4996', '/wd4819']
else:
    extra_compile_args = ["-std=c++11", "-Wall", "-O3"]

extra_compile_args += ["-DUSE_PYBIND11", "-undefined dynamic_lookup"]

ext_modules = [
    Pybind11Extension(
        "cigse.cpp._CXX_SEGY",
        [
            "cigse/cpp/segywrap.cpp",
            "cigse/cpp/segyrw.cpp",
        ],
        extra_compile_args=extra_compile_args,
    ),
]

setup(name=package_name,
      version=version,
      long_description=open('README.rst').read(),
      long_description_content_type='text/x-rst',
      author='Jintao Li',
      url='https://github.com/JintaoLee-Roger/cigse',
      license='MIT',
      install_requires=['numpy'],
      python_requires=">=3.6",
      setup_requires=['pybind11'],
      ext_modules=ext_modules,
      cmdclass={"build_ext": build_ext},
      packages=find_packages(exclude=['docs', 'python', 'tools']),
      include_package_data=True,
      exclude_package_data={'cigse': ['*.cpp', '*.txt', 'setup.py']})
