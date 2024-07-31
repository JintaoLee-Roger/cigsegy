import os, sys
import subprocess
from pathlib import Path


def install_package(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])


try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext
except ImportError:
    install_package("pybind11")
    from pybind11.setup_helpers import Pybind11Extension, build_ext

from setuptools import setup
from distutils.ccompiler import get_default_compiler

fmt_root = ''
for v in sys.argv:
    if v.startswith('--fmt_root='):
        sys.argv.remove(v)
        # fmt_root = v.removeprefix('--fmt_root=')
        fmt_root = v[11:]
        break

cwd = Path(__file__).resolve().parent

package_name = "cigsegy"
git_hash = "unknown"

version_path = Path(__file__).parent / "VERSION.txt"
if not version_path.exists():
    raise FileNotFoundError("VERSION.txt file not found")
version = version_path.read_text().strip()

if not version:
    raise RuntimeError("Failed to parse version from VERSION")


def get_extensions():
    # add segy
    ext_modules = []
    sources = ['src/segy.cpp', 'python/PySegy.cpp']
    include_dirs = ['src/include']
    if fmt_root:
        include_dirs.append(str(Path(fmt_root) / 'include'))
    if get_default_compiler() == 'msvc':
        extra_compile_args = ['/wd4244', '/wd4996', '/wd4819']
    else:
        extra_compile_args = ["-std=c++11", "-Wall", "-O3"]

    extra_compile_args.append("-DUSE_PYBIND11")

    extra_link_args = []
    # extra_link_args = ['-lfmt']

    ext_modules.append(
        Pybind11Extension(name=f'{package_name}.{package_name}',
                          sources=[str(s) for s in sources],
                          include_dirs=include_dirs,
                          extra_compile_args=extra_compile_args,
                          extra_link_args=extra_link_args))

    return ext_modules


setup(name=package_name,
      version=version,
      long_description=open('README.rst').read(),
      long_description_content_type='text/x-rst',
      author='Jintao Li',
      url='https://github.com/JintaoLee-Roger/cigsegy',
      license='MIT',
      install_requires=['numpy'],
      python_requires=">=3.6",
      ext_modules=get_extensions(),
      cmdclass={"build_ext": build_ext},
      packages=['cigsegy'],
      package_dir={'cigsegy': 'python'},
      include_package_data=True,
      exclude_package_data={'cigsegy': ['*.cpp', '*.txt', 'setup.py']})
