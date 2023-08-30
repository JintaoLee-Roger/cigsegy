import os, sys
import subprocess

from pathlib import Path

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
version = "1.1.6"
git_hash = "unknown"

try:
    git_hash = (subprocess.check_output(["git", "rev-parse", "HEAD"],
                                        cwd=cwd).decode().strip())
except (FileNotFoundError, subprocess.CalledProcessError):
    pass


def write_version_file():
    path = 'python/version.py'
    with open(path, "w") as f:
        f.write(f'__version__ = "{version}"\n')
        f.write(f'git_version = "{git_hash}"\n')


write_version_file()


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
    extra_link_args = []
    # extra_link_args = ['-lfmt']

    ext_modules.append(
        Pybind11Extension(name=f'{package_name}.{package_name}',
                          sources=[str(s) for s in sources],
                          include_dirs=include_dirs,
                          extra_compile_args=extra_compile_args,
                          extra_link_args=extra_link_args))

    return ext_modules


setup(
    name=package_name,
    version=version,
    description=
    'A tool for segy-format file reading and segy-format creating from binary file',
    author='roger',
    url='https://github.com/JintaoLee-Roger/cigsegy',
    license='MIT',
    # install_requires=['numpy'],
    python_requires=">=3.6",
    ext_modules=get_extensions(),
    cmdclass={"build_ext": build_ext},
    packages=['cigsegy'],
    package_dir={'cigsegy': 'python'},
    include_package_data=True,
    exclude_package_data={'cigsegy': ['*.cpp', '*.txt', 'setup.py']})
