# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import subprocess
import sys
import shutil
from pathlib import Path

DIR = Path(__file__).parent.resolve()

project = 'cigsegy'
copyright = '2023, Roger Lee'
author = 'Roger Lee'
version = 'v1.1.6'
language = 'en'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["breathe", "sphinx.ext.autodoc", "sphinx.ext.napoleon"]

breathe_projects = {"cigsegy": "_doxygenxml/xml/"}
breathe_default_project = "cigsegy"
breathe_domain_by_extension = {"h": "cpp"}

templates_path = ['_templates']
exclude_patterns = [
    '_build', 'Thumbs.db', '.DS_Store', '_doxygenxml', 'cigsegy'
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'
# html_static_path = ['_static']


def generate_doxygen_xml(app):
    build_dir = os.path.join(app.confdir, "_doxygenxml")
    if not os.path.exists(build_dir):
        os.mkdir(build_dir)

    try:
        subprocess.call(["doxygen", "--version"])
        retcode = subprocess.call(["doxygen"], cwd=app.confdir)
        if retcode < 0:
            sys.stderr.write(f"doxygen error code: {-retcode}\n")
    except OSError as e:
        sys.stderr.write(f"doxygen execution failed: {e}\n")


def prepare(app):

    # copy readme
    with open(DIR.parent / "README.rst") as f:
        contents = f.read()
        contents = contents.replace('docs/assets/logo.svg', 'assets/logo.svg')

    with open(DIR / "readme.rst", "w") as f:
        f.write(contents)

    sourcedir = DIR.parent / 'python'
    targetdir = DIR / 'cigsegy'
    if not os.path.exists(targetdir):
        os.mkdir(targetdir)

    # copy 'python' dir as 'cigsegy'
    with open(sourcedir / '__init__.py') as f:
        contents = f.read()
        contents = contents.replace('from .cigsegy', 'from .cigsegyc')

    with open(targetdir / '__init__.py', 'w') as f:
        f.write(contents)

    with open(sourcedir / 'tools.py') as f:
        contents = f.read()
        contents = contents.replace('from .cigsegy', 'from .cigsegyc')
        contents = contents.replace('create_by_sharing_header,',
                                    'create_by_sharing_header)')
        contents = contents.replace(
            '_load_prestack3D, kBinaryHeaderHelp, kTraceHeaderHelp)', '')
        contents = contents.replace('from . import utils', '')

    with open(targetdir / 'tools.py', 'w') as f:
        f.write(contents)

    with open(sourcedir / 'utils.py') as f:
        contents = f.read()
        contents = contents.replace('from .cigsegy', 'from .cigsegyc')
        contents = contents.replace(', kTraceHeaderHelp', '')

    with open(targetdir / 'utils.py', 'w') as f:
        f.write(contents)

    with open(sourcedir / 'cigsegy.pyi') as f:
        contents = f.read()
        contents = contents.replace('import cigsegy', '')

    with open(targetdir / 'cigsegyc.py', 'w') as f:
        f.write(contents)

    with open(sourcedir / 'plot.py') as f:
        contents = f.read()
        contents = contents.replace('from .cigsegy', 'from .cigsegyc')

    with open(targetdir / 'plot.py', 'w') as f:
        f.write(contents)


def clean_up(app, exception):  # noqa: ARG001
    os.remove(DIR / 'readme.rst')
    shutil.rmtree(DIR / 'cigsegy')


def setup(app):
    # Add hook for building doxygen xml when needed
    app.connect("builder-inited", generate_doxygen_xml)

    # Copy the readme in
    app.connect("builder-inited", prepare)

    # Clean up the generated readme
    app.connect("build-finished", clean_up)
