# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import warnings
import subprocess
import sys
import shutil
from pathlib import Path
from sphinx.deprecation import RemovedInSphinx80Warning

DIR = Path(__file__).parent.resolve()

warnings.filterwarnings(
    "ignore",
    message=r"Sphinx 8 will drop support for representing paths as strings.*",
    category=RemovedInSphinx80Warning,
    module=r"breathe\.project",
)
os.environ["CIGSEGY_BUILDING_DOCS"] = "1"

project = 'cigsegy'
copyright = '2024, Jintao Li'
author = 'Jintao Li'
verison_path = DIR.parent / "VERSION.txt"
version =  verison_path.read_text().strip()
language = 'en'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["breathe", "sphinx.ext.autodoc", "sphinx.ext.napoleon"]
autodoc_typehints = "none"

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


def update_doxyfile():
    doxyfile_in = DIR / 'Doxyfile.in'
    doxyfile = DIR / 'Doxyfile'

    with open(doxyfile_in, 'r', encoding='utf-8') as infile:
        content = infile.read()

    print(f'Overide version: {version}')
    content = content.replace('@VERSION@', version)

    with open(doxyfile, 'w', encoding='utf-8') as outfile:
        outfile.write(content)


def generate_doxygen_xml(app):
    update_doxyfile()
    build_dir = os.path.join(app.confdir, "_doxygenxml")
    if not os.path.exists(build_dir):
        os.mkdir(build_dir)

    if shutil.which("doxygen") is None:
        print("doxygen not found; skip C++ XML generation.")
        return

    try:
        subprocess.call(["doxygen", "--version"])
        retcode = subprocess.call(["doxygen"], cwd=app.confdir)
        if retcode < 0:
            sys.stderr.write(f"doxygen error code: {-retcode}\n")
    except OSError as e:
        sys.stderr.write(f"doxygen execution failed: {e}\n")


def clean_up(app, exception):  # noqa: ARG001
    for path in (DIR / 'Doxyfile',):
        if path.exists():
            os.remove(path)
    # shutil.rmtree(DIR / 'cigsegy')


def setup(app):
    # Add hook for building doxygen xml when needed
    app.connect("builder-inited", generate_doxygen_xml)

    # Clean up generated files.
    app.connect("build-finished", clean_up)
