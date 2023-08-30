@echo off

set "envs=py37 py38 base py310 py311"

for %%e in (%envs%) do (
    echo Executing in %%e environment
    conda activate %%e
    python -m pip install --upgrade pip
    pip install pybind11
    pip wheel . --build-option="--fmt_root=D:\fmt"
    conda deactivate
)