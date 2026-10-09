# pyNR toolchain with spack (any machine)

pyNR is Python with Numba: "compiling" means Numba compiling the kernels on their first call. The minimal files to
go from spack to a working, compiled pyNR:

| file | role |
|---|---|
| `spack.yaml` | spack environment with Python 3.11 (built with spack's compiler) |
| `setup.sh` | virtual environment from that Python, `pip install -e ".[dev]"`, `pynr thorns`, tests (Numba compiles) |
| `env.sh` | activates the virtual environment in a new shell |

Generated and git-ignored: `env/spack.lock`, `env/.spack-env/` (spack), `.venv/` (the virtual environment).

```bash
spack -e env install          # once: Python 3.11 in env/.spack-env (or: spack env create pynr env/spack.yaml && spack -e pynr install)
env/setup.sh                  # once: .venv + pyNR + tests (~1-2 min, then ~20 s of Numba compilation)
source env/env.sh             # every shell
pynr run par/gauge_wave.par
```

Without spack, `env/setup.sh` uses `python3` from the PATH (≥ 3.10), or `PYTHON=/path/to/python env/setup.sh`.
Conda users: see the installation page of the documentation.
