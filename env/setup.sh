#!/usr/bin/env bash
# From spack to a working pyNR: virtual environment from spack's Python, pyNR installed in it, kernels compiled.
#
#   env/setup.sh [extras]          extras of pyproject.toml (default: dev = viz, notebook, docs, test)
#   environment: PYTHON (an interpreter to use instead of spack's), VENV (default <repo>/.venv), SKIP_TESTS=1
#
# Python is taken, in this order, from $PYTHON; the spack environment created in env/ (env/.spack-env/view); the
# named spack environment "pynr"; python3 on the PATH (any Python >= 3.10 works; spack only pins the version).
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
EXTRAS=${1:-dev}
VENV=${VENV:-$ROOT/.venv}
SROOT=${SPACK_ROOT:-$(command -v spack >/dev/null 2>&1 && dirname "$(dirname "$(command -v spack)")" || echo "$HOME/sw/spack")}
for p in "${PYTHON:-}" "$ROOT/env/.spack-env/view/bin/python3" "$SROOT/var/spack/environments/pynr/.spack-env/view/bin/python3" \
         "$(command -v python3 || true)"; do
  [ -n "$p" ] && [ -x "$p" ] && { PY=$p; break; }
done
"$PY" -c 'import sys; assert sys.version_info >= (3, 10), sys.version' || { echo "setup.sh: $PY is older than 3.10" >&2; exit 1; }
echo "== Python: $PY ($("$PY" --version))"
[ -x "$VENV/bin/python" ] || "$PY" -m venv "$VENV"
"$VENV/bin/python" -m pip install --quiet --upgrade pip
echo "== pip install -e \".[$EXTRAS]\" into $VENV"
"$VENV/bin/python" -m pip install --quiet -e "$ROOT[$EXTRAS]"
"$VENV/bin/pynr" thorns > /dev/null && echo "== pynr command works"
if [ "${SKIP_TESTS:-0}" != 1 ] && "$VENV/bin/python" -c "import pytest" 2>/dev/null; then
  echo "== tests (the first run compiles the Numba kernels, ~20 s)"
  (cd "$ROOT" && "$VENV/bin/python" -m pytest -q -x)
fi
echo "Done. In every new shell:  source env/env.sh"
