# Activate pyNR's virtual environment (made by env/setup.sh).   source env/env.sh
# Safe to source under "set -eu". VENV (default <repo>/.venv) chooses another venv.
_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)
_VENV=${VENV:-$_ROOT/.venv}
if [ -f "$_VENV/bin/activate" ]; then
  set +u; . "$_VENV/bin/activate"; set -u 2>/dev/null || true
  # Numba's compiled-kernel cache lives next to the sources (__pycache__); first runs compile, later runs reuse it
else
  echo "env.sh: no virtual environment at $_VENV; run env/setup.sh first" >&2
fi
unset _ROOT _VENV
