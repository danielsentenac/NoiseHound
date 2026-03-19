#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${1:-$HOME/NOISEHOUND}"
VENV_DIR="${VENV_DIR:-$HOME/.venvs/noisehound}"

python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

python -m pip install --upgrade pip setuptools wheel
python -m pip install -e "$WORKDIR"

python - <<'PY'
import importlib

try:
    import gwpy.io.gwf  # noqa: F401
    for module_name in ("LDAStools.frameCPP", "framel", "lalframe"):
        try:
            importlib.import_module(module_name)
            print(f"Verified GWF backend: {module_name}")
            break
        except Exception:
            continue
    else:
        raise RuntimeError("no supported gwf backend import found")
except Exception as error:
    print("warning: GWF backend verification failed:", error)
    print("The environment installed, but reading .gwf files may still fail.")
PY

printf 'Environment ready in %s\n' "$VENV_DIR"
