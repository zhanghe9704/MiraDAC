#!/usr/bin/env bash
#
# setup_symengine.sh — build the pinned SymEngine and symengine.py against it.
#
# The only supported way to get a SymEngine that MiraDAC and symengine.py can
# share. All versions and options come from cmake/symengine_pin.txt.
#
# Usage:
#   scripts/setup_symengine.sh [--prefix DIR] [--python PYTHON] [--print-env]
#
#   --prefix DIR     SymEngine install prefix
#                    (default: $HOME/.local/opt/symengine-<version>-<commit[:8]>)
#   --python PYTHON  Python that gets symengine.py (default: <repo>/.venv/bin/python)
#   --print-env      print "export SymEngine_DIR=..." on stdout, for
#                    eval "$(scripts/setup_symengine.sh --print-env)"
#
# Steps:
#   1. read the pin
#   2. skip the SymEngine build if the prefix stamp matches the pin and the
#      installed library's sha256
#   3. download, verify, build and install SymEngine at the pinned commit
#   4. write the stamp <prefix>/share/symengine/miradac-pin.txt
#   5. download and verify the symengine.py sdist
#   6. build and install symengine.py from it against the prefix
#   7. verify the installed symengine.py links the prefix's libsymengine
#   8. with --print-env, print the environment
# Steps 5-6 are skipped when the installed symengine.py already passes step 7.
#
# All progress goes to stderr; stdout carries only the --print-env output.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PIN_FILE="$ROOT/cmake/symengine_pin.txt"
PREFIX=""
PY="$ROOT/.venv/bin/python"
PRINT_ENV=0

while [ $# -gt 0 ]; do
  case "$1" in
    --prefix)    PREFIX="$2"; shift ;;
    --python)    PY="$2"; shift ;;
    --print-env) PRINT_ENV=1 ;;
    -h|--help)   sed -n '2,28p' "$0"; exit 0 ;;
    *) echo "setup_symengine: unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

log() { echo "==> $*" >&2; }
die() { echo "setup_symengine: error: $*" >&2; exit 1; }

# --- 1. Read the pin ---------------------------------------------------------
[ -f "$PIN_FILE" ] || die "pin file not found: $PIN_FILE"
declare -A PIN
OPTIONS=()
while IFS='=' read -r key value; do
  case "$key" in ''|\#*) continue ;; esac
  PIN[$key]="$value"
  case "$key" in SYMENGINE_*) ;; *) OPTIONS+=("$key") ;; esac
done < "$PIN_FILE"
for key in SYMENGINE_COMMIT SYMENGINE_VERSION SYMENGINE_TARBALL_SHA256 \
           SYMENGINE_PY_VERSION SYMENGINE_PY_SDIST_SHA256; do
  [ -n "${PIN[$key]:-}" ] || die "$PIN_FILE has no $key"
done
COMMIT="${PIN[SYMENGINE_COMMIT]}"
VERSION="${PIN[SYMENGINE_VERSION]}"
PY_VERSION="${PIN[SYMENGINE_PY_VERSION]}"
SOVERSION="${VERSION%.*}"

PREFIX="${PREFIX:-$HOME/.local/opt/symengine-$VERSION-${COMMIT:0:8}}"
mkdir -p "$PREFIX"
PREFIX="$(cd "$PREFIX" && pwd)"
STAMP="$PREFIX/share/symengine/miradac-pin.txt"
LIB="$PREFIX/lib/libsymengine.so.$VERSION"
[ -x "$PY" ] || die "Python not found: $PY (create it with: uv venv .venv --python 3.13)"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/setup_symengine.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

sha256() { sha256sum "$1" | cut -d' ' -f1; }
fetch() {  # fetch URL FILE SHA256
  log "downloading $1"
  curl -fsSL -o "$2" "$1" || die "download failed: $1"
  local got; got="$(sha256 "$2")"
  [ "$got" = "$3" ] || die "sha256 mismatch for $1: expected $3, got $got"
}

# --- 2. Is the installed SymEngine already the pinned one? ------------------
stamp_ok() {
  [ -f "$STAMP" ] && [ -f "$LIB" ] || return 1
  local key
  for key in SYMENGINE_COMMIT SYMENGINE_VERSION "${OPTIONS[@]}"; do
    grep -qxF "$key=${PIN[$key]}" "$STAMP" || return 1
  done
  grep -qxF "LIB_SHA256=$(sha256 "$LIB")" "$STAMP"
}

if stamp_ok; then
  log "SymEngine $VERSION ($COMMIT) already installed in $PREFIX, skipping build"
else
  # --- 3. Build and install SymEngine ---------------------------------------
  fetch "https://github.com/symengine/symengine/archive/$COMMIT.tar.gz" \
        "$WORK/symengine.tar.gz" "${PIN[SYMENGINE_TARBALL_SHA256]}"
  tar -xzf "$WORK/symengine.tar.gz" -C "$WORK"
  CMAKE_ARGS=(-S "$WORK/symengine-$COMMIT" -B "$WORK/build-symengine"
              -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_INSTALL_LIBDIR=lib
              -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF)
  for key in "${OPTIONS[@]}"; do CMAKE_ARGS+=("-D$key=${PIN[$key]}"); done
  command -v ninja >/dev/null && CMAKE_ARGS+=(-G Ninja)
  log "configuring SymEngine: ${CMAKE_ARGS[*]}"
  cmake "${CMAKE_ARGS[@]}" >&2 || die "SymEngine configure failed"
  cmake --build "$WORK/build-symengine" -j"$(nproc 2>/dev/null || echo 4)" >&2 \
    || die "SymEngine build failed"
  cmake --install "$WORK/build-symengine" >&2 || die "SymEngine install failed"
  [ -f "$LIB" ] || die "expected library not installed: $LIB"

  # --- 4. Stamp ---------------------------------------------------------------
  mkdir -p "$(dirname "$STAMP")"
  {
    echo "SYMENGINE_COMMIT=$COMMIT"
    echo "SYMENGINE_VERSION=$VERSION"
    for key in "${OPTIONS[@]}"; do echo "$key=${PIN[$key]}"; done
    echo "LIB_SHA256=$(sha256 "$LIB")"
  } > "$STAMP"
  stamp_ok || die "stamp verification failed: $STAMP"
  log "installed SymEngine $VERSION ($COMMIT) into $PREFIX"
fi

# --- 7. Verification of the installed symengine.py --------------------------
py_ok() {
  local wrapper
  wrapper="$("$PY" -c 'import symengine.lib.symengine_wrapper as w; print(w.__file__)' 2>/dev/null)" \
    || { echo "symengine.py is not importable" >&2; return 1; }
  readelf -d "$wrapper" | grep -q "NEEDED.*\[libsymengine.so.$SOVERSION\]" \
    || { echo "$wrapper does not need libsymengine.so.$SOVERSION" >&2; return 1; }
  readelf -d "$wrapper" | grep -E 'R(UN)?PATH' | grep -qF "$PREFIX/lib" \
    || { echo "$wrapper has no RUNPATH into $PREFIX/lib" >&2; return 1; }
  "$PY" -c "import symengine; assert symengine.__version__ == '$PY_VERSION', symengine.__version__" \
    || { echo "symengine.__version__ is not $PY_VERSION" >&2; return 1; }
}

if py_ok 2>/dev/null; then
  log "symengine.py $PY_VERSION already built against $PREFIX, skipping install"
else
  # --- 5. symengine.py sdist --------------------------------------------------
  SDIST="$WORK/symengine-$PY_VERSION.tar.gz"
  fetch "https://files.pythonhosted.org/packages/source/s/symengine/symengine-$PY_VERSION.tar.gz" \
        "$SDIST" "${PIN[SYMENGINE_PY_SDIST_SHA256]}"
  SDIST_COMMIT="$(tar -xzOf "$SDIST" "symengine-$PY_VERSION/symengine_version.txt" | tr -d '[:space:]')"
  [ "$SDIST_COMMIT" = "$COMMIT" ] \
    || die "symengine.py $PY_VERSION pins SymEngine $SDIST_COMMIT, not $COMMIT"

  # --- 6. Build symengine.py against the prefix --------------------------------
  "$PY" -m pip --version >/dev/null 2>&1 || "$PY" -m ensurepip >&2 \
    || die "cannot bootstrap pip in $PY"
  "$PY" -m pip install cython setuptools >&2 || die "cannot install cython/setuptools"
  log "building symengine.py $PY_VERSION against $PREFIX"
  # symengine.py's CMake finds cython on PATH, so put the target Python's bin first.
  PATH="$(dirname "$PY"):$PATH" CMAKE_PREFIX_PATH="$PREFIX" "$PY" -m pip install --no-build-isolation --no-binary symengine \
    --no-cache-dir --force-reinstall --no-deps "$SDIST" >&2 \
    || die "symengine.py build failed"
  py_ok || die "symengine.py verification failed (see above)"
  log "installed symengine.py $PY_VERSION into $PY"
fi

# --- 8. Environment -----------------------------------------------------------
if [ "$PRINT_ENV" = 1 ]; then
  echo "export SymEngine_DIR=$PREFIX/lib/cmake/symengine"
fi
