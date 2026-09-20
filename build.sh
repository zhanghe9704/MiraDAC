#!/usr/bin/env bash
#
# build.sh — configure, build, and test MiraDAC (the `da` library).
#
# The library has two build variants:
#   * numerical only  (WITH_SYMBOLIC=OFF) — no external dependencies.
#   * numerical + symbolic (WITH_SYMBOLIC=ON) — requires SymEngine (and its
#     dependency, GMP). Symbolic support adds da::SDA and NDA<->SDA interop.
#
# Usage:
#   ./build.sh                       # numerical only
#   ./build.sh --symbolic            # numerical + symbolic (SymEngine auto-found)
#   ./build.sh --symbolic \
#       --symengine-dir /path/to/lib/cmake/symengine \
#       --gmp-dir /path/to/gmp/prefix     # point at a non-standard SymEngine/GMP
#   ./build.sh --no-check-env        # compile out the cross-env guard (DA_CHECK_ENV=0)
#   ./build.sh --build-dir build --jobs 8 --no-tests
#
# Environment overrides (alternative to flags):
#   SYMENGINE_DIR, GMP_PREFIX, BUILD_DIR, JOBS
#
set -euo pipefail

SYMBOLIC=OFF
CHECK_ENV=1
BUILD_DIR="${BUILD_DIR:-build}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
RUN_TESTS=1
SYMENGINE_DIR="${SYMENGINE_DIR:-}"
GMP_PREFIX="${GMP_PREFIX:-}"
BUILD_TYPE="${BUILD_TYPE:-Release}"

while [ $# -gt 0 ]; do
  case "$1" in
    --symbolic)        SYMBOLIC=ON ;;
    --no-check-env)    CHECK_ENV=0 ;;
    --build-dir)       BUILD_DIR="$2"; shift ;;
    --symengine-dir)   SYMENGINE_DIR="$2"; shift ;;
    --gmp-dir)         GMP_PREFIX="$2"; shift ;;
    --jobs)            JOBS="$2"; shift ;;
    --no-tests)        RUN_TESTS=0 ;;
    -h|--help)         sed -n '2,30p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

ROOT="$(cd "$(dirname "$0")" && pwd)"

CMAKE_ARGS=(
  -S "$ROOT"
  -B "$BUILD_DIR"
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
  -DWITH_SYMBOLIC="$SYMBOLIC"
  -DDA_CHECK_ENV="$CHECK_ENV"
)

# CMake >= 4 needs a policy floor for older transitive projects (e.g. SymEngine).
CMAKE_ARGS+=( -DCMAKE_POLICY_VERSION_MINIMUM=3.5 )

if [ "$SYMBOLIC" = "ON" ]; then
  [ -n "$SYMENGINE_DIR" ] && CMAKE_ARGS+=( -DSymEngine_DIR="$SYMENGINE_DIR" )
  if [ -n "$GMP_PREFIX" ]; then
    # Support both plain and multiarch include/lib layouts.
    for inc in "$GMP_PREFIX/include" "$GMP_PREFIX/include/x86_64-linux-gnu"; do
      [ -f "$inc/gmp.h" ] && CMAKE_ARGS+=( -DCMAKE_INCLUDE_PATH="$inc" )
    done
    for lib in "$GMP_PREFIX/lib" "$GMP_PREFIX/lib/x86_64-linux-gnu"; do
      [ -e "$lib/libgmp.so" ] && CMAKE_ARGS+=( -DCMAKE_LIBRARY_PATH="$lib" )
    done
  fi
fi

echo "==> cmake ${CMAKE_ARGS[*]}"
cmake "${CMAKE_ARGS[@]}"
echo "==> cmake --build $BUILD_DIR -j$JOBS"
cmake --build "$BUILD_DIR" -j"$JOBS"

if [ "$RUN_TESTS" = "1" ]; then
  # Tests read reference data files from the source test/ directory, so run there.
  echo "==> running tests (from $ROOT/test)"
  RUNTIME_LIBS=""
  if [ "$SYMBOLIC" = "ON" ]; then
    # Make SymEngine/GMP shared libs findable at runtime if they are non-standard.
    [ -n "$SYMENGINE_DIR" ] && RUNTIME_LIBS="$(cd "$SYMENGINE_DIR/../.." && pwd)/lib:$RUNTIME_LIBS"
    [ -n "$GMP_PREFIX" ] && RUNTIME_LIBS="$GMP_PREFIX/lib:$GMP_PREFIX/lib/x86_64-linux-gnu:$RUNTIME_LIBS"
  fi
  ( cd "$ROOT/test" && LD_LIBRARY_PATH="$RUNTIME_LIBS${LD_LIBRARY_PATH:-}" "$ROOT/$BUILD_DIR/test/run_tests" )
fi

echo "==> done."
