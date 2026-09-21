# Standalone consumer smoke test.
#
# NOT part of the main build. It is configured separately against an install
# tree to prove that find_package(da) exports the compile definitions the
# headers need:
#
#   * DA_CHECK_ENV must arrive as 1, or the cross-environment guard is
#     silently compiled out of consumer code.
#   * DA_WITH_SYMBOLIC must arrive for a symbolic install, or DAEnv is
#     missing pool_e_ in the consumer's view (256 bytes vs 352) — an
#     ODR/ABI mismatch against the library.
#
# Usage:
#   cmake -S test/consumer -B <bld> -Dda_DIR=<prefix>/lib/cmake/da
#   cmake --build <bld> && <bld>/consumer_smoke

cmake_minimum_required(VERSION 3.16)
project(da_consumer_smoke CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(da REQUIRED)

add_executable(consumer_smoke main.cc)
target_link_libraries(consumer_smoke PRIVATE da::daStatic)

enable_testing()
add_test(NAME consumer_smoke COMMAND consumer_smoke)
