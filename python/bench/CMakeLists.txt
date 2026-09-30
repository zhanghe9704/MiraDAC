# python/bench/CMakeLists.txt — C++ baseline for the Python benchmark (DA_BUILD_BENCH).

add_executable(bench_cpp bench_cpp.cpp)
target_compile_options(bench_cpp PRIVATE -Wall -Wextra)
target_link_libraries(bench_cpp PRIVATE daStatic)
