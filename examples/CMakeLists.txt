# examples/CMakeLists.txt — Stage 5 numeric examples

add_executable(examples examples.cc)
target_include_directories(examples PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(examples PRIVATE -Wall -Wextra)
target_link_libraries(examples PRIVATE daStatic)

add_executable(example_complex_da example_complex_da.cc)
target_include_directories(example_complex_da PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_compile_options(example_complex_da PRIVATE -Wall -Wextra)
target_link_libraries(example_complex_da PRIVATE daStatic)

# Symbolic examples - only build when WITH_SYMBOLIC is ON
if(WITH_SYMBOLIC)
    foreach(_example example_1_symbolic example_2_symbolic example_3_symbolic)
        add_executable(${_example} ${_example}.cc)
        target_include_directories(${_example} PRIVATE
            ${CMAKE_SOURCE_DIR}/include
            ${DA_SYMENGINE_INCLUDES})
        target_compile_options(${_example} PRIVATE -Wall -Wextra)
        target_link_libraries(${_example} PRIVATE daStatic ${DA_SYMENGINE_LIBS})
    endforeach()
endif()

# Stage 9: interop example (WITH_SYMBOLIC only)
if(WITH_SYMBOLIC)
    add_executable(example_interop example_interop.cc)
    target_include_directories(example_interop PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${DA_SYMENGINE_INCLUDES})
    target_compile_options(example_interop PRIVATE -Wall -Wextra)
    target_link_libraries(example_interop PRIVATE daStatic ${DA_SYMENGINE_LIBS})
endif()

# Stage C6: symbolic complex DA example (WITH_SYMBOLIC only)
if(WITH_SYMBOLIC)
    add_executable(example_complex_symbolic example_complex_symbolic.cc)
    target_include_directories(example_complex_symbolic PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${DA_SYMENGINE_INCLUDES})
    target_compile_options(example_complex_symbolic PRIVATE -Wall -Wextra)
    target_link_libraries(example_complex_symbolic PRIVATE daStatic ${DA_SYMENGINE_LIBS})
endif()

# Stage 9: benchmark (numerical, requires ref/tpsa static lib)
set(REF_LIB "/tmp/libref_tpsa.a")
if(EXISTS ${REF_LIB})
    add_executable(benchmark_composition benchmark_composition.cc)
    target_include_directories(benchmark_composition PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ref/tpsa/include)
    target_compile_options(benchmark_composition PRIVATE -Wall -Wextra -O2)
    target_link_libraries(benchmark_composition PRIVATE daStatic ${REF_LIB})
endif()
