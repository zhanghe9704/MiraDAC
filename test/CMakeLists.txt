# test/CMakeLists.txt

set(TEST_SOURCES
    catch_main.cc
    test_smoke.cc
    test_layout.cc
    test_pool.cc
    test_env.cc
    test_numeric.cc
    test_symbolic.cc
    test_multienv.cc
    test_symbolic_cd.cc
)

add_executable(run_tests ${TEST_SOURCES})

target_include_directories(run_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
)

target_compile_options(run_tests PRIVATE -Wall -Wextra)

target_link_libraries(run_tests PRIVATE daStatic)
if(WITH_SYMBOLIC)
    target_include_directories(run_tests PRIVATE ${DA_SYMENGINE_INCLUDES})
    target_link_libraries(run_tests PRIVATE ${DA_SYMENGINE_LIBS})
endif()

add_test(NAME run_tests COMMAND run_tests)

# ASan build (optional target) - uses same TEST_SOURCES (includes test_multienv.cc)
add_executable(run_tests_asan ${TEST_SOURCES})
target_include_directories(run_tests_asan PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
)
target_compile_options(run_tests_asan PRIVATE -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer)
target_link_options(run_tests_asan PRIVATE -fsanitize=address)
target_link_libraries(run_tests_asan PRIVATE daStatic)
if(WITH_SYMBOLIC)
    target_include_directories(run_tests_asan PRIVATE ${DA_SYMENGINE_INCLUDES})
    target_link_libraries(run_tests_asan PRIVATE ${DA_SYMENGINE_LIBS})
endif()

# test_kernels: bit-exact comparison against reference engine.
# Requires the ref lib at /tmp/libref_tpsa.a. Build it first:
#   g++ -std=c++17 -O2 -c -I<repo>/ref/tpsa/include \
#       <repo>/ref/tpsa/src/tpsa_extend.cc -o /tmp/libref_build/tpsa_extend.o
#   ar rcs /tmp/libref_tpsa.a /tmp/libref_build/tpsa_extend.o
set(REF_LIB "/tmp/libref_tpsa.a")
if(EXISTS ${REF_LIB})
    add_executable(test_kernels
        catch_main.cc
        test_kernels.cc
    )
    target_include_directories(test_kernels PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/ref/tpsa/include
    )
    target_compile_options(test_kernels PRIVATE -Wall -Wextra)
    target_link_libraries(test_kernels PRIVATE daStatic ${REF_LIB})
    add_test(NAME test_kernels COMMAND test_kernels)
endif()
