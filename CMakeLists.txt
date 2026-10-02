cmake_minimum_required(VERSION 3.16)
project(da VERSION 1.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Options ---
option(WITH_SYMBOLIC "Build symbolic (SDA) support" ON)
option(DA_IGNORE_SYMENGINE_PIN
    "Accept any SymEngine and only warn when it is not the pinned build (benchmark comparisons only)" OFF)
option(DA_BUILD_TESTS "Build C++ tests" ON)
option(DA_BUILD_EXAMPLES "Build examples" ON)
option(DA_BUILD_BENCH "Build the C++ baseline benchmark for the Python binding" OFF)
option(DA_BUILD_CAPI "Build the C API library miradac_c" OFF)
set(DA_CHECK_ENV 1 CACHE STRING "Enable env guard (1=on, 0=off)")

# A Python build (scikit-build-core) links daStatic into a shared module, so it
# needs PIC, and it builds nothing but the library and the module.
if(SKBUILD)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    set(DA_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(DA_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
endif()

# The C API links daStatic into a shared library, so it needs PIC too.
if(DA_BUILD_CAPI)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# NOTE: DA_CHECK_ENV and DA_WITH_SYMBOLIC are applied below as PUBLIC compile
# definitions on the library targets, NOT with add_compile_definitions().
# add_compile_definitions() is build-tree only: it is not exported as a usage
# requirement, so a downstream find_package(da) consumer compiled the headers
# with the env guard silently disabled and — for a symbolic install — with a
# DAEnv missing pool_e_ (256 bytes instead of 352), an ODR/ABI mismatch.

# --- Symbolic support ---
if(WITH_SYMBOLIC)
    # Finds SymEngine at the version in cmake/symengine_pin.txt and checks the
    # build stamp written by scripts/setup_symengine.sh.
    include(${CMAKE_SOURCE_DIR}/cmake/CheckSymEnginePin.cmake)
    set(DA_SYMENGINE_INCLUDES ${SYMENGINE_INCLUDE_DIRS})
    set(DA_SYMENGINE_LIBS ${SYMENGINE_LIBRARIES})
else()
    set(DA_SYMENGINE_INCLUDES "")
    set(DA_SYMENGINE_LIBS "")
endif()

# --- Library sources ---
set(DA_SOURCES
    src/layout.cpp
    src/env.cpp
    src/engine.cpp
    src/davector.cpp
    src/functions.cpp
    src/base.cpp
    src/symbolic.cpp
)

# --- Shared library ---
add_library(daShared SHARED ${DA_SOURCES})
target_include_directories(daShared
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        ${DA_SYMENGINE_INCLUDES}
)
if(WITH_SYMBOLIC)
    target_link_libraries(daShared PUBLIC ${DA_SYMENGINE_LIBS})
endif()
set_target_properties(daShared PROPERTIES
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    # Keep the directories of linked libraries outside the system paths (the
    # pinned SymEngine prefix) in the installed RUNPATH, so the installed
    # library finds libsymengine without LD_LIBRARY_PATH.
    INSTALL_RPATH_USE_LINK_PATH ON
)
target_compile_definitions(daShared PUBLIC DA_CHECK_ENV=${DA_CHECK_ENV})
if(WITH_SYMBOLIC)
    target_compile_definitions(daShared PUBLIC DA_WITH_SYMBOLIC)
endif()

# --- Static library ---
add_library(daStatic STATIC ${DA_SOURCES})
target_include_directories(daStatic
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        ${DA_SYMENGINE_INCLUDES}
)
if(WITH_SYMBOLIC)
    target_link_libraries(daStatic PUBLIC ${DA_SYMENGINE_LIBS})
endif()
target_compile_definitions(daStatic PUBLIC DA_CHECK_ENV=${DA_CHECK_ENV})
if(WITH_SYMBOLIC)
    target_compile_definitions(daStatic PUBLIC DA_WITH_SYMBOLIC)
endif()

# --- Install rules ---
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS daShared daStatic
    EXPORT daTargets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/da/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/da
)

install(EXPORT daTargets
    FILE      daTargets.cmake
    NAMESPACE da::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/da
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/daConfigVersion.cmake"
    VERSION       ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/daConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/daConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/da
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/daConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/daConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/da
)

# --- Tests ---
if(DA_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# --- Examples ---
if(DA_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# --- Benchmark (C++ baseline for the Python binding) ---
if(DA_BUILD_BENCH)
    add_subdirectory(python/bench)
endif()

# --- C API (for the Julia binding) ---
if(DA_BUILD_CAPI)
    add_subdirectory(capi)
endif()

# --- Python module ---
if(SKBUILD)
    add_subdirectory(python)
endif()
